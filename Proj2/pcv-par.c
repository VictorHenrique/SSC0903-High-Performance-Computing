#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <omp.h>

#define MAIN 0

typedef struct _best_path{
    int *path, cost;
} BEST_PATH;

int factorial(int n) {
    int res = 1;
    #pragma omp simd reduction(*: res)
    for (int i = n; i > 0; i--) res *= i;
    return res;
}

void swap(int *list, int a, int b) {
    int old_a = list[a];
    list[a] = list[b];
    list[b] = old_a;
}

int get_cost(int *path, int *adj_matrix, int n, int rank) {
    int cost = 0;
    for (int i = 1; i <= n; i++) {
        int row = path[i - 1], col = path[i];
        cost += adj_matrix[row * n + col];  // equivalente a adj_matrix[i][j] 
    }
    
    return cost;
}

// Permutacao na parte interna
void solve_tsp(int *base, int size, BEST_PATH *best_path, int n, int *adj_matrix, int rank) {
    if (n == 1) {
        int cost = get_cost(base, adj_matrix, size, rank);
        if (best_path->cost > cost) {
            best_path->cost = cost;
            for(int i = 0; i <= size; i++)
                best_path->path[i] = base[i];
        }

        return;
    }

    // SIMD
    for (int i = 1; i < n; i++) {
        solve_tsp(base, size, best_path, n - 1, adj_matrix, rank);
        if (n % 2) 
            swap(base, 1, n - 1);
        else    
            swap(base, i, n - 1);
    }
}   

void read_adj(int *adj_matrix, int n) {
    int n_sqr = n*n;
    for (int i = 0; i < n_sqr; i++) 
        scanf("%d", &adj_matrix[i]);
}

int *create_cities(int n) {
    int *city = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) 
        city[i] = i;

    return city;
}

int *create_bases(int n, int src, int *city) {
    int *base = malloc((n*n - 1) * sizeof(int));
    for (int i = 0, row = 0; i < n; i++) {
        if (city[i] == src) continue;
        
        int path_start = (row++)*(n+1), path_end = path_start + n;
        base[path_start] = base[path_end] = src;
        base[path_end - 1] = city[i];
        for (int j = 0, k = path_start + 1; j < n; j++)
            if (city[j] != src && city[j] != city[i]) base[k++] = city[j]; 
    }

    return base;
}

int main(int argc, char *argv[]) {

    /* MPI */
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    BEST_PATH optimal = {NULL, INT_MAX};

    double start, end;

    int rank, nprocs, n, bases_per_proc, max, tag = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int *adj_matrix, *path, *bases;
    if (rank == MAIN) {
        int src;
        scanf("%d", &n);
        MPI_Bcast(&n, 1, MPI_INT, MAIN, MPI_COMM_WORLD);
        
        int *city = create_cities(n);
        adj_matrix = malloc(n*n * sizeof(int));
        read_adj(adj_matrix, n);

        scanf("%d", &src);

        start = MPI_Wtime();
        bases = create_bases(n, src, city);

        /* Divisão da carga de trabalho */
        bases_per_proc = (n-1) / nprocs;
        int rem = (n-1) % nprocs;
        int cur_pos = rem ? bases_per_proc + 1: bases_per_proc;
        cur_pos *= n + 1;
        for (int i = 1; i < nprocs; i++) {
            int node_workload = i < rem ? bases_per_proc + 1 : bases_per_proc;
            node_workload *= n+1;

            MPI_Isend(bases + cur_pos, node_workload, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
            MPI_Isend(adj_matrix, n*n, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &request);
            
            cur_pos += node_workload;
        }   

        if (rem) bases_per_proc++;

        free(city);
    } else {
        MPI_Bcast(&n, 1, MPI_INT, MAIN, MPI_COMM_WORLD);
        int max = factorial(n-1), rem = (n-1) % nprocs;

        bases_per_proc = (n-1) / nprocs;
        if (rank < rem) 
            bases_per_proc++;

        int bases_length = bases_per_proc * (n + 1);
        bases = malloc(bases_length * sizeof(int));
        MPI_Recv(bases, bases_length, MPI_INT, MAIN, tag, MPI_COMM_WORLD, &status);

        adj_matrix = malloc(n*n * sizeof(int));
        MPI_Recv(adj_matrix, n*n, MPI_INT, MAIN, tag + 1, MPI_COMM_WORLD, &status);
    }
    
    /* Gerando as permutacoes */
    int num_of_perms = factorial(n-2), arr_length = bases_per_proc * num_of_perms * (n+1);
    optimal.path = malloc((n+1) * sizeof(int));

    /* Cada nó recebeu  *bases_per_proc* bases, cada uma gerando (n-2) 
      permutacoes, que são arrays de tamanho n+1  */
    path = malloc(arr_length * sizeof(int));
    for (int i = 0; i < bases_per_proc; i++) 
        solve_tsp(&bases[i*(n+1)], n, &optimal, n-1, adj_matrix, rank);
    
    if (rank == MAIN) {
        int *cost = malloc(nprocs * sizeof(int));
        int **proc_path = malloc(nprocs * sizeof(int *));
        for (int i = 1; i < nprocs; i++) {
            proc_path[i] = malloc((n+1) * sizeof(int));
            MPI_Recv(&(cost[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(proc_path[i], n+1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
        }

        for (int i = 1; i < nprocs; i++) {
            if (cost[i] < optimal.cost) {
                optimal.cost = cost[i];

                for (int j = 0; j <= n; j++)
                    optimal.path[j] = proc_path[i][j];
            }
        }
        
        end = MPI_Wtime();

        printf("Best path: ");
        for (int i = 0; i <= n; i++) 
            printf("%d ", optimal.path[i]);
        printf("\nCost: %d\n", optimal.cost);

        printf("Tempo de resposta sem considerar E/S, em segundos: %.3fs\n", end - start);

        /* Liberando Memória */
        for (int i = 1; i < nprocs; i++) 
            free(proc_path[i]);
        free(proc_path), free(cost);
    } else {
        tag = rank;
        MPI_Send(&(optimal.cost), 1, MPI_INT, MAIN, tag, MPI_COMM_WORLD);
        MPI_Send(optimal.path, n+1, MPI_INT, MAIN, tag, MPI_COMM_WORLD);
    }

    
    free(optimal.path), free(adj_matrix), free(path), free(bases);
    MPI_Finalize();

    return 0;
}