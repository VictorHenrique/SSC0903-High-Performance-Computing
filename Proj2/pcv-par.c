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

// Permutacao na parte interna
void heap_permutation(int *base, int size, int *perm, int n, int *idx) {
    if (n == 1) {
        // SIMD
        #pragma omp simd 
        for (int i = 0; i <= size; i++) 
            perm[(*idx) * (size+1) + i] = base[i];
        (*idx)++;

        return;
    }

    // SIMD
    for (int i = 1; i < n; i++) {
        heap_permutation(base, size, perm, n - 1, idx);
        if (n % 2) 
            swap(base, 1, n - 1);
        else    
            swap(base, i, n - 1);
    }
}   

void generate_permutations(int *city, int n, int src, int *permutations) {
    // SIMD
    int *base = malloc((n + 1) * sizeof(int));
    base[0] = base[n] = src;
    
    for (int i = 0, j = 1; i < n; i++) 
        if (city[i] != src) base[j++] = city[i]; 

    int idx = 0;
    heap_permutation(base, n, permutations, n, &idx);

    free(base);
}

int *read_adj(int n) {
    int n_sqr = n*n, *adj_matrix = malloc(n_sqr * sizeof(int));
    for (int i = 0; i < n_sqr; i++) 
        scanf("%d", &adj_matrix[i]);
    return adj_matrix;
}

int *create_cities(int n) {
    int *city = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) 
        city[i] = i;

    return city;
}

int get_cost(int *path, int *adj_matrix, int n, int rank) {
    int cost = 0;
    // #pragma omp simd reduction(+: cost)
    for (int i = 1; i <= n; i++) {
        int row = path[i - 1], col = path[i];
        cost += adj_matrix[row * n + col];  // equivalente a adj_matrix[i][j] 
    }
    
    return cost;
}

int main(int argc, char *argv[]) {

    /* MPI */
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    BEST_PATH optimal = {NULL, INT_MAX};

    double start, end;

    int rank, nprocs, n, n_paths, max, tag = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int *adj_matrix, *path;
    if (rank == MAIN) {
        int src;
        scanf("%d", &n);
        MPI_Bcast(&n, 1, MPI_INT, MAIN, MPI_COMM_WORLD);

        max = factorial(n-1);
        int arr_length = max * (n+1); 
        path = malloc(arr_length * sizeof(int));
        
        int *city = create_cities(n);
        adj_matrix = read_adj(n);
        
        scanf("%d", &src);

        start = MPI_Wtime();

        generate_permutations(city, n, src, path);

        /* Divisão da carga de trabalho */
        n_paths = max / nprocs;
        int rem = max % nprocs;
        int cur_pos = rem ? n_paths + 1: n_paths;
        cur_pos *= n + 1;
        for (int i = 1; i < nprocs; i++) {
            int node_workload = i < rem ? n_paths + 1 : n_paths;
            node_workload *= n+1;
            
            MPI_Isend(path + cur_pos, node_workload, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
            MPI_Isend(adj_matrix, n*n, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &request);
            
            cur_pos += node_workload;
        }   

        free(city);
    } else {
        MPI_Bcast(&n, 1, MPI_INT, MAIN, MPI_COMM_WORLD);
        int max = factorial(n-1), rem = max % nprocs;

        n_paths = max / nprocs;
        if (rank < rem) 
            n_paths++;

        int arr_length = n_paths * (n + 1);
        path = malloc(arr_length * sizeof(int));

        MPI_Recv(path, arr_length, MPI_INT, MAIN, tag, MPI_COMM_WORLD, &status);

        adj_matrix = malloc(n*n * sizeof(int));
        MPI_Recv(adj_matrix, n*n, MPI_INT, MAIN, tag + 1, MPI_COMM_WORLD, &status);
    }

    // #pragma omp declare reduction(best_path : struct _best_path : omp_out = omp_in.cost > omp_out.cost ? omp_in : omp_out)
    // #pragma omp simd reduction(best_path: optimal)
    for (int i = 0; i < n_paths; i++) {
        int path_start = i * (n+1);
        int cost = get_cost(&path[path_start], adj_matrix, n, rank);
        if (optimal.cost > cost) {
            optimal.cost = cost;
            optimal.path = &path[path_start]; 
        }
    }

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
                optimal.path = proc_path[i];
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

    free(path), free(adj_matrix);
    MPI_Finalize();

    return 0;
}