#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <omp.h>

#define MAIN 0
#define LIST_END -1
#define LIST_START 0

const long long fact[16] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000};

typedef struct _best_path{
    int *path, cost;
} BEST_PATH;

long long get_factorial(int n) {
    long long res = 1;
    for (int i = 2; i <= n; i++)
        res*=i;

    return res;
}

int get_cost(int *path, int *adj_matrix, int n) {
    int cost = adj_matrix[path[0]] + adj_matrix[n * path[n-2]];
    #pragma omp simd reduction(+: cost)
    for (int i = 1; i < n - 1; i++) {
        int row = path[i - 1], col = path[i];
        cost += adj_matrix[row * n + col];  
    }
    
    return cost;
}

void get_kth_perm(int n, long long k, int *route) {
    n--;
    int *city = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) city[i] = i+1;

    int *perm = malloc(n * sizeof(int));
    int perm_i = 0;

    k--;
    for (int i = 1, j = 0; i <= n; i++, j++) {
        long long factorial = n - i> 15 ? get_factorial(n - i) : fact[n - i];
        int idx = k/factorial;
        
        perm[j] = city[idx];
        for (int l = idx; l < n - 1; l++)
            city[l] = city[l+1];
        
        k -= idx*factorial;
    }

	for (int i = 0; i < n; i++)
        route[i] = perm[i];

    free(perm);
    free(city);
}

void swap(int *route, int i, int j) {
    int temp = route[i];
    route[i] = route[j];
    route[j] = temp;
}

void next_permutation(int *route, int n) {
    int i = n - 2, j;
    while (i >= 0 && route[i + 1] <= route[i]) i--;
    
    if (i >= 0) {
        j = n - 1;
        while (route[j] <= route[i]) j--;
        swap(route, i, j); 
    }
    int k = i + 1, l = n - 1;
    for (; k < l; k++, l--) 
        swap(route, k, l);
}

void solve_tsp(long long starting_idx, int size, long long num_of_perms, BEST_PATH *best_path, int *adj_matrix) {
    int *route = malloc(size * sizeof(int));

    get_kth_perm(size, starting_idx, route);
    int cost = get_cost(route, adj_matrix, size);
    if (cost < best_path->cost) {
        best_path->cost = cost;
        best_path->path[0] = best_path->path[size] = 0;
        #pragma omp simd
        for(int i = 0; i < size - 1; i++)
            best_path->path[i+1] = route[i];
    }

    for (long long i = 0; i < num_of_perms - 1; i++) {
        next_permutation(route, size - 1);

        int cost = get_cost(route, adj_matrix, size);
        if (cost < best_path->cost) {
            best_path->cost = cost;
            #pragma omp simd
            for(int i = 1; i < size; i++)
                best_path->path[i] = route[i-1];
        }
    }
    free(route);
}

int main(int argc, char *argv[]) {

    /* MPI */
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    BEST_PATH optimal = {NULL, INT_MAX};

    double start, end;

    int rank, nprocs, n, tag = 0, seed = 12;
    long long max, rem, node_workload;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) return 1;
    else n = atoi(argv[1]);
    srand(seed);

    int *adj_matrix, n_sqrd = n*n;

    adj_matrix = malloc(n_sqrd * sizeof(int));
    for (int i = 0, j = 0; i < n_sqrd; i++) {
        if (j * n + j == i) {           // Elemento da diagonal
            adj_matrix[i] = 0;
            j++; 
        } else {
            adj_matrix[i] = rand() % 100;
            if (adj_matrix[i] == 0) adj_matrix[i] = (n+1) * 101;  // custo "infinito" -> suficiente para descartar o caminho
        }
    }
    
    /* Medindo o tempo */
    if (rank == MAIN) 
        start = MPI_Wtime();  
    
    /* Carga de trabalho do processo principal */
    max = n > 16 ? get_factorial(n - 1) : fact[n - 1];
    rem = max % nprocs;
    node_workload = max / nprocs;
    if (rem) node_workload++;

    /* Gerando as permutacoes */
    long long starting_idx = -1;         
    if (rank < rem)
        starting_idx = rank * node_workload + 1;
    else if (max > nprocs)
        starting_idx = (rem * (node_workload + 1)) + ((rank - rem) * node_workload) + 1;

    /* Melhor caminho e custo retornados em optimal */
    optimal.path = malloc((n+1) * sizeof(int));
    optimal.path[0] = optimal.path[n] = 0;
    if (starting_idx >= 0) 
        solve_tsp(starting_idx, n, node_workload, &optimal, adj_matrix);  

    if (rank == MAIN) {
        int num_of_receives = max < nprocs ? max : nprocs;
        int *cost = malloc(num_of_receives * sizeof(int));
        int **proc_path = malloc(num_of_receives * sizeof(int *));
        for (int i = 1; i < num_of_receives; i++) {
            proc_path[i] = malloc((n+1) * sizeof(int));
            MPI_Recv(&(cost[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(proc_path[i], n+1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
        }

        for (int i = 1; i < num_of_receives; i++) {
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
        for (int i = 1; i < num_of_receives; i++) 
            free(proc_path[i]);
        free(proc_path), free(cost);
    } else if (starting_idx != -1) {
        tag = rank;
        MPI_Send(&(optimal.cost), 1, MPI_INT, MAIN, tag, MPI_COMM_WORLD);
        MPI_Send(optimal.path, n+1, MPI_INT, MAIN, tag, MPI_COMM_WORLD);
    }
    
    free(optimal.path), free(adj_matrix);
    MPI_Finalize();

    return 0;
}
