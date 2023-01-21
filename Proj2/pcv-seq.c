#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

typedef struct _best_path{
    int *path, cost;
} BEST_PATH;

int factorial(int n) {
    int res = 1;
    for (int i = n; i > 0; i--) res *= i;

    return res;
}

void swap(int *list, int a, int b) {
    int old_a = list[a];
    list[a] = list[b];
    list[b] = old_a;
}

int get_cost(int *path, int **adj_matrix, int n) {
    int cost = 0;
    for (int i = 1; i <= n; i++)
        cost += adj_matrix[path[i-1]][path[i]];

    return cost;
}

void heap_permutation(int *base, int size, int n, BEST_PATH *optimal, int **adj_matrix) {
    if (n == 1) {
        int cost = get_cost(base, adj_matrix, size);
        if (cost < optimal->cost) {
            optimal->cost = cost;
            for (int i = 0; i <= size; i++)
                optimal->path[i] = base[i];
        }
        return;
    }

    for (int i = 1; i < n; i++) {
        heap_permutation(base, size, n - 1, optimal, adj_matrix);
        if (n % 2) 
            swap(base, 1, n - 1);
        else    
            swap(base, i, n - 1);
    }
}   

void solve_tsp(int *city, int **adj_matrix, int n, int src, BEST_PATH *optimal) {
    int *base = malloc((n + 1) * sizeof(int));
    base[0] = base[n] = src;

    for (int i = 0, j = 1; i < n; i++) 
        if (city[i] != src) base[j++] = city[i]; 

    heap_permutation(base, n, n, optimal, adj_matrix);

    free(base);
}

void free_mem(int **adj_matrix, int *city, int n, int max, BEST_PATH optimal) {
    for (int i = 0; i < n; i++)
        free(adj_matrix[i]);
    free(adj_matrix);
    free(city), free(optimal.path);
}

int main(int argc, char *argv[]) {
    int n, seed = 12;
    if (argc > 2) return 1;
    else n = atoi(argv[1]);
    srand(seed);

    int **adj_matrix = malloc(n * sizeof(int *));
    for (int i = 0, j = 0; i < n; i++) {
        adj_matrix[i] = malloc(n * sizeof(int));
        for (int j = 0; j < n; j++) {
            if (j == i) {
                adj_matrix[i][j] = 0;
            } else {
                adj_matrix[i][j] = rand() % 100;
                if (adj_matrix[i][j] == 0) adj_matrix[i][j] = (n+1) * 101;  // custo "infinito" -> suficiente para descartar o caminho
            }
        }
    }
    
    int *city = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) 
        city[i] = i;

    int src = 0, max = factorial(n-1); 
    double seq_start = omp_get_wtime();
    
    BEST_PATH route = {malloc((n+1) * sizeof(int)), INT_MAX};
    solve_tsp(city, adj_matrix, n, src, &route);
    
    double seq_end = omp_get_wtime();

    printf("Best path: ");
    for (int i = 0; i <= n; i++) 
        printf("%d ", route.path[i]);
    printf("\nCost: %d\n", route.cost);
    
    printf("Tempo de resposta sem considerar E/S, em segundos: %.3fs\n", seq_end - seq_start);
    free_mem(adj_matrix, city, n, max, route);

    return 0;
}
