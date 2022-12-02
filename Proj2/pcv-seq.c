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

void heap_permutation(int *base, int size, int **perm, int n, int *idx) {
    if (n == 1) {
        for (int i = 0; i <= size; i++) 
            perm[*idx][i] = base[i];
        (*idx)++;
        return;
    }

    for (int i = 1; i < n; i++) {
        heap_permutation(base, size, perm, n - 1, idx);
        if (n % 2) 
            swap(base, 1, n - 1);
        else    
            swap(base, i, n - 1);
    }
}   

int **generate_permutations(int *city, int n, int src) {
    int max = factorial(n - 1);
    int **permutations = malloc(max * sizeof(int *));

    for (int i = 0; i < max; i++) 
        permutations[i] = malloc((n + 1) * sizeof(int));

    int *base = malloc((n + 1) * sizeof(int));
    base[0] = base[n] = src;

    for (int i = 0, j = 1; i < n; i++) 
        if (city[i] != src) base[j++] = city[i]; 

    int idx = 0;
    heap_permutation(base, n, permutations, n, &idx);

    free(base);
    
    return permutations;
}

int solve_tsp(int *path, int **adj_matrix, int n) {
    int cost = 0;
    for (int i = 1; i <= n; i++)
        cost += adj_matrix[path[i-1]][path[i]];

    return cost;
}

void free_mem(int **adj_matrix, int *city, int **permutations, int n, int max) {
    for (int i = 0; i < n; i++)
        free(adj_matrix[i]);
    free(adj_matrix);
    free(city);

    for (int i = 0; i < max; i++)
        free(permutations[i]);
    free(permutations);
}

int main() {
    int n;
    scanf("%d", &n);

    int *city = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) 
        city[i] = i;

    int **adj_matrix = malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        adj_matrix[i] = malloc(n * sizeof(int));
        
        adj_matrix[i][i] = 0;
        for (int j = 0; j < n; j++) 
            scanf("%d", &adj_matrix[i][j]);
    }

    int src, max = factorial(n-1); 
    scanf("%d", &src);

    double seq_start = omp_get_wtime();
    int **path = generate_permutations(city, n, src);
    
    BEST_PATH route = {NULL, INT_MAX};
    for (int i = 0; i < max; i++) {
        int cost = solve_tsp(path[i], adj_matrix, n); 
        if (route.cost > cost) {
            route.cost = cost;
            route.path = path[i];
        }
    }
    double seq_end = omp_get_wtime();

    printf("Best path: ");
    for (int i = 0; i <= n; i++) 
        printf("%d ", route.path[i]);
    printf("\nCost: %d\n", route.cost);
    
    printf("Tempo de resposta sem considerar E/S, em segundos: %.3fs\n", seq_end - seq_start);
    
    free_mem(adj_matrix, city, path, n, max);

    return 0;
}
