#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define BUCKET_LEN 101

typedef struct _cities {
    int *grade_frequencies; // per city
    double mean, sd, median;
    int min, max;
} CITIES;

typedef struct _region {
    CITIES *city;
    int *grade_frequencies; // per region
    double mean, sd, median;
    int min, max;
} REGIONS;


/* Geração pseudo-aleatória das notas de cada aluno a partir de uma seed */
void get_grades(int *grades, int R, int C, int A, int seed) {
    srand(seed);
    for (int i = 0, student = 0; i < R; i++) 
        for (int j = 0; j < C; j++) 
            for (int k = 0; k < A; k++) 
                grades[student++] = rand() % BUCKET_LEN;    
}

/* Counting Sort */
void counting_sort(int *grade, int *bucket, int start, int end) {
    for (int i = start; i <= end; i++)
        bucket[grade[i]]++;

    for (int i = 0, j = start; i < BUCKET_LEN; i++) 
        for (int k = bucket[i]; k; k--)  
            grade[j++] = i; 
}

void add_buckets(int *f_bucket, int *s_bucket) {
    for (int i = 0; i < BUCKET_LEN; i++)
        f_bucket[i] += s_bucket[i];
}

void sort_from_bucket(int *grade, int *bucket, int start, int end) {
    for (int i = 0, j = start; i < BUCKET_LEN; i++) 
        for (int k = bucket[i]; k; k--) 
            grade[j++] = i;
}

/* frequencia * nota */
double mean(int *grade_freq, int num_of_elems) {
    double mean = 0;
    for (int i = 0; i < BUCKET_LEN; i++)
        mean += grade_freq[i] * i;
    
    return mean / num_of_elems;
}

/* frequencia * (nota - media)^2 */
double standard_deviation(int *bucket, double mean, int num_of_elems) {
    double sum = 0;
    for (int i = 0; i < BUCKET_LEN; i++)
            sum += bucket[i] * pow(i - mean, 2);

    return sqrt(sum / num_of_elems);
}

double grades_median(int *grades, int first_element, int n_elements){
    double median;
    if (!(n_elements % 2))
        median = (double)(grades[first_element - 1 + (n_elements / 2)] + grades[first_element + n_elements / 2]) / 2;
    else 
        median = grades[first_element + n_elements / 2];

    return median;
}

REGIONS *initialize_regions(int R, int C) {
    REGIONS *br = malloc(sizeof(REGIONS) * R);
    for (int i = 0; i < R; i++) {
        br[i].city = calloc(C, sizeof(CITIES));
        br[i].grade_frequencies = calloc(BUCKET_LEN, sizeof(int));
        for (int j = 0; j < C; j++)
            br[i].city[j].grade_frequencies = calloc(BUCKET_LEN, sizeof(int));
    }        

    return br;
}

/* Prints */
void print_cities(REGIONS *brazil, int R, int C) {
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n",
                i, j, 
                brazil[i].city[j].min,  
                brazil[i].city[j].max, 
                brazil[i].city[j].median, 
                brazil[i].city[j].mean, 
                brazil[i].city[j].sd
            );
        }
        printf("\n");
    }
    printf("\n");
}

void print_regions(REGIONS *brazil, int R) {
    for (int i = 0; i < R; i++) {
        printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", 
            i, brazil[i].min, 
            brazil[i].max, 
            brazil[i].median, 
            brazil[i].mean, 
            brazil[i].sd
        );
    }
    printf("\n");
}

/* Liberando a memória */
void free_memory(REGIONS *arr, int R, int C, int A) {
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            free(arr[i].city[j].grade_frequencies);
        }
        free(arr[i].grade_frequencies);
        free(arr[i].city);
    }            
    free(arr);
}

/* Cada regiao tem R*C elementos (de 0 até R*c -1) */
int main() {
    int R, C, A, seed;
    scanf("%d%d%d%d", &R, &C, &A, &seed); 
    
    int *grades = malloc(sizeof(int) *  R * C * A);
    get_grades(grades, R, C, A, seed);
    double best_city[3] = {0}, best_region[2] = {0};

    /* Prints por cidade de regiao */
    REGIONS *brazil = initialize_regions(R, C);

    double seq_start = omp_get_wtime();
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            // (i * C * A) + (j * A) = primeira pos da regiao + primeira posicao da cidade dentro da regiao
            int region_start = A * (i * C + j); 
            counting_sort(grades, brazil[i].city[j].grade_frequencies, region_start, region_start + A - 1);

            brazil[i].city[j].min = grades[region_start];
            brazil[i].city[j].max = grades[region_start + A - 1];
            brazil[i].city[j].median = grades_median(grades, region_start, A);
            brazil[i].city[j].mean = mean(brazil[i].city[j].grade_frequencies, A); 
            brazil[i].city[j].sd = standard_deviation(brazil[i].city[j].grade_frequencies, brazil[i].city[j].mean, A);

            if (brazil[i].city[j].mean > best_city[2]) {
                best_city[0] = i;
                best_city[1] = j;
                best_city[2] = brazil[i].city[j].mean;
            }
        }
    }

    /* Print por regiao */
    for (int i = 0; i < R; i++) {
        brazil[i].mean = 0;
        for (int j = 0; j < C; j++) {
            add_buckets(brazil[i].grade_frequencies, brazil[i].city[j].grade_frequencies);
            brazil[i].mean += brazil[i].city[j].mean;
        } 
        brazil[i].mean /= C;

        int region_start = i * C * A, region_end = region_start - 1 + (C * A);
        sort_from_bucket(grades, brazil[i].grade_frequencies, region_start, region_end);

        brazil[i].min = grades[region_start];
        brazil[i].max = grades[region_end];
        brazil[i].median = grades_median(grades, region_start, C*A);            
        brazil[i].sd = standard_deviation(brazil[i].grade_frequencies, brazil[i].mean, C*A);

        // printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, min, max, median, brazil[i].mean, sd);

        if (brazil[i].mean > best_region[1]) {
            best_region[0] = i;
            best_region[1] = brazil[i].mean;
        }
    }   

    /* Print do Brasil */
    double mean = 0;
    int bucket[BUCKET_LEN] = {0};
    for (int i = 0; i < R; i++) 
        add_buckets(bucket, brazil[i].grade_frequencies), mean += brazil[i].mean;    
    mean /= R;

    sort_from_bucket(grades, bucket, 0, R*C*A - 1);
    int min = grades[0], max = grades[R*C*A - 1];
    double median = grades_median(grades, 0, R * C * A);
    double sd = standard_deviation(bucket, mean, R*C*A);

    double seq_end = omp_get_wtime();
    
    // print_cities(brazil, R, C);
    // print_regions(brazil, R);

    printf("Brasil: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n\n", min, max, median, mean, sd);

    /* Melhor cidade e melhor regiao */
    printf("Melhor região: Região %.0f\nMelhor cidade: Região %.0f, Cidade %.0f\n", best_region[0], best_city[0], best_city[1]);


    printf("Tempo: %fs\n", seq_end - seq_start);

    free_memory(brazil, R, C, A);
    free(grades);

    return 0;
}
 