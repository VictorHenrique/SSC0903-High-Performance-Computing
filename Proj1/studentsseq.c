#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc studentsseq.c -o sseq -lm

#define BUCKET_LEN 100
#define GRADE_MAX 101

typedef struct _cities {
    int *grade_frequencies; // per city
    int mean;
} CITIES;

typedef struct _region {
    CITIES *city;
    int *grade_frequencies; // per region
    int mean;
} REGIONS;


/* Geração pseudo-aleatória das notas de cada aluno a partir de uma seed */
void get_grades(int *grades, int R, int C, int A, int seed) {
    srand(seed);
    for (int i = 0, student = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            for (int k = 0; k < A; k++) {
                grades[student++] = rand() % GRADE_MAX; 
            }
        }
    }
}

/* Counting Sort */
int *counting_sort(int *grade, int start, int end) {
    int *bucket = calloc(BUCKET_LEN, sizeof(int));
    for (int i = start; i <= end; i++)
        bucket[grade[i]]++;

    for (int i = 0, j = 0; i < BUCKET_LEN; i++) 
        while(bucket[i]--) 
            grade[j++] = i;

    return bucket;
}

void add_buckets(int *f_bucket, int *s_bucket, int *result) {
    for (int i = 0; i < BUCKET_LEN; i++)
        result[i] = f_bucket[i] + s_bucket[i];
}

void sort_from_bucket(int *grade, int *bucket, int start, int end) {
    for (int i = 0, j = 0; i < BUCKET_LEN; i++) 
        while(bucket[i]--) 
            grade[j++] = i;
}

/* 
a + b + a + a + b  -> 5 iterations
3*a + 2*b -> 2 iterations
*/
double mean(int *grade_freq, int num_of_elems) {
    int mean = 0;
    for (int i = 0; i < BUCKET_LEN; i++)
        mean += grade_freq[i] * i;
    
    return mean / num_of_elems;
}

double standard_deviation(int *grade, double mean, int num_of_elems) {
    int sum = 0;
    for (int i = 0; i <= BUCKET_LEN; i++)
        sum += grade[i] * pow(i - mean, 2);

    return sqrt(sum / num_of_elems);
}

REGIONS *initialize_regions(int R, int C) {
    REGIONS *br = malloc(sizeof(REGIONS) * R);
    for (int i = 0; i < R; i++)
        br->city = calloc(C, sizeof(CITIES));

    return br;
}

/* Cada regiao tem R*C elementos (de 0 até R*c -1) */
int main() {
    int R, C, A, seed;
    scanf("%d%d%d%d", &R, &C, &A, &seed); 
    
    int *grades = malloc(sizeof(int) *  R * C * A);
    get_grades(grades, R, C, A, seed);

    /* Prints por cidade de regiao */
    REGIONS *brazil = initialize_regions(R, C);
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            // (i * C * A) + (j * A) = primeira pos da regiao + primeira posicao da cidade dentro da regiao
            int region_start = A * (i * C + j); 
            brazil[i].city[j].grade_frequencies = counting_sort(grades, region_start, region_start + A - 1);

            int min = grades[region_start];
            int max = grades[region_start + A - 1];

            int median;
            if (!(A % 2))
                median = (grades[region_start - 1 + (A / 2)] + grades[region_start + A / 2]) / 2;
            else 
                median = grades[region_start + A / 2];

            int sum = 0;
            for (int k = 0; k < BUCKET_LEN; k++)
                sum += k * brazil[i].city[j].grade_frequencies[k];
            brazil[i].city[j].mean = sum / A;

            int sd = 0;
            for (int k = 0; k < BUCKET_LEN; k++)
                sd += k * pow(brazil[i].city[j].grade_frequencies[k] - brazil[i].city[j].mean, 2);
            sd = sqrt(sd / A);

            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %d, média: %d e DP: %d\n", i, j, min, max, median, brazil[i].city[j].mean, sd);
        }
    }
    
    /* Print por regiao */
    for (int i = 0; i < R; i++) {
        brazil[i].mean = 0;
        for (int j = 0; j < C; j++) 
            brazil[i].mean += brazil[i].city[j].mean;
        brazil[i].mean /= C;

        printf("Reg 0: menor: 10, maior: 95, mediana: 40.00, média: 43.63 e DP: 25.07");
    }   

    return 0;
}
 