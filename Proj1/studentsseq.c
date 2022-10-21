#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc studentsseq.c -o sseq -lm

#define BUCKET_LEN 101
#define GRADE_MAX 101

typedef struct _cities {
    int *grade_frequencies; // per city
    float mean;
} CITIES;

typedef struct _region {
    CITIES *city;
    int *grade_frequencies; // per region
    float mean;
} REGIONS;


/* Geração pseudo-aleatória das notas de cada aluno a partir de uma seed */
void get_grades(int *grades, int R, int C, int A, int seed) {
    srand(seed);
    for (int i = 0, student = 0; i < R; i++) 
        for (int j = 0; j < C; j++) 
            for (int k = 0; k < A; k++) 
                grades[student++] = rand() % GRADE_MAX;    
}

/* Counting Sort */
void counting_sort(int *grade, int *bucket, int start, int end) {
    for (int i = start; i <= end; i++)
        bucket[grade[i]]++;

    for (int i = 0, j = start; i < BUCKET_LEN; i++) 
        for (int k = bucket[i]; k; k--)  
            grade[j++] = i; 
        // while(bucket[i]--) 
        //     grade[j++] = i;
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

double median(int *grades, int first_element, int n_elements){
    int median;
    if (!(n_elements % 2))
        median = (grades[first_element - 1 + (n_elements / 2)] + grades[first_element + n_elements / 2]) / 2;
    else 
        median = grades[first_element + n_elements / 2];

    return median;
}

double standard_deviation(int *grade, double mean, int num_of_elems) {
    int sum = 0;
    for (int i = 0; i <= BUCKET_LEN; i++)
        sum += grade[i] * pow(i - mean, 2);

    return sqrt(sum / num_of_elems);
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

    /* Prints por cidade de regiao */
    REGIONS *brazil = initialize_regions(R, C);
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            // (i * C * A) + (j * A) = primeira pos da regiao + primeira posicao da cidade dentro da regiao
            int region_start = A * (i * C + j); 
            counting_sort(grades, brazil[i].city[j].grade_frequencies, region_start, region_start + A - 1);

            int min = grades[region_start];
            int max = grades[region_start + A - 1];

            float median;
            if (!(A % 2))
                median = (float)(grades[region_start - 1 + (A / 2)] + grades[region_start + A / 2]) / 2;
            else 
                median = grades[region_start + A / 2];

            float sum = 0;
            for (int k = 0; k < BUCKET_LEN; k++)
                sum += k * brazil[i].city[j].grade_frequencies[k];
            brazil[i].city[j].mean = sum / A;

            float sd = 0;
            for (int k = 0; k < BUCKET_LEN; k++)
                sd += brazil[i].city[j].grade_frequencies[k] * pow(k - brazil[i].city[j].mean, 2);
            sd = sqrt(sd / A);

            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, j, min, max, median, brazil[i].city[j].mean, sd);
        }
        printf("\n");
    }
    
    /* Print por regiao */
    printf("\n\n");
    for (int i = 0; i < R; i++) {
        brazil[i].mean = 0;
        for (int j = 0; j < C; j++) {
            add_buckets(brazil[i].grade_frequencies, brazil[i].city[j].grade_frequencies);
            brazil[i].mean += brazil[i].city[j].mean;
        } 
        brazil[i].mean /= C;

        int region_start = i * C * A, region_end = region_start - 1 + (C * A);
        sort_from_bucket(grades, brazil[i].grade_frequencies, region_start, region_end);

        int min = grades[region_start], max = grades[region_end];
        float median;
        if (!((C * A)% 2))
            median = (grades[region_start - 1 + (C * A / 2)] + grades[region_start + (C * A / 2)]) / 2;
        else 
            median = grades[region_start + (C * A / 2)];
            
        float sd = 0;
        for (int k = 0; k < BUCKET_LEN; k++)
            sd += brazil[i].grade_frequencies[k] * pow(k - brazil[i].mean, 2);
        sd = sqrt(sd / A);

        printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, min, max, median, brazil[i].mean, sd);
    }   

    /* Print do Brasil */
    printf("\n\n");
    float mean = 0;
    int bucket[BUCKET_LEN] = {0};
    for (int i = 0; i < R; i++) 
        add_buckets(bucket, brazil[i].grade_frequencies), mean += brazil[i].mean;    
    mean /= R;

    sort_from_bucket(grades, bucket, 0, R*C*A - 1);
    int min = grades[0], max = grades[R*C*A - 1];
    
    float median;
    if (!((R * C * A)% 2))
            median = (grades[(R * C * A) / 2 - 1] + grades[(R * C * A / 2)]) / 2;
        else 
            median = grades[(R * C * A) / 2];

    float sd = 0;
    for (int k = 0; k < BUCKET_LEN; k++)
            sd += bucket[k] * pow(k - mean, 2);
        sd = sqrt(sd / R);

    printf("Brasil: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f", min, max, median, mean, sd);

    /* Melhor cidade e melhor regiao */
    printf("\n\n");
    float best_city[3] = {0}, best_region[2] = {0};
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            if (brazil[i].city[j].mean > best_city[2]) {
                best_city[0] = i;
                best_city[1] = j;
                best_city[2] = brazil[i].city[j].mean;
            }
        }

        if (brazil[i].mean > best_region[1]) {
            best_region[0] = i;
            best_region[1] = brazil[i].mean;
        }
    }
    printf("Melhor região: Região %.0f\nMelhor cidade: Região %.0f, Cidade %.0f\n", best_region[0], best_city[0], best_city[1]);


    free_memory(brazil, R, C, A);
    free(grades);

    return 0;
}
 