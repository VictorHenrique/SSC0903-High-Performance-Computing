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

typedef struct _best_city {
    double mean;
    int c_idx, r_idx;
} BEST_CITY;

typedef struct _best_region {
    double mean;
    int r_idx;
} BEST_REGION;

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
    #pragma omp parallel for reduction(+:bucket[:BUCKET_LEN])
    for (int i = start; i <= end; i++) {
        bucket[grade[i]]++;
    }

    int s_pos[BUCKET_LEN] = {start}, e_pos[BUCKET_LEN] = { start + bucket[0] };
    for (int i = 1; i < BUCKET_LEN; i++) {
        e_pos[i] = e_pos[i-1] + bucket[i];
        s_pos[i] = e_pos[i-1];
    }
    
    for (int i = 0; i < BUCKET_LEN; i++) {
        int j;
        #pragma omp simd private(j) 
        for (j = s_pos[i]; j < e_pos[i]; j++)
            grade[j] = i;
    }
}


void sort_from_bucket(int *grade, int *bucket, int *s_pos, int *e_pos) {
    int i;
    #pragma omp parallel for private(i) num_threads(4) 
    for (i = 0; i < BUCKET_LEN; i++) {
        int j;
        #pragma omp simd private(j) 
        for (j = s_pos[i]; j < e_pos[i]; j++)
            grade[j] = i;
    }
}

/* frequencia * nota */
double mean(int *grade_freq, int num_of_elems) {
    double mean = 0;
    int i;

    #pragma omp simd reduction(+: mean)
    for (i = 0; i < BUCKET_LEN; i++)
        mean += grade_freq[i] * i;
    
    return mean / num_of_elems;
}

/* frequencia * (nota - media)^2 */
double standard_deviation(int *bucket, double mean, int num_of_elems) {
    double sum = 0;
    #pragma omp simd reduction(+:sum)
    for (int i = 0; i < BUCKET_LEN; i++){
        double temp = i - mean;
        sum += bucket[i] * (temp * temp);
    }
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

    /* Prints por cidade de regiao */
    REGIONS *brazil = initialize_regions(R, C);

    double par_start = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < R; i++) {

        #pragma omp parallel for 
        for (int j = 0; j < C; j++) {
            // (i * C * A) + (j * A) = primeira pos da regiao + primeira posicao da cidade dentro da regiao
            int region_start = (i * A * C) + (j * A); 
		    counting_sort(grades, brazil[i].city[j].grade_frequencies, region_start, region_start + A - 1);

            brazil[i].city[j].min = grades[region_start];
            brazil[i].city[j].max = grades[region_start + A - 1];
            brazil[i].city[j].median = grades_median(grades, region_start, A);

            brazil[i].city[j].mean = mean(brazil[i].city[j].grade_frequencies, A);
            brazil[i].city[j].sd = standard_deviation(brazil[i].city[j].grade_frequencies, brazil[i].city[j].mean, A);
        }
    }
    
    /* https://stackoverflow.com/questions/66664531/reducing-the-max-value-and-saving-its-index */
    #pragma omp declare reduction(best_city : struct _best_city : omp_out = omp_in.mean > omp_out.mean ? omp_in : omp_out)
    BEST_CITY best_city;

    #pragma omp parallel reduction(best_city: best_city)
    for (int i = 0; i < R; i++) {
        #pragma omp parallel for  
        for (int j = 0; j < C; j++) {
            if (brazil[i].city[j].mean > best_city.mean) {
                best_city.r_idx = i;
                best_city.c_idx = j;
                best_city.mean = brazil[i].city[j].mean;
            }     
        } 
    }  
    
    /* Print por regiao */
    #pragma omp parallel for
    for (int i = 0; i < R; i++) { 

        #pragma omp task
        {
        double mean = 0;
        best_city.mean = 0;    
        #pragma omp parallel for reduction(+: mean) 
        for (int j = 0; j < C; j++) 
            mean += brazil[i].city[j].mean;
        brazil[i].mean = mean / C;
        }    

        #pragma omp task 
        {
        int *g_freq = brazil[i].grade_frequencies;
        #pragma omp parallel for reduction(+: g_freq[:BUCKET_LEN]) 
        for (int j = 0; j < C; j++) {
            for (int k = 0; k < BUCKET_LEN; k++)
                g_freq[k] += brazil[i].city[j].grade_frequencies[k];    
        }
        } 
             
        int region_start = i * C * A, region_end = region_start - 1 + (C * A);

        /* Definindo as posições de inserção ordenada */
        int s_pos[BUCKET_LEN] = {region_start}, e_pos[BUCKET_LEN] = { region_start + brazil[i].grade_frequencies[0] };
        for (int j = 1; j < BUCKET_LEN; j++) {
            e_pos[j] = e_pos[j-1] + brazil[i].grade_frequencies[j];
            s_pos[j] = e_pos[j-1];
        }
        
        sort_from_bucket(grades, brazil[i].grade_frequencies, s_pos, e_pos);

        brazil[i].min = grades[region_start];
        brazil[i].max = grades[region_end];
        brazil[i].median = grades_median(grades, region_start, C*A);            
        brazil[i].sd = standard_deviation(brazil[i].grade_frequencies, brazil[i].mean, C*A);
    }  

    #pragma omp declare reduction(best_region : struct _best_region : omp_out = omp_in.mean > omp_out.mean ? omp_in : omp_out)
    BEST_REGION best_region;

    #pragma omp parallel reduction(best_region: best_region)
    for (int i = 0; i < R; i++) {
        if (brazil[i].mean > best_region.mean) {
            best_region.r_idx = i;
            best_region.mean = brazil[i].mean;
        }   
    }  

    /* Print do Brasil */
    double mean = 0;
    int bucket[BUCKET_LEN] = {0};
    #pragma omp parallel reduction(+: bucket[:BUCKET_LEN])
    for (int i = 0; i < R; i++) {
        
        #pragma omp for
        for (int j = 0; j < BUCKET_LEN; j++)
            bucket[j] += brazil[i].grade_frequencies[j];    
    }   

    #pragma omp parallel for reduction(+: mean)
    for (int i = 0; i < R; i++) 
        mean += brazil[i].mean;    
    mean /= R;

    int s_pos[BUCKET_LEN] = {0}, e_pos[BUCKET_LEN] = {bucket[0]};
    for (int i = 1; i < BUCKET_LEN; i++) {
        e_pos[i] = e_pos[i-1] + bucket[i];
        s_pos[i] = e_pos[i-1];
    }
    sort_from_bucket(grades, bucket, s_pos, e_pos);

    int min = grades[0], max = grades[R*C*A - 1];
    double median = grades_median(grades, 0, R * C * A);
    double sd = standard_deviation(bucket, mean, R*C*A);

    double par_end = omp_get_wtime();

    // print_cities(brazil, R, C);
    // print_regions(brazil, R);

    printf("Brasil: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n\n", min, max, median, mean, sd);

    /* Melhor cidade e melhor regiao */
    printf("Melhor região: Região %d\n", best_region.r_idx);
    printf("Melhor cidade: Região %d, Cidade %d\n", best_city.r_idx, best_city.c_idx);

    printf("Tempo: %fs\n", par_end - par_start);

    free_memory(brazil, R, C, A);
    free(grades);

    return 0;
}
 
