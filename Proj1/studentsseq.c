#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc studentsseq.c -o sseq -lm

/* Os dados retornados
para cada nível são: a menor e a maior nota, a mediana, a média aritmética simples, e o
desvio padrão */

void read_entries(int *grades, int R, int C, int A) {
    for (int i = 0, student = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            for (int k = 0; k < A; k++) {
                scanf("%d", &grades[student++]);
            }
        }
    }
}

int max(int a, int b) {
    return a > b ? a : b;
}

int min(int a, int b) {
    return a < b ? a : b;
}

int *find_min_max(int *grade, int start, int end) {
    int min_max[2] = {-1, 101};
    for (int i = start; i <= end; i++)
        min_max[0] = min(min_max[0], grade[i]), min_max[1] = max(min_max[1], grade[i]);
    
    return min_max;
}

double mean(int *grade, int start, int end) {
    int mean = 0;
    for (int i = start; i <= end; i++)
        mean += grade[i];
    
    return mean / (end - start + 1);
}

double standard_deviation(int *grade, int start, int end, double mean) {
    int sum = 0;
    for (int i = start; i <= end; i++)
        sum += pow(grade[i] - mean, 2);

    return sqrt(sum / (end - start + 1));
}

void bucket_sort(int *grade, int start, int end) {
    int bucket[101] = {0};
    for (int i = start; i <= end; i++) 
        bucket[grade[i]]++;
    
    for (int i = 0, j = 0; i < 101; i++) 
        while(bucket[i]--) 
            grade[j++] = i;
}

int main() {
    int R, C, A, seed;
    scanf("%d%d%d%d", &R, &C, &A, &seed); 
    
    int *grades = malloc(sizeof(int) *  R * C * A);
    read_entries(grades, R, C, A);

    free(grades);

    return 0;
}
