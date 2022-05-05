#include "doolittle.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/*
 * Print the matrix
 */
void matrix_print(const char *name, struct matrix m) {
    printf("%s = \n", name);
    for(int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.n; j++) {
            if (j == 0) 
                printf("\t|%7.2lf", m.m[i][j]);
            else
                printf("%8.2lf", m.m[i][j]);
        }
        printf("|\n");
    }
}

/*
 * Print the vector
 */
void vector_print(const char *name, double *x, int n) {
    printf("%s = \n", name);
    for (int i = 0; i < n; i++) {
        printf("\t|%7.2lf|\n", x[i]);
    }
} 

/*
 * 将矩阵A进行LU分解
 */
void decom_LU(struct matrix A, struct matrix *L, struct matrix *U) {
    // 1. 计算 U1i 和 Li1.
    for (int i = 0; i < A.n; i++) {
        U->m[0][i] = A.m[0][i];
        L->m[i][0] = A.m[i][0] / U->m[0][0];
    }

    // 2. 计算U的第r行, L的第r列元素 (r = 1,2,..,n-1)
    for (int r = 1; r < A.n; r++) {
        for (int i = r; i < A.n; i++) {
            double Uri = A.m[r][i];
            double Lir = A.m[i][r];
            for (int k = 0; k <= r -1; k++) {
                Uri -= (L->m[r][k] * U->m[k][i]);
                Lir -= (L->m[i][k] * U->m[k][r]);
            }
            U->m[r][i] = Uri;
            L->m[i][r] = Lir / U->m[r][r];
        }
    }

}

/*
 * 经过 LU分解, 等式 'Ax = b' 转化为了 'LUx = b',
 * 令 'y = Ux', 先通过等式 'Ly = b', 解出 y,
 * 然后通过 'Ux = y' 解出 x. 
 * 返回列向量 x
 */
double* solve_equation(struct matrix L, struct matrix U, double *b) {
    int n = L.n;

    // 1. Ly = b, get y
    double *y = (double *)malloc(n * sizeof(double));
    y[0] = b[0];
    for (int i = 1; i < n; i++) {
        double yi = b[i];
        for (int k = 0; k <= i - 1; k++) {
            yi -= (L.m[i][k] * y[k]);
        }
        y[i] = yi;
    }

    // 2. Ux = y, get x
    double *x = (double *)malloc(n * sizeof(double));
    x[n-1] = y[n-1] / U.m[n-1][n-1];
    
    for (int i = n - 2; i >=0 ; i--) {
        double xi = y[i];
        for (int k = i + 1; k < n; k++) {
            xi -= (U.m[i][k] * x[k]);
        }
        x[i] = xi / U.m[i][i];
    }

    free(y);

    return x;
}

/*
 * Init matrix
 */
bool matrix_init(struct matrix **m, int n) {
    *m = (struct matrix*)malloc(sizeof(struct matrix));
    if (NULL == *m) {
        return false;
    }

   (*m)->n = n;

   (*m)->m = (double **)malloc(n * sizeof(double *));
   if (NULL == (*m)->m) {
       free(*m);
       return false;
   }
   for (int i = 0; i < n; i++) {
       (*m)->m[i] = (double *)malloc(n * sizeof(double));
       if (NULL == (*m)->m[i]) {
           free(*m);
           free((*m)->m);
           return false;
       }
   }

   return true;
}
void matrix_free(struct matrix *m) {
    for (int i = 0; i < m->n; i++) {
        free(m->m[i]);
    }
    free(m->m);
    free(m);
}

/*
 * Main process
 * 使用 doolittle 方法求解非线性齐次方程组
 * Input : 矩阵A, 列向量b
 * Return: 解x 
 */
double *doolittle(struct matrix *A, double *b) {
    struct matrix *L, *U;
    double *x;

    if (false == matrix_init(&L, 3)) {
        printf("Allocate memory error\n");
        exit(-1);
    }
    if (false == matrix_init(&U, 3)) {
        printf("Allocate memory error\n");
        exit(-1);
    }

    system("clear");
    printf("============ Initial Conditions ==================\n");
    matrix_print("A", *A);
    vector_print("b", b, A->n);

    printf("============ LU Decomposition ====================\n");
    decom_LU(*A, L, U); // LU分解

    matrix_print("L", *L);
    matrix_print("U", *U);

    printf("============ Doolittle Process -==================\n");
    x = solve_equation(*L, *U, b);
    
    vector_print("x", x, L->n);

    printf("\n============ Process End =========================\n");
    printf("\n\n\n");
    
    matrix_free(L);
    matrix_free(U);
    matrix_free(A);

    return x;
}

    
int main() {
    struct matrix *A;
    double b[3] = {14, 18, 20};
    double *x = NULL;
    
    if (false == matrix_init(&A, 3)) {
        printf("Allocate memory error\n");
        exit(-1);
    }
    
    A->m[0][0] = 1;
    A->m[0][1] = 2;
    A->m[0][2] = 3;
    A->m[1][0] = 2;
    A->m[1][1] = 5;
    A->m[1][2] = 2;
    A->m[2][0] = 3;
    A->m[2][1] = 1;
    A->m[2][2] = 5;

    x = doolittle(A, b);

    free(x); // User free
    return 0;
}
