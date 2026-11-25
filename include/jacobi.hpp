#pragma once

#include <cmath>

/**
* @brief Подготавливает данные для Ax=f.
* @param m Число строк и столбцов матрицы B / длина вектора x.
* @param B Матрица уравнения.
* @param g Вектор правой части.
* @param x Вектор для начального приближения решения.
*/
void init(int m, double* B, double* g, double* x) 
{
    for(int i = 0; i < m; i ++) {
        g[i] /= B[i*m+i];
        x[i] = 1.;
        for(int j = 0; j < m; j ++) {
            if (i==j) continue;
            B[i * m + j] /= B[i*m+i];
        }
        B[i*m+i] = 0.0;
    }
}

/**
* @brief Использует базовые данные для Ax=f (решение - вектор 1, 2, ... m).
* @param m Длина вектора x.
* @param B Матрица A с элементами, поделёнными на элементы с
главной диагонали.
* @param g Вектор f с элементами, поделёнными на элементы с
главной диагонали.
* @param x Вектор x.
*/
void init_test(int m, double* B, double* g, double* x) 
{
    for(int i = 0; i < m; i ++) {
        g[i] = ((double)m+1)/4.+(1.-1./(double)(2*m))*(i+1);
        x[i] = 1.;
        for(int j = 0; j < m; j ++) {
            B[i * m + j] = 1./(double)(2 * m);
        }
        B[i * m + i] = 0.0;
    }
}

/**
* @brief Рассчитывает евклидово расстояние между векторами.
* @param u,v вектор.
* @param m длина векторов u,v.
* @return Значение расстояния
*/
double evalDiff(double* u, double* v, int m) 
{
    int i;
    double a = 0.0;
    for(i = 0; i < m; i ++) {
        double b;
        b = v[i] - u[i];
        a += b * b;
    }
    return sqrt(a);
}

/**
* @brief Умножает матрицу на вектор.
* @param B Матрица.
* @param g Вектор.
* @param xold Прошлое приближение вектора x.
* @param x Вектор x.
* @param k Число строк матрицы B, длина вектора f (g). Размер чанка при распараллеливании.
* @param m Число столбцов матрицы B, длина вектора x.
*/
void matvec(double* B, double* g, double* xold, 
    double* x, int k, int m) 
{
    int i;
    for(i = 0; i < k; i ++) {
        double a = 0.;
        double* row = B + i * m;
        for(int j = 0; j < m; j ++)
            a += row[j] * xold[j];
        x[i] = -a + g[i];
    }
}