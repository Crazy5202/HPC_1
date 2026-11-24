#pragma once

#include <fstream>
#include <iomanip>

/**
 * @brief Считывает матрицу из файла.
 * @param file Файл с матрицей.
 * @param B Инициализированная переменная для хранения вектора.
 * @param m Число строк и столбцов матрицы B.
 */
void load_matrix(char* file, double* B, int m)
{
    std::ifstream f(file);

    if (!f)
        throw std::runtime_error("Matrix file not found or could not be opened");

    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            f >> B[i * m + j];

    f.close();
}

/**
 * @brief Считывает вектор из файла.
 * @param file Файл с вектором.
 * @param g Инициализированная переменная для хранения вектора.
 * @param m Длина вектора g.
 */
void load_vector(char* file, double* g, int m)
{
    std::ifstream f(file);

    if (!f)
        throw std::runtime_error("Vector file not found or could not be opened");


    for (int i = 0; i < m; i++)
        f >> g[i];

    f.close();
}

/**
 * @brief Записывает вектор в файл.
 * @param file Файл для записи вектора.
 * @param g Переменная c вектором.
 * @param m Длина вектора g.
 */
void save_vector(std::string file, double* g, int m)
{
    std::ofstream f(file);
    f << std::setprecision(15);
    f << m << std::endl;

    for (int i = 0; i < m; i++)
        f << g[i] << (i < m - 1 ? ' ' : '\n');

    f.close();
}