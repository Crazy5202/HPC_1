#pragma once

#include <iostream>

#include <mpi.h>
#include <unistd.h>

/**
* Захватывает процесс в ловушку для дебага. 
Нужно остановить исполнение и сменить i на другое значение для продолжения.
*/
void waitForAttach() 
{
    volatile int i = 0;
    printf("PID %d ready for attach\n", getpid());
    fflush(stdout);
    while (i == 0) {
        sleep(5);
    }
}

/**
* @brief Считывает вектор из файла.
* @param file Файл с вектором.
* @param g Неинициализированная переменная для сохранения вектора.
* @return Число элементов m
*/
int load_vector(char* file, double* g) 
{
    
}

/**
* @brief Считывает матрицу из файла.
* @param file Файл с матрицей.
* @param m Число строк и столбцов матрицы B / длина вектора x.
* @param B Неинициализированная переменная для сохранения матрицы.
*/
void load_matrix(char* file, int m, double* B) 
{
    
}