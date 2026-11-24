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