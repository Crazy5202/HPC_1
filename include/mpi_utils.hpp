#include <iostream>

#include <mpi.h>
#include <unistd.h>

void waitForAttach() 
{
    volatile int i = 0;
    printf("PID %d ready for attach\n", getpid());
    fflush(stdout);
    while (i == 0) {
        sleep(5);
    }
}