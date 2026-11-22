#include <iostream>

#include <mpi.h>
#include <unistd.h>

void waitForAttach() {
    volatile int i = 0;
    printf("PID %d ready for attach\n", getpid());
    fflush(stdout);
    while (0 == i)
        sleep(5);
}

int main(int argc, char** argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank of the current process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    waitForAttach();

    if (rank == 0) // Process 0 (Sender)
    {
        // Sending an array of integers with a variable length
        int data[] = {10, 20, 30, 40, 50}; // Example data
        int count = 5; // Number of elements
        MPI_Send(data, count, MPI_INT, 1, 0, MPI_COMM_WORLD);
        std::cout << "Process 0 sent " << count << " numbers to Process 1." << std::endl;
    }
    else if (rank == 1) // Process 1 (Receiver)
    {
        MPI_Status status;

        // Probe to check for an incoming message
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

        // Determine the number of integers in the message
        int count;
        MPI_Get_count(&status, MPI_INT, &count);

        // Receive the message
        int data[count]; // Dynamically allocate array
        MPI_Recv(data, count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Print the received data
        std::cout << "Process 1 received " << count << " numbers from Process 0: ";
        for (int i = 0; i < count; i++)
        {
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}