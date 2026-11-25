# HPC_1

`CTRL + SHIFT + B` - build using tasks.json

`mpirun -n <n> ./<program>.0` - запуск программы `program` с `n` процессами

Варианты запуска `mpi_complex_main.cpp`:

`mpirun -n 2 ./mpi_complex_main.o 0.000001 100 5 ./input/matrix.txt ./input/vector.txt`

`mpirun -n 2 ./mpi_complex_main.o 0.000001 100 5`