#include <mpi_utils.hpp>
#include <jacobi.hpp>

int main(int argc, char* argv[]) 
{
    if (argc < 4) {
        std::cerr << "Ошибка: укажите аргументы <m>, <eps>, <MAXITERS>.\n";
        return 1;
    }

    int np, rk; // Число процессов и ранг текущего процесса

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np); // узнаём число и ранг
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    int m,MAXITERS,chunk; // Длина x, максимальное число итераций, размер блока для вычислений на отдельном процессе
    double eps, t; // Погрешность, время вычисления

    if(rk == 0) { // на root
        t = MPI_Wtime(); // запускаем таймер на root

        m = atoi(argv[1]); // считываем гиперпараметры для вычислений из консоли
        eps = atof(argv[2]);
        MAXITERS = atoi(argv[3]);

        chunk = m / np; // вычисляем размер блока для вычислений на отдельном процессе
    }
    
    MPI_Bcast(&m, 1, MPI_INT, 0, // передаём гиперпараметры для вычислений с root на другие процессы
        MPI_COMM_WORLD);
    MPI_Bcast(&eps,1,MPI_DOUBLE,0,
        MPI_COMM_WORLD);
    MPI_Bcast(&chunk,1,MPI_INT,0,
        MPI_COMM_WORLD);
    
    double *B,*g,*x, // Для Bx=g - матрица B, вектор g, векторы x и xold текущего и прошлого приближений
        *xold = (double*)malloc(m * sizeof(double)); // xold в отличие от x приходится выделять на каждом процессе

    if(rk == 0) { // только для root
        B = (double*)malloc(m*m*sizeof(double)); // выделяем память
        g = (double*)malloc(m * sizeof(double));
        x = (double*)malloc(m * sizeof(double));

        init_test(m, B, g, xold); // заполняем B, g и начальное приближение xold
    }

    double *Bloc, *gloc; // Локальные части B и g для использования всеми процессами в вычислениях
    int *chunks, *disps; // Массивы размеров блоков и их смещений
    
    chunks = (int*) malloc(np * sizeof(int));
    disps = (int*) malloc(np * sizeof(int));

    for(int i = 0; i < np; i++) { // настройка раздачи матрицы B
        disps[i] = i * chunk * m; // смещение i-го блока матрицы: chunk строк по m элементов
        if(i == (np - 1)) // последний процесс
            chunks[i] = (m - (np - 1) * chunk)*m; // остаток строк матрицы
        else
            chunks[i] = chunk * m; // chunk строк по m элементов
    }

    Bloc = (double*)malloc(chunks[rk] * sizeof(double)); // выделяем память

    MPI_Scatterv(B,chunks,disps, // раздаём матрицу B из root в процессы
        MPI_DOUBLE,Bloc,chunks[rk],MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    for(int i = 0; i < np; i++) { // настройка раздачи вектора g
        disps[i] = i * chunk; // смещение i-го блока вектора: chunk элементов
        if(i == (np - 1)) // последний процесс
            chunks[i] = m - (np - 1) * chunk; // остаток элементов вектора
        else
            chunks[i] = chunk; // chunk элементов
    }

    gloc = (double*)malloc(chunks[rk] * sizeof(double)); // выделяем память

    MPI_Scatterv(g,chunks,disps,MPI_DOUBLE, // раздаём вектор g из root в процессы
        gloc,chunks[rk],MPI_DOUBLE,0,
        MPI_COMM_WORLD);

    int i = 0;
    double diff, *xloc = (double*)malloc(chunks[rk] * sizeof(double));

    do {
        MPI_Bcast(xold, m, MPI_DOUBLE,0,MPI_COMM_WORLD); // root рассылает всем предыдущее приближение xold

        matvec(Bloc,gloc,xold,xloc,chunks[rk], m); // каждый процесс считает свою часть xloc = -Bloc*xold + gloc

        MPI_Gatherv(xloc,chunks[rk],MPI_DOUBLE, // собираем локальные части xloc на root в полный вектор x
            x,chunks,disps,MPI_DOUBLE,0,
            MPI_COMM_WORLD);

        if(rk == 0) { // только root
            diff = evalDiff(xold, x, m); // считаем евклидову норму разности xold и x
            printf("diff = %lf, eps = %lf\n", diff, eps);
            memcpy(xold, x, m * sizeof(double)); // копируем новое приближение x в xold для следующей итерации
        }

        MPI_Bcast(&diff,1,MPI_DOUBLE,0, // root рассылает значение diff всем процессам
            MPI_COMM_WORLD);

        ++i; // увеличиваем счётчик выполненных итераций
        
    } while((diff >= eps) && (i < MAXITERS)); // продолжаем, пока не достигнута точность и не превышен лимит итераций

    if(rk == 0) { // только на root
        t = MPI_Wtime() - t; // вычисляем затраченное время, печатаем
        printf("%d iterations consumed %lf sec\n", i, t);
        for (int i = 0; i<5; ++i) printf("%lf ", x[i]);
        fflush(stdout);
    }

    MPI_Finalize(); // завершаем работу MPI
}