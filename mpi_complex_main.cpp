#include <mpi_utils.hpp>
#include <jacobi.hpp>

#include <mpi_utils.hpp>
#include <general_utils.hpp>
#include <jacobi.hpp>

int main(int argc, char* argv[]) 
{
    if (argc != 4 and argc != 6) {
        std::cerr << "Ошибка: укажите аргументы <eps>, <MAXITERS>, <m>. Опционально добавляются: <filepath_mat> + <filepath_vec>.\n";
        return 1;
    }

    int np, rk; // Число процессов и ранг текущего процесса

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np); // узнаём число и ранг
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    int m = 0, MAXITERS = 0, chunk = 0; // Длина x, максимальное число итераций, размер блока для вычислений на отдельном процессе
    double eps = 0.0, t = 0.0, // Погрешность, время вычисления
        *B = NULL,*g = NULL,*x = NULL; // Для Bx=g - матрица B, вектор g, векторы x и xold текущего и прошлого приближений

    if(rk == 0) { // на root
        t = MPI_Wtime(); // запускаем таймер на root

        eps = atof(argv[1]);
        MAXITERS = atoi(argv[2]);

        m = atoi(argv[3]);

        x = (double*)malloc(m * sizeof(double)); // выделяем память
        g = (double*)malloc(m * sizeof(double));
        B = (double*)malloc(m*m*sizeof(double));

        if (argc == 6) { // Получаем m из файлов, если они указаны
            load_matrix(argv[4], B, m);
            load_vector(argv[5], g, m);
        }

        // for(int i = 0; i < m; i ++) {
        //     for(int j = 0; j < m; j ++) {
        //         printf("%lf ", B[i*m+j]);
        //     }
        //     printf("\n");
        // }

        // for(int i = 0; i < m; i ++) {
        //     printf("%lf ", g[i]);
        // }

        // printf("\n");

        chunk = m / np; // вычисляем размер блока для отдельного процесса
    }

    MPI_Bcast(&m, 1, MPI_INT, 0, // передаём гиперпараметры и chunk с root на другие процессы
        MPI_COMM_WORLD);
    MPI_Bcast(&MAXITERS, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
    MPI_Bcast(&eps, 1, MPI_DOUBLE, 0,
        MPI_COMM_WORLD);
    MPI_Bcast(&chunk, 1, MPI_INT, 0,
        MPI_COMM_WORLD);

    double *xold = (double*)malloc(m * sizeof(double)); // xold выделяется на каждом процессе, т.к. нужен при умножении

    if (rk==0) { // только на root
        if (argc == 6) { // Иницилизируем реальные данные либо тестовые
            init(m, B, g, xold);
        } else { // Иницилизируем реальные тестовые данные
            init_test(m, B, g, xold);
        }
    }
    
    double *Bloc = NULL, *gloc = NULL; // Локальные части B и g для использования всеми процессами в вычислениях
    int *chunks = NULL, *disps = NULL; // Массивы размеров блоков и их смещений
    
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
        MPI_DOUBLE,Bloc,chunks[rk],MPI_DOUBLE,0,
        MPI_COMM_WORLD);

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
    double diff = 0.0, *xloc = (double*)malloc(chunks[rk] * sizeof(double));

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
        t = MPI_Wtime() - t; // вычисляем затраченное время
        printf("%d iterations consumed %lf sec\n", i, t);

        save_vector("./result.txt", x, m); // Сохраняем результат
    }

    MPI_Finalize(); // завершаем работу MPI
}