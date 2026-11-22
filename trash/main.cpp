#include <mpi_utils.hpp>
#include <jacobi.hpp>

#include <chrono>

int main(int argc, char* argv[]) 
{
    auto start = std::chrono::steady_clock::now();

    std::ios::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    if (argc < 4) {
        std::cerr << "Ошибка: укажите аргументы <m>, <eps>, <MAXITERS>.\n";
        return 1;
    }

    int m, I = 0, MAXITERS;
    double *B, *g, *x, *xold, 
        eps, diff, t;

    m = atoi(argv[1]);
    eps = atof(argv[2]);
    MAXITERS = atoi(argv[3]);

    B = (double*)malloc(m * m * sizeof(double));
    g = (double*)malloc(m * sizeof(double));
    xold = (double*)malloc(m * sizeof(double));
    x = (double*)malloc(m * sizeof(double));

    init(m, B, g, xold);

    do {
        matvec(B, g, xold, x, m, m);
        diff = evalDiff(xold, x, m);
        I ++;
        printf("diff=%lf, eps = %lf\n", diff, eps);
        memcpy(xold, x, m * sizeof(double));
    } while ((diff >= eps) && (I <= MAXITERS));

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    printf("%d iterations consumed %lf seconds\n", I, (duration.count()/1e6));
}