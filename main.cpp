#include <mpi_utils.hpp>
#include <jacobi.hpp>

int main(int argc, char* argv[]) 
{
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

    t = time(NULL);
    do {
        matvec(B, g, xold, x, m, m);
        diff = evalDiff(xold, x, m);
        I ++;
        printf("diff=%lf, eps = %lf\n", diff, eps);
        memcpy(xold, x, m * sizeof(double));
    } while ((diff >= eps) && (I <= MAXITERS));
    t = time(NULL) - t;

    printf("%d iterations consumed %lf seconds\n", I, t);
}