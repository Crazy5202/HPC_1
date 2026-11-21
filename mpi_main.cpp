#include <mpi_utils.hpp>
#include <jacobi.hpp>

int main(int argc, char* argv[]) 
{
    if (argc < 4) {
        std::cerr << "Ошибка: укажите аргументы <m>, <eps>, <MAXITERS>.\n";
        return 1;
    }

    int m,np,rk,chunk,i,I = 0,MAXITERS,
        *chunks,*disps;
    double *B,*Bloc,*g,*gloc,*x,*xloc,*xold,t,
        diff,eps;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    if(rk == 0) {
        m = atoi(argv[1]);
        eps = atof(argv[2]);
        MAXITERS = atoi(argv[3]);
        chunk = m / np;
        B = (double*)malloc(m*m*sizeof(double));
    }
    
    MPI_Bcast(&m, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
    MPI_Bcast(&eps,1,MPI_DOUBLE,0,
        MPI_COMM_WORLD);
    MPI_Bcast(&chunk,1,MPI_INT,0,
        MPI_COMM_WORLD);
    chunks = (int*) malloc(np * sizeof(int));
    disps = (int*) malloc(np * sizeof(int));
    g = (double*)malloc(m * sizeof(double));
    x = (double*)malloc(m * sizeof(double));
    xold = (double*)malloc(m *
        sizeof(double));
    if(rk == 0) {
        init(m, B, g, xold);
        t = MPI_Wtime();
    }
    for(i = 0; i < np; i++) {
        disps[i] = i * chunk * m;
        if(i == (np - 1))
            chunks[i] = (m - (np - 1) * chunk)*m;
        else
            chunks[i] = chunk * m;
    }
    Bloc = (double*)malloc(chunks[rk] * sizeof(double));
    MPI_Scatterv(B,chunks,disps,
        MPI_DOUBLE,Bloc,chunks[rk],MPI_DOUBLE,
        0, MPI_COMM_WORLD);
    for(i = 0; i < np; i++) {
        disps[i] = i * chunk;
        if(i == (np - 1))
            chunks[i] = m - (np - 1) * chunk;
        else
            chunks[i] = chunk;
    }
    gloc = (double*)malloc(chunks[rk] * sizeof(double));
    MPI_Scatterv(g,chunks,disps,MPI_DOUBLE,
        gloc,chunks[rk],MPI_DOUBLE,0,
        MPI_COMM_WORLD);
    xloc=(double*)malloc(chunks[rk]*sizeof(double));

    do {
        MPI_Bcast(xold,m,MPI_DOUBLE,0,MPI_COMM_WORLD);
        matvec(Bloc,gloc,xold,xloc, chunks[rk],m);
        MPI_Gatherv(xloc,chunks[rk],MPI_DOUBLE,
            x,chunks,disps,MPI_DOUBLE,0,
            MPI_COMM_WORLD);
        if(rk == 0) {
            diff = evalDiff(xold, x, m);
            printf("diff = %lf, eps = %lf\n",diff, eps);
            memcpy(xold, x, m * sizeof(double));
        }
        MPI_Bcast(&diff,1,MPI_DOUBLE,0,
            MPI_COMM_WORLD);
        I ++;
    } while((diff >= eps) && (I <= MAXITERS));

    if(rk == 0) {
        t = MPI_Wtime() - t;
        printf("%d iterations consumed %lf sec\n", I, t);
    }

    MPI_Finalize();
}