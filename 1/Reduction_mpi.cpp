#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include <mpi.h>

#define epsilon 1.e-8
#define MASTER 0        //MASTER进程的rank编号

using namespace std;

int main(int argc, char *argv[]) {

    int M, N;

    string T, P, Db;
    M = atoi(argv[1]);
    N = atoi(argv[2]);	

    double elapsedTime, elapsedTime2;
    timeval start, end, end2;

    if (argc > 3) {
        T = argv[3];
        if (argc > 4) {
            P = argv[4];
            if (argc > 5) {
                Db = argv[5];
            }
        }
    }
    // cout<<T<<P<<endl;

    //double **U_t;
    double *U_t;
    //double alpha, beta, gamma, **Alphas, **Betas, **Gammas;
    double alpha, beta, gamma, *Alphas, *Betas, *Gammas;

    int acum = 0;
    int temp1, temp2;

    /*
     * 这里需要说明:
     * 我们假设这个算法没有编写错误,不修改原算法
     * 所以按照这个算法的隐藏条件,M必须小于等于N才行,否则会引起崩溃
     * */
    if (M > N)
    {
        cout << "M must smaller than N, or the original algorithm will go crash!" << endl;
        return 0;
    }

    //声明为数组对于数据传递更方便
    U_t = new double[N * N];
    Alphas = new double[N * N];
    Betas = new double[N * N];
    Gammas = new double[N * N];
    printf("Arrayes[%d^2] initialized!\n", N);


    MPI_Init(&argc, &argv);
    int rank, rank_size;
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    printf("rank: %d size: %d\n", rank, rank_size);

    //Read from file matrix, if not available, app quit
    //Already transposed
    //Only MASTER rank can read the file and initialize U_t
    if (rank == MASTER)
    {
        cout << "Preparing to read file 'matrix'...\n";
        ifstream matrixfile("matrix");
        if(!(matrixfile.is_open())){
            cout<<"Error: file not found"<<endl;
            return 0;
        }

        for(int i = 0; i < M; i++){
            for(int j =0; j < N; j++){
                //matrixfile >> U_t[i][j];
                matrixfile >> U_t[i * N + j];
            }
        }

        matrixfile.close();
        printf("Read file 'matrix' completed!\n");
    }

    //分发数据
    printf("@@@ %d @@@ Before Broadcast...\n", rank);
    MPI_Bcast(U_t, N * N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    printf("@@@ %d @@@ Broadcasting...\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);    //同步所有进程,确保每个进程都拿到了数据之后在进行下一步
    printf("@@@ %d @@@ Broadcasted!\n", rank);

    //每个rank计算几行,乘回去会大于M
    int rowPerRank = (M + rank_size - 1) / rank_size;
    printf("@@@ %d @@@ rowPerRank = %d\n", rank, rowPerRank);

    MPI_Request *requests_a;
    MPI_Request *requests_b;
    MPI_Request *requests_g;
    requests_a = new MPI_Request[rank_size];
    requests_b = new MPI_Request[rank_size];
    requests_g = new MPI_Request[rank_size];
    //主进程异步收集数据
    if (rank == MASTER)
    {
        printf("@@@ %d @@@ before set Irecv\n", rank);
        int _size = rowPerRank * M;     //每次接收的数据量/每个进程会发送的数据量
        for(int i = 1; i < rank_size; i++)      //自己的不用收集,所以i从1开始
        {
            printf("@@@ %d @@@ setting Irecv for %d\n", rank, i);
            MPI_Irecv(&Alphas[i * _size], _size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &requests_a[i-1]);
            MPI_Irecv(&Betas[i * _size], _size, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &requests_b[i-1]);
            MPI_Irecv(&Gammas[i * _size], _size, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &requests_g[i-1]);
            printf("@@@ %d @@@ Irecv for %d set!\n", rank, i);
        }
        printf("@@@ %d @@@ set Irecv!\n", rank);
    }

    /* Reductions 所有进程都参与计算*/

    printf("@@@ %d @@@ Start computing...\n", rank);
    gettimeofday(&start, NULL);
    /*
     * 把最外层循环拆开(即i循环)
     * 这样可以对应将Alphas/Betas/Gammas按照行分成若干部分,便于最后汇总数据
     * */
    for (int i = rank * rowPerRank; i < (rank + 1) * rowPerRank; i++)
    {
        if (i < M)
        {
            for(int j = 0; j < M; j++)
            {
                alpha = 0.0;
                beta = 0.0;
                gamma = 0.0;
                for (int k = 0; k < N; k++)
                {
                    alpha += U_t[i * N + k] * U_t[i * N + k];
                    beta += U_t[j * N + k] * U_t[j * N + k];
                    gamma += U_t[i * N + k] * U_t[j * N + k];
                }
                //三个结果矩阵视作M*M的矩阵,为了让行与行之间能在内存中连续存储,读取会更高效
                Alphas[i * M + j] = alpha;
                Betas[i * M + j] = beta;
                Gammas[i * M + j] = gamma;
            }
        }
    }
    printf("@@@ %d @@@ End computing!\n", rank);

    //除了主进程之外所有进程都要把数据发送给主进程,使用同步send函数即可;主进程则需要等待所有进程传输完毕
    if (rank != MASTER)
    {
        printf("@@@ %d @@@ before send data\n", rank);
        MPI_Request a, b, c;
        int _size = rowPerRank * M;
        printf("@@@ %d @@@ _size = %d\n", rank, _size);
        if (rank != rank_size - 1)      //如果不是最后一个进程,最后一个进程因为可能会有一部分是大于M的行
        {
            MPI_Isend(&Alphas[rank * _size], _size, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &a);
            MPI_Isend(&Betas[rank * _size], _size, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD, &b);
            MPI_Isend(&Gammas[rank * _size], _size, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD, &c);
        }
        else if (rank == rank_size - 1)     //最后一个进程
        {
            int __size = (M - rank * rowPerRank) * M;
            MPI_Isend(&Alphas[rank * _size], __size, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &a);
            MPI_Isend(&Betas[rank * _size], __size, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD, &b);
            MPI_Isend(&Gammas[rank * _size], __size, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD, &c);
        }
        printf("@@@ %d @@@ sending data\n", rank);

        //等待发送结束就可以退出了
        MPI_Wait(&a, MPI_STATUS_IGNORE);
        MPI_Wait(&b, MPI_STATUS_IGNORE);
        MPI_Wait(&c, MPI_STATUS_IGNORE);

        //从属进程退出
        printf("@@@ %d @@@ exit!\n", rank);
        MPI_Finalize();
        return 0;
    }

    //等待主进程收集完所有数据之后继续
    if (rank_size > 1)
    {
        MPI_Waitall(rank_size - 1, requests_a, MPI_STATUS_IGNORE);
        MPI_Waitall(rank_size - 1, requests_b, MPI_STATUS_IGNORE);
        MPI_Waitall(rank_size - 1, requests_g, MPI_STATUS_IGNORE);
    }

    gettimeofday(&end, NULL);


    //Output time and iterations
    if (T == "-t" || P == "-t") {
        elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
        elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
        cout << "Time: " << elapsedTime << " ms." << endl << endl;
    }

    // Output the matrixes for debug
    if (T == "-p" || P == "-p") {
        cout << "Alphas" << endl << endl;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << Alphas[i * M + j] << "  ";
            }
            cout << endl;
        }

        cout << endl << "Betas" << endl << endl;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << Betas[i * M + j] << "  ";
            }
            cout << endl;
        }

        cout << endl << "Gammas" << endl << endl;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << Gammas[i * M + j] << "  ";
            }
            cout << endl;
        }
    }

    //Generate files for debug purpouse
    if (Db == "-d" || T == "-d" || P == "-d")
    {
        ofstream Af;
        //file for Matrix A
        Af.open("AlphasMPI.mat");

        //这一行要注释掉才能通过验证，因为验证程序不会读取MPI输出的文件的这两个参数
        //Af << M << "  " << N;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Af << " " << Alphas[i * M + j];
            }
            Af << "\n";
        }

        Af.close();

        ofstream Uf;

        //File for Matrix U
        Uf.open("BetasMPI.mat");

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Uf << " " << Betas[i * M + j];
            }
            Uf << "\n";
        }
        Uf.close();

        ofstream Vf;
        //File for Matrix V
        Vf.open("GammasMPI.mat");

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vf << " " << Gammas[i * M + j];
            }
            Vf << "\n";
        }


        Vf.close();

        ofstream Sf;
    }

    delete[] Alphas;
    delete[] Betas;
    delete[] Gammas;
    delete[] U_t;

    MPI_Finalize();
    return 0;
}
