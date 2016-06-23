#include<cstdlib>
#include<fstream>
#include<sys/time.h>
#include<memory.h>

using namespace std;

#define MAX_ITERATION_NUM 1000
#define EPSILON 1e-6

int main(int argc, char* argv[])
{
    int N = atoi(argv[1]);     //demension of matrix

    timeval start, end;
    double elapsedTime;

    double *b = new double[N];
    double **A = new double*[N];
    for(int i = 0; i < N; i++)
        A[i] = new double[N];
    double *p = new double[N];
    double *r = new double[N];
    double *Ap = new double[N];
    double *x = new double[N];
    double r2, beta, alpha, res;


    //read data from file
    ifstream input("input.dat");

    if(!input.is_open())
    {
        cout << "File 'input.dat' can't be opened.\n";
        return 0;
    }

    //read in A
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            input >> A[i][j];
        }
    }

    //read in b
    for(int i = 0; i < N; i++)
    {
        input >> b[i];
    }

    //initialize x_0, r_0, p_0, etc
    memset(x, 0, sizeof(double) * N);
    memcpy(r, b, sizeof(double) * N);
    memcpy(p, b, sizeof(double) * N);
    r2 = 0

    //start of iteration
    gettimeofday(&start, NULL);
    for(int it = 0; it < MAX_ITERATION_NUM; it++)
    {
        //compute Ap
        memset(Ap, 0, sizeof(double) * N);
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                Ap[i] += A[i][j] * p[j];
            }
        }

        //compute alpha
        double pr = 0;
        double pAp = 0;
        for(int i = 0; i < N; i++)
        {
            pr += p[i] * r[i];
            pAp += p[i] * Ap[i];
        }
        alpha = pr / pAp;

        //compute x_k+1, r_k+1
        for(int i = 0; i < N; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        //compute res
        res = 0;
        for(int i = 0; i < N; i++)
        {
            res += abs(r[i]);
        }
        if(res < EPSILON)
            break;

        //compute r
        double _r2 = r2;
        r2 = 0;
        for(int i = 0; i < N; i++)
        {
            r2 += r[i] * r[i];
        }

        //compute beta
        beta = r2 / _r2;

        //compute p_k+1
        for(int i = 0; i < N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
    }
    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000;
    printf("Serial Time: %fms\n", elapsedTime);

    //release memory
    delete[] b;
    for(int i = 0; i < N; i++)
        delete[] A[i];
    delete[] A;
    delete[] p;
    delete[] r;
    delete[] Ap;
    delete[] x;
}
