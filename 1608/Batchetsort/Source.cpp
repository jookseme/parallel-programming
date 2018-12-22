#include <iostream>
#include <time.h>
#include <mpi.h>


using namespace std;

int* CreateMatrix(int n, int m)
{
    int* matrix = new int[n * m];
    for (int i = 0; i < n * m; i++)
        matrix[i] = rand() % 10;
    return matrix;
}


void OutputMatrix(int n, int m, int* matrixToOutput)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << matrixToOutput[i*m + j] << " ";

        cout << endl;
    }
}

int main(int argc, char **argv)
{
    int procNum = 0, procRank = 0;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int n = 5;
    int m = 5;
    int s = 5;

    int* m1 = new int[n*m];
    int* m2 = new int[m*s];
    int* m3 = new int[n*s];
    int* m4 = new int[n*s];

    if (procRank == 0)
    {
        cout << procNum << " " << procRank << endl;
        srand(time(0));

        m1 = CreateMatrix(n, m);
        m2 = CreateMatrix(m, s);
        m3 = CreateMatrix(n, s);
        m4 = CreateMatrix(n, s);

        OutputMatrix(m, n, m1);
        cout << endl;
        OutputMatrix(n, s, m2);
        cout << endl;

        startTime = MPI_Wtime();


        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < s; j++)
            {
                int tmp = 0;
                for (int k = 0; k < m; k++)
                {
                    tmp += m1[i * m + k] * m2[k * s + j];
                }
                m3[i * s + j] = tmp;
            }
        }
        endTime = MPI_Wtime();

        cout << "Linear result matrix: " << endl;
        OutputMatrix(n, s, m3);
        cout << endl << "Linear time:" << endTime - startTime << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    startTime = MPI_Wtime();

    int* send_quantity = new int[procNum]; 
    int* send_displs = new int[procNum]; 

    int* rec_quantity = new int[procNum]; 
    int* rec_displs = new int[procNum]; 
    int rowCount; 
    int* matrix1par, *matrix2par;

    if (procRank == 0)
    {
        int ch = n / procNum;
        int r = n % procNum;
        cout << "ch: " << ch << " r: " << r;
        if (procRank < r) rowCount = ch + 1;
        else rowCount = ch;

        for (int i = 0; i < procNum; i++)
        {
            send_quantity[i] = ch * m;
            send_displs[i] = i * ch * m;
            rec_quantity[i] = ch * s;
            rec_displs[i] = i * ch * s;

        }
        r = n - (procNum - 1) * ch;
        send_quantity[procNum - 1] = r * m;
        rec_quantity[procNum - 1] = r * s;
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(m2, m * s, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(send_quantity, procNum, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(send_displs, procNum, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(rec_quantity, procNum, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(rec_displs, procNum, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);




    matrix1par = new int[send_quantity[procRank]];
    matrix2par = new int[rec_quantity[procRank]];


    MPI_Scatterv(m1, send_quantity, send_displs, MPI_INT, matrix1par, send_quantity[procRank], MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < rowCount; i++)
    {
        for (int k = 0; k < s; k++)
        {
            matrix2par[i * s + k] = 0;
            for (int j = 0; j < m; j++)
            {
                matrix2par[i * s + k] += matrix1par[i * m + j] * m2[j * s + k];
            }
        }
    }


    MPI_Gatherv(matrix2par, rec_quantity[procRank], MPI_INT, m4, rec_quantity, rec_displs, MPI_INT, 0, MPI_COMM_WORLD);
    endTime = MPI_Wtime();

    if (procRank == 0)
    {
        endTime = MPI_Wtime();
        cout << endl << "MPI result martrix:" << endl;
        OutputMatrix(n, s, m4);
        cout << endl << "MPI time: " << endTime - startTime << endl;
    }

    delete[] send_displs;
    delete[] send_quantity;
    delete[] rec_displs;
    delete[] rec_quantity;

    delete[] m1;
    delete[] m2;
    delete[] m3;
    delete[] m4;
    MPI_Finalize();
}