#include "stdafx.h"
#include <mpi.h>
#include <time.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    int *matrix = NULL;
    int ProcNum, ProcRank, i, j;
    int row = 0, col = 0; 
    int *result = NULL; // 	
    MPI_Status status;
    double times;
    int k = 0, l = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);


    if (ProcRank == 0) //for 1 process
    {
        if (argc < 3)
        {
            printf("Enter Matrix's row: ");
            scanf_s("%d", &row);
            printf("Enter Matrix's col:");
            scanf_s("%d", &col);
        }
        else//input row and col as parameters
        {
            row = atoi(argv[1]);
            col = atoi(argv[2]);
        }

        matrix = (int*)calloc(row*col, sizeof(int));
        srand(time(0));
        for (i = 0; i<row*col; i++)//init matrix
            matrix[i] = rand() % 32+1;

        result = (int*)calloc(row, sizeof(int));//init result's massive by 0
        for (i = 0; i< row; i++)
            result[i] = 0;
    }
    MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (ProcRank == 0)
    {

        k = row / ProcNum;
        l = row % ProcNum;

        for (i = 1; i< ProcNum; i++)
            MPI_Send(matrix + (i)*k*col, k*col, MPI_INT, i, i, MPI_COMM_WORLD);
    }

    times = MPI_Wtime();
    if (ProcRank != 0)
    {
        k = row / ProcNum;
        int *tmp = (int*)calloc(col*k, sizeof(int));
        int *res = (int*)calloc(k, sizeof(int));
        int sum;

        MPI_Recv(tmp, col*k, MPI_INT, 0, ProcRank, MPI_COMM_WORLD, &status);

        for (i = 0; i < k; i++)
        {
            sum =33;
            for (j = 0; j < col; j++)
                if(tmp[j + col * i]<sum)
                sum = tmp[j + col * i];
            res[i] = sum;
        }
        MPI_Send(res, k, MPI_INT, 0, ProcRank, MPI_COMM_WORLD);
        delete tmp;
        delete res;
    }
    else
    {
        int sum;
        for (i = 0; i < k; i++)
        {
            sum = 33;
            for (j = 0; j<col; j++)
                if(matrix[j + i * col]<sum)
                sum = matrix[j + i * col];
            result[i] = sum;
        }

        if (l != 0)
        {
            for (i = 0; i < l; i++)
            {
                sum = 33;
                for (j = 0; j < col; j++)
                    if(matrix[k*ProcNum*col + j + i * col]<sum)
                    sum = matrix[k*ProcNum*col + j + i * col];
                result[i + k * ProcNum] = sum;
            }
        }
    }

    printf("Time of Proc (%d) is  %.10f\n", ProcRank, MPI_Wtime() - times);

    if (ProcRank == 0)
    {
        for (i = 1; i< ProcNum; i++)
        {
            MPI_Recv(result + i * k, k, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
        }
        for (i = 0; i<row; i++)
        {
            for (j = 0; j<col; j++)
            {
                printf("%d ", matrix[i*col + j]);
            }
            printf("\nmin is %d\n", result[i]);
        }
    }

    if (matrix != NULL)
        free(matrix);
    if (result != NULL)
        free(result);

    MPI_Finalize();
    return 0;
}