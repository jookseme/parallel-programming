// #include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

void MultiplicateMatrixByVector(double **matrix, int matrixSize, double* unitMatrix, double* vector)
{
    for(int i = 0; i < matrixSize; i++)
    {
        vector[i] = 0;
        for (int j = 0; j < matrixSize; j++)
        {
            vector[i] += matrix[i][j] * unitMatrix[j];
        }
    }
}

void PrintMatrix(double **matrix, int matrixSize)
{
    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void PrintVector(double *vector, int vectorSize)
{
    for (int i = 0; i < vectorSize; i++)
    {
        cout << vector[i] << " ";
    }
    cout << endl;
}

void Jacobi(double **matrix, int matrixSize, double* vector, double* unitMatrix)  //TODO: change args order
{
    const double eps = 0.001;
    double norm, start, end;
    double *tempUnitMatrix = new double[matrixSize];

    for (int i = 0; i < matrixSize; i++)
    {
        tempUnitMatrix[i] = unitMatrix[i];
    }
    int cnt = 0;
    
    do
    {
        for (int i = 0; i < matrixSize; i++)
        {
            tempUnitMatrix[i] = vector[i];
            for (int j = 0; j < matrixSize; j++)
            {
                if (i != j)
                {
                    tempUnitMatrix[i] -= matrix[i][j] * unitMatrix[j];
                }
            }
            tempUnitMatrix[i] /= matrix[i][i];
        }

        norm = abs(unitMatrix[0] - tempUnitMatrix[0]);

        for (int k = 0; k < matrixSize; k++)
        {
            if (abs(unitMatrix[k] - tempUnitMatrix[k]) > norm)
            {
                norm = abs(unitMatrix[k] - tempUnitMatrix[k]);
            }
            unitMatrix[k] = tempUnitMatrix[k];
        }
        cnt++;
        cout << "iteration #" << cnt << endl;
    }
    while (norm > eps);

    cout << "Iterations amount: " << cnt << endl;
    delete[] tempUnitMatrix;
}

void GenerateMatrixAndVector(double **&matrix, int &matrixSize, double *&vector)
{
    srand(time(NULL));
    matrix = new double *[matrixSize];
    for (int i = 0; i < matrixSize; i++)
    {
        matrix[i] = new double[matrixSize];
    }
    vector = new double[matrixSize];
    
    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            // matrix[i][j] = rand() % 100 + 50;
            matrix[i][j] = rand() % 10;
        }
        // vector[i] = rand() % 100 + 50;
        vector[i] = rand() % 10;
    }
}

int main(int argc, char *argv[])
{
    double **matrix, *vector, *unitMatrix, *anotherVector;
    int size = 4;

    GenerateMatrixAndVector(matrix, size, vector);

    cout << "Matrix\n";
    PrintMatrix(matrix, size);

    cout << "Vector B\n";
    PrintVector(vector, size);

    anotherVector = new double[size];
    for (int i = 0; i < size; i++)
    {
        anotherVector[i] = 1.0;
    }
    unitMatrix = new double[size];
    cout << endl << vector[1] << endl;

    Jacobi(matrix, size, vector, anotherVector);
    cout << "Result | x\n";
    PrintVector(anotherVector, size);

    cout << "A * x = b | b\n";
    MultiplicateMatrixByVector(matrix, size, anotherVector, unitMatrix);
    PrintVector(unitMatrix, size);

    delete[] matrix;
    delete[] vector;
    delete[] unitMatrix;
    delete[] anotherVector;

    return 0;
}
