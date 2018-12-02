#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

using namespace std;

void ParallelJacobiMethod(double *coefficientMatrix, double *resultingValuesPerProc, double *freeMembersColumn, int dimension,
                          double EPS, int MAX_ITERATIONS, int rowsPerProc, int procRank);
void LinearJacobiMethod(double *coefficientMatrix, double *resultingValuesPerProc, double *freeMembersColumn, int dimension,
                        double EPS, int MAX_ITERATIONS);

void GenerateMatrix(double *&matrix, int dimension);
void GenerateVector(double *&vector, int dimension);

void PrintLinearEquationSystem(double *coefficientMatrix, double *freeMembersColumn, int dimension);
void PrintSolutions(double *solutions, int dimension);

// Нахождение текущей точности вычислений
double FindMatrixNorm(double *previousResultingValues, double *currentResultingValues, int dimension);
// Проверка эквивалентности полученных результатов
bool AreSolutionsEqual(double *linearSolutions, double *parallelSolutions, int dimension, double EPS);

int main(int argc, char *argv[])
{
    const double EPS = 0.001;       // Точность вычислений
    const int MAX_ITERATIONS = 100; // Предельное число шагов в методе Якоби
    int dimension;                  // Разменость матрицы
    int rowsPerProc;                // Число строк матрицы, отсылаемых каждому процессу
    int procNumber, procRank;
    double startTime = 0.0, endTime = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNumber);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    if (procRank == 0)
    {
        dimension = (argc == 2) ? atoi(argv[1]) : 1000;
        rowsPerProc = dimension / procNumber;
    }

    MPI_Bcast(&dimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowsPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double *coefficientMatrix = new double[dimension * dimension];
    double *freeMembersColumn = new double[dimension];

    if (procRank == 0)
    {
        srand(time(NULL));
        GenerateMatrix(coefficientMatrix, dimension);
        GenerateVector(freeMembersColumn, dimension);

        if (dimension < 10)
        {
            PrintLinearEquationSystem(coefficientMatrix, freeMembersColumn, dimension);
        }
    }

    double *elementsPerProc = new double[rowsPerProc * dimension];
    double *freeMembersPerProc = new double[rowsPerProc];

    startTime = MPI_Wtime();

    // Рассылка частей матрицы коэффициентов и свободных членов каждому процессу
    MPI_Scatter(coefficientMatrix, rowsPerProc * dimension, MPI_DOUBLE, elementsPerProc, rowsPerProc * dimension,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(freeMembersColumn, rowsPerProc, MPI_DOUBLE, freeMembersPerProc, rowsPerProc, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    double *resultingValuesPerProc = new double[dimension];
    // Решение системы линейных уравнений с помощью паралелльного метода Якоби
    ParallelJacobiMethod(elementsPerProc, resultingValuesPerProc, freeMembersPerProc, dimension, EPS, MAX_ITERATIONS,
                         rowsPerProc, procRank);

    double *resultingValues = new double[dimension];
    // Сбор результатов
    MPI_Gather(resultingValuesPerProc, rowsPerProc, MPI_DOUBLE, resultingValues, rowsPerProc, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    endTime = MPI_Wtime();

    delete[] resultingValuesPerProc;
    delete[] elementsPerProc;
    delete[] freeMembersPerProc;

    if (procRank == 0)
    {
        cout << "Parallel Jacobi method: " << setprecision(3) << endTime - startTime << " sec" << endl;

        if (dimension < 10)
        {
            cout << "Parallel solutions: ";
            PrintSolutions(resultingValues, dimension);
        }

        double *linearResultingValues = new double[dimension];

        startTime = MPI_Wtime();

        LinearJacobiMethod(coefficientMatrix, linearResultingValues, freeMembersColumn, dimension, EPS, MAX_ITERATIONS);

        endTime = MPI_Wtime();

        delete[] coefficientMatrix;
        delete[] freeMembersColumn;

        cout << "Linear Jacobi method: " << setprecision(3) << endTime - startTime << " sec" << endl;

        if (dimension < 10)
        {
            cout << "Linear solutions:   ";
            PrintSolutions(linearResultingValues, dimension);
        }

        // Проверка эквивалентности линейного и параллельного решений
        if (AreSolutionsEqual(linearResultingValues, resultingValues, dimension, EPS))
        {
            cout << "Solutions are equal" << endl;
        }
        else
        {
            cout << "Solutions are not equal" << endl;
        }

        delete[] linearResultingValues;
    }

    delete[] resultingValues;

    MPI_Finalize();
    return 0;
}

void ParallelJacobiMethod(double *coefficientMatrix, double *resultingValuesPerProc, double *freeMembersColumn,
                          int dimension, double EPS, int MAX_ITERATIONS, int rowsPerProc, int procRank)
{
    double *previousResultingValues = new double[dimension];
    double *currentResultingValues = new double[dimension];

    // На первом шаге метода Якоби за текущее решение берется столбец свободных членов
    MPI_Allgather(freeMembersColumn, rowsPerProc, MPI_DOUBLE, currentResultingValues, rowsPerProc,
                  MPI_DOUBLE, MPI_COMM_WORLD);

    int iterationNumber = 0;
    int startIndex = 0;
    do
    {
        iterationNumber++;

        // Своппинг currentResultingValues и previousResultingValues
        double *temp;
        temp = previousResultingValues;
        previousResultingValues = currentResultingValues;
        currentResultingValues = temp;

        for (int i = 0; i < rowsPerProc; i++)
        {
            startIndex = i + procRank * rowsPerProc;
            resultingValuesPerProc[i] = freeMembersColumn[i];
            for (int j = 0; j < startIndex; j++)
            {
                resultingValuesPerProc[i] -= coefficientMatrix[i * dimension + j] * previousResultingValues[j];
            }
            for (int j = startIndex + 1; j < dimension; j++)
            {
                resultingValuesPerProc[i] -= coefficientMatrix[i * dimension + j] * previousResultingValues[j];
            }
            // j Не может быть равным startIndex, так как на этом шаге текущий элемент coefficientMatrix
            // принадлежит главной диагонали, а значит по условиям метода на него нужно делить
            resultingValuesPerProc[i] /= coefficientMatrix[i * dimension + startIndex];
        }

        // Сбор полученных на текущем шаге решений в единый массив currentResultingValues
        MPI_Allgather(resultingValuesPerProc, rowsPerProc, MPI_DOUBLE, currentResultingValues,
                      rowsPerProc, MPI_DOUBLE, MPI_COMM_WORLD);

        // Выход из цикла происходит при достижении желаемой точности вычислений
    } while ((FindMatrixNorm(currentResultingValues, previousResultingValues, dimension) >= EPS) &&
             (iterationNumber < MAX_ITERATIONS));

    delete[] previousResultingValues;
    delete[] currentResultingValues;
}

void LinearJacobiMethod(double *coefficientMatrix, double *resultingValuesPerProc, double *freeMembersColumn,
                        int dimension, double EPS, int MAX_ITERATIONS)
{
    double *previousResultingValues = new double[dimension];
    double *currentResultingValues = new double[dimension];

    // На первом шаге метода Якоби за текущее решение берется столбец свободных членов
    for (int i = 0; i < dimension; i++)
    {
        currentResultingValues[i] = freeMembersColumn[i];
    }

    int iterationNumber = 0;
    do
    {
        iterationNumber++;
        // Своппинг currentResultingValues и previousResultingValues
        double *temp = previousResultingValues;
        previousResultingValues = currentResultingValues;
        currentResultingValues = temp;

        for (int i = 0; i < dimension; i++)
        {
            resultingValuesPerProc[i] = freeMembersColumn[i];
            for (int j = 0; j < i; j++)
            {
                resultingValuesPerProc[i] -= coefficientMatrix[i * dimension + j] * previousResultingValues[j];
            }
            for (int j = i + 1; j < dimension; j++)
            {
                resultingValuesPerProc[i] -= coefficientMatrix[i * dimension + j] * previousResultingValues[j];
            }
            // j Не может быть равным startIndex, так как на этом шаге текущий элемент coefficientMatrix
            // принадлежит главной диагонали, а значит по условиям метода на него нужно делить
            resultingValuesPerProc[i] /= coefficientMatrix[i * dimension + i];
        }

        // Сбор полученных на текущем шаге решений в единый массив currentResultingValues
        for (int i = 0; i < dimension; i++)
        {
            currentResultingValues[i] = resultingValuesPerProc[i];
        }

        // Выход из цикла происходит при достижении желаемой точности вычислений
    } while ((FindMatrixNorm(currentResultingValues, previousResultingValues, dimension) >= EPS) &&
             (iterationNumber < MAX_ITERATIONS));

    delete[] previousResultingValues;
    delete[] currentResultingValues;
}

double FindMatrixNorm(double *previousResultingValues, double *currentResultingValues, int dimension)
{
    double sum = 0.0;

    for (int i = 0; i < dimension; i++)
    {
        sum += abs(previousResultingValues[i] - currentResultingValues[i]);
    }

    return sum;
}

void PrintSolutions(double *solutions, int dimension)
{
    for (int i = 0; i < dimension; i++)
    {
        cout << setprecision(4) << solutions[i] << " ";
    }
    cout << endl;
}

void GenerateMatrix(double *&matrix, int dimension)
{
    matrix = new double[dimension * dimension];

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            // Элементы главной диагонали матрицы должны быть больше суммы других элементов строки
            // Иначе условие сходимости может не выполняться, и в таком случае метод Якоби применять нельзя
            matrix[dimension * i + j] = (i == j) ? (rand() % 30 + 10 * (dimension - 1))
                                                 : (rand() % 10 + 1);
        }
    }
}

void GenerateVector(double *&vector, int dimension)
{
    vector = new double[dimension];

    for (int i = 0; i < dimension; i++)
    {
        vector[i] = rand() % 30 + 1;
    }
}

void PrintLinearEquationSystem(double *coefficientMatrix, double *freeMembersColumn,
                               int dimension)
{
    cout << "Linear equation system" << endl;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            cout << setw(4) << setprecision(4);
            cout << coefficientMatrix[dimension * i + j];
            cout << "x" << j + 1;
            (j != dimension - 1) ? cout << "   + " : cout << " ";
        }
        cout << "= " << freeMembersColumn[i] << endl;
    }
}

bool AreSolutionsEqual(double *linearSolutions, double *parallelSolutions, int dimension, double EPS)
{
    for (int i = 0; i < dimension; i++)
    {
        if ((parallelSolutions[i] - linearSolutions[i]) > EPS)
        {
            return false;
        }
    }

    return true;
}
