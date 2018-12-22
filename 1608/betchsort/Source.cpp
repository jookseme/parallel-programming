#include <iostream>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

int* create_array(int size)
{
    int* array = new int[size];
    for (int i = 0; i < size; i++)
    {
        array[i] = rand() % 1000;
    }
    return array;
}

void compexch(int* a, int* b)
{
    if (b < a)
        swap(a, b);
}
int* get_part_of_array(int* arr, int start_pos, int end_pos)
{
    int* part_arr = new int[end_pos - start_pos];
    int j = 0;
    for (int i = start_pos; i < end_pos; i++)
    {
        part_arr[j] = arr[i];
        j++;
    }
    return part_arr;
}
int* merge_arrays(int* arr1, int* arr2, int block_size)
{
    int* merged_arr = new int[block_size * 2];
    for (int i = 0; i < block_size; i++)
    {
        merged_arr[i] = arr1[i];
        merged_arr[i + block_size] = arr2[i];
    }
    return merged_arr;
}
int* merge_arrays(int* arr1, int* arr2, int arr1_block_size, int arr2_block_size)
{
    int* merged_arr = new int[arr1_block_size + arr2_block_size];
    int j = 0;
    for (int i = 0; i < arr1_block_size; i++)
    {
        merged_arr[j] = arr1[i];
        j++;
    }
    for (int i = 0; i < arr2_block_size; i++)
    {
        merged_arr[j] = arr2[i];
        j++;
    }
    return merged_arr;
}

int getMax(int* arr, int n)
{
    int max = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}
void countSort(int* arr, int n, int exp)
{
    int *output = new int(n); // output array 
    int i, count[10] = { 0 };

    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
    for (i = n - 1; i >= 0; i--)
    {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    for (i = 0; i < n; i++)
        arr[i] = output[i];
}
void radixsort(int arr[], int n)
{
    int m = getMax(arr, n);
    for (int exp = 1; m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}



void print_array(int* array, int length)
{
    cout << "A: ";
    for (int i = 0; i < length; i++)
    {
        cout << array[i] << " ";
    }
    cout << endl;
}
/* Схема:
На ROOT_NODE создается исходный массив, который рассылается по нодам, с размером блока ( размер_массива / количество_нод )
На ROOT_NODE также остается только массив с полученным размером блока.
Потом уже работает алгоритм чет-нечетной перестановки.
Узлы меняются друг с другом массивами.
Сортируются.
Младший нод оставляет себе меньшую часть. Старший большую. */

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);
    int procNum = 0, procRank = 0;
    double startTime, endTime;

    int *a = new int[n];
    int* recv_array;
    int* tmp_recv_array;
    int** send_array;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (procRank == 0)
    {
        /*if (n%procNum != 0)
        {
            cout << "N % ProcNum !=0";
            return -1;
        }*/
        a = create_array(n);
        print_array(a, n);
        int* tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = a[i];
        startTime = MPI_Wtime();
        radixsort(tmp, n);
        endTime = MPI_Wtime();
        print_array(tmp, n);
        cout << "Linear time: " << endTime - startTime;
        startTime = MPI_Wtime();
    }

    int nodes_count = procNum;

    if (n < nodes_count)
    {
        MPI_Finalize();
        cout << "ERROR!!! Array size should be greater or equal to nodes count! " << endl;
        return 0;
    }

    int current_node = procRank;

    if (nodes_count > n)
        nodes_count = n - 1;
    const int LAST = nodes_count - 1;
    const int ROOT = 0;
    int block_size = n / nodes_count;
    bool diff_slices;
    if (n % nodes_count != 0)
        diff_slices = true;
    else
        diff_slices = false;


    cout << "current node is " << current_node << endl;
    if (current_node != ROOT && current_node != LAST)
    {
        recv_array = new int[block_size];
        MPI_Recv(recv_array, block_size, MPI_INT, ROOT, 0, MPI_COMM_WORLD, &status);
    }
    /* Для LAST */
    else if (current_node == LAST && current_node != ROOT)
    {
        if (diff_slices)
        {
            recv_array = new int[block_size + n % nodes_count];
            MPI_Recv(recv_array, block_size + n % nodes_count, MPI_INT, ROOT, 0, MPI_COMM_WORLD, &status);
        }
        else
        {
            recv_array = new int[block_size];
            MPI_Recv(recv_array, block_size, MPI_INT, ROOT, 0, MPI_COMM_WORLD, &status);
        }
    }
    /* Для ROOT */
    else
    {
        /* Делим массив на блоки для отправки */
        send_array = new int*[nodes_count];
        for (int i = ROOT; i < LAST; i++)
            send_array[i] = new int[block_size];

        if (diff_slices)
            send_array[LAST] = new int[block_size + n % nodes_count];
        else
            send_array[LAST] = new int[block_size];

        for (int i = ROOT; i < LAST; i++)
            for (int j = 0; j < block_size; j++)
                send_array[i][j] = a[i*block_size + j];

        if (diff_slices)
            for (int j = 0; j < block_size + n % nodes_count; j++)
                send_array[LAST][j] = a[block_size*(LAST)+j];
        else
            for (int j = 0; j < block_size; j++)
                send_array[LAST][j] = a[(LAST)*block_size + j];
        /* ----- */

        if (nodes_count != 1)
        {
            for (int i = 1; i < LAST; i++)
            {
                MPI_Send(send_array[i], block_size, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            if (diff_slices)
                MPI_Send(send_array[LAST], block_size + n % nodes_count, MPI_INT, LAST, 0, MPI_COMM_WORLD);
            else
                MPI_Send(send_array[LAST], block_size, MPI_INT, LAST, 0, MPI_COMM_WORLD);
        }

        recv_array = new int[block_size];
        recv_array = send_array[0];

    }
    /* Ждем пока каждый нод получит свою долю */
    MPI_Barrier(MPI_COMM_WORLD);

    for (int j = 0; j < nodes_count; j++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (j == current_node)
            if (diff_slices && current_node == LAST)
            {
                cout << "node " << current_node << " count " << block_size + n % nodes_count << endl;
                for (int i = 0; i < block_size + n % nodes_count; i++)
                    cout << "node " << current_node << " array " << i << " is " << recv_array[i] << endl;
            }
            else
            {
                for (int i = 0; i < block_size; i++)
                    cout << "node " << current_node << " array " << i << " is " << recv_array[i] << endl;
            }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (nodes_count == 1)
        radixsort(recv_array, n);

    /* Алгоритм чет-нечетной перестановки */
    /* Обмен массивами между нодами, слияния, сортировка и отсечение */
    for (int n = 0; n < nodes_count; n++)
    {
        /* Четная итерация */
        if (n % 2 == 0)
        {
            for (int j = 0; j < LAST; j += 2)
            {
                if (current_node == j)
                {
                    /* Отправление соседней ноде своего массива */  
                    MPI_Send(recv_array, block_size, MPI_INT, j + 1, 0, MPI_COMM_WORLD);
                    /* Если есть блок с остатком и соседний нод последний */
                    /* Если это не так, то в текущем ноде всегда будет количество элементов без остатка */
                    if (diff_slices && j + 1 == LAST)
                    {
                        tmp_recv_array = new int[block_size + n % nodes_count];
                        MPI_Recv(tmp_recv_array, block_size + n % nodes_count, MPI_INT, j + 1, 0, MPI_COMM_WORLD, &status);
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size, block_size + n % nodes_count);
                        radixsort(merged_array, block_size + block_size + n % nodes_count);
                        recv_array = get_part_of_array(merged_array, 0, block_size);
                    }
                    else
                    {
                        tmp_recv_array = new int[block_size];
                        MPI_Recv(tmp_recv_array, block_size, MPI_INT, j + 1, 0, MPI_COMM_WORLD, &status);
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size);
                        radixsort(merged_array, block_size * 2);
                        recv_array = get_part_of_array(merged_array, 0, block_size);
                    }
                }
                else if (current_node == j + 1)
                {
                    tmp_recv_array = new int[block_size];
                    /* Если текущий нод последний и есть остаток, то отправляем соседу массив с размером блока с остатком */
                    if (diff_slices && current_node == LAST)
                        MPI_Send(recv_array, block_size + n % nodes_count, MPI_INT, j, 0, MPI_COMM_WORLD);
                    else
                        MPI_Send(recv_array, block_size, MPI_INT, j, 0, MPI_COMM_WORLD);

                    MPI_Recv(tmp_recv_array, block_size, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
                    if (diff_slices && current_node == LAST)
                    {
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size + n % nodes_count, block_size);
                        radixsort(merged_array, block_size * 2 + n % nodes_count);
                        recv_array = get_part_of_array(merged_array, block_size, 2 * block_size + n % nodes_count);
                    }
                    else
                    {
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size);
                        radixsort(merged_array, block_size * 2);
                        recv_array = get_part_of_array(merged_array, block_size, 2 * block_size);
                    }
                }
            }
        }
        /* Нечетная итерация */
        else
        {
            for (int j = 1; j < LAST; j += 2)
            {
                if (current_node == j)
                {
                    MPI_Send(recv_array, block_size, MPI_INT, j + 1, 0, MPI_COMM_WORLD);
                    if (diff_slices && j + 1 == LAST)
                    {
                        tmp_recv_array = new int[block_size + n % nodes_count];
                        MPI_Recv(tmp_recv_array, block_size + n % nodes_count, MPI_INT, j + 1, 0, MPI_ANY_TAG, &status);
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size, block_size + n % nodes_count);
                        radixsort(merged_array, block_size * 2 + n % nodes_count);
                        recv_array = get_part_of_array(merged_array, 0, block_size);
                    }
                    else
                    {
                        tmp_recv_array = new int[block_size];
                        MPI_Recv(tmp_recv_array, block_size, MPI_INT, j + 1, 0, MPI_COMM_WORLD, &status);
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size);
                        radixsort(merged_array, block_size * 2);
                        recv_array = get_part_of_array(merged_array, 0, block_size);
                    }
                }
                else if (current_node == j + 1)
                {
                    tmp_recv_array = new int[block_size];
                    if (diff_slices && current_node == LAST)
                        MPI_Send(recv_array, block_size + n % nodes_count, MPI_INT, j, 0, MPI_COMM_WORLD);
                    else
                        MPI_Send(recv_array, block_size, MPI_INT, j, 0, MPI_COMM_WORLD);

                    MPI_Recv(tmp_recv_array, block_size, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
                    if (diff_slices && current_node == LAST)
                    {
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size + n % nodes_count, block_size);
                        radixsort(merged_array, block_size * 2 + n % nodes_count);
                        recv_array = get_part_of_array(merged_array, block_size, 2 * block_size + n % nodes_count);
                    }
                    else
                    {
                        int* merged_array = merge_arrays(recv_array, tmp_recv_array, block_size);
                        radixsort(merged_array, block_size * 2);
                        recv_array = get_part_of_array(merged_array, block_size, 2 * block_size);
                    }
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int j = 0; j < nodes_count; j++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (j == current_node)
            if (diff_slices && current_node == LAST)
            {
                cout << " count " << block_size + n % nodes_count << endl;
                for (int i = 0; i < block_size + n % nodes_count; i++)
                    cout<<"node " << current_node << " Sorted array " << i<<" is " << recv_array[i] << endl;
            }
            else
            {
                for (int i = 0; i < block_size; i++)
                {
                    cout << "node " << current_node << " Sorted array " << i << " is " << recv_array[i] << endl;
                }
            }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (procRank == 0)
        cout << " Paralel time : " << MPI_Wtime() - startTime << endl << endl;
    MPI_Finalize();
    return 0;
}