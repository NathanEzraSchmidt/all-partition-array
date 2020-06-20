#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include "ap_array.h"

/*
619 903 1382 nate_ezra_1981
123
321171184
*/

void printa(int* arr, int n)
{
    for (int i = 0; i < n; ++i)
        printf("%d ", arr[i]);
    printf("\n");
}

void get_cumsum_arr(int* output, int* input, int len)
{
    int x = 0;
    for (int i = 0; i < len; ++i){
        output[i] = x;
        x += input[i];
    }
    output[len] = x;
}

void get_cumsum_arr_2(int* output, int* input, int len, int constant)
{
    int x = 0;
    for (int i = 0; i < len; ++i){
        output[i] = x + constant;
        x += input[i];
    }
    output[len] = x + constant;
}

void reverse_matrix(struct V* block)
{
    struct V* b = block;
    int cols = b -> cols + 2;
    int rows = b -> rows;
    int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD];
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            matrix[j][i] = b -> matrix[j][cols-1-i];
        }
    }

    for (int i = 0; i < rows; ++i)
        memcpy(b -> matrix[i], matrix[i], cols*sizeof(int));

}

int enumerate_combinations(int* output, int n, int k)
{
    int x = 0;
    int a[PARTITIONS_NUM];
    for (int i = 0; i < PARTITIONS_NUM; ++i)
        a[i] = -1;
    for (int i = 0; i < k; ++i)
        a[i] = i;
    int pos = k - 1;
    while (1) {
        memcpy(output, a, sizeof(a));
        output += PARTITIONS_NUM;
        ++x;
        while (a[pos] < n - 1) {
            ++a[pos];
            memcpy(output, a, sizeof(a));
            output += PARTITIONS_NUM;
            ++x;
        }
        while (a[pos] == n - k + pos) {
            --pos;
            if (pos == -1)
                return x;
        }
        ++a[pos];
        while (pos < k - 1) {
            ++pos;
            a[pos] = a[pos-1] + 1;
        }
    }
}

int enumerate_compositions(int* output, int pt_num)
{
    // compositions of pt_num usings 2s and 3s, stored in output and return len
    int* combs = malloc(2000000*sizeof(int));

    //output should be initialized to zeros
    // k = {2,3};
    // n = pt_num
    int P = PARTITIONS_NUM;
    int p = pt_num % 2;
    int total_len = 0;
    for (int i = 2-p; i <= pt_num / 3; i += 2) {
        int three = i;
        int x = (pt_num - 3*three) / 2;
        int two = x;
        int n = two + three;
        int len = enumerate_combinations(combs, n, three);

        for (int j = 0; j < len; ++j) {
            int* comb = combs + j*P;
            int composition[P];
            for (int jj = 0; jj < P; ++jj)
                composition[jj] = -1;
            for (int jj = 0; jj < n; ++jj)
                composition[jj] = 2;
            for (int jj = 0; jj < three; ++jj) {
                composition[comb[jj]] = 3;
            }

            memcpy(output+total_len*P, composition, P*sizeof(int));
            ++total_len;
        }
    }

    if (p == 0) {
        for (int i = 0; i < P; ++i)
            output[total_len*P+i] = -1;
        for (int i = 0; i < pt_num / 2; ++i) {
            output[total_len*P+i] = 2;
        }
    }

    free(combs);
    if (p == 0)
        return total_len + 1;
    else
        return total_len;

}

void shuffle(int* arr, int n)
{
    for (int i = 0; i < n; ++i){
        int r = i + rand() % (n - i);
        int x = arr[i];
        arr[i] = arr[r];
        arr[r] = x;
    }
}

void reorder_arr(int* input_arr, int* input_order, int row_len, int len)
{ // reorders input arr  using order of positions in input order and stores result in output

    // if input_order[n] = x, then output[n] = input_arr[x]

    int output[PARTITIONS_NUM*MAX_MOD];
    int n = row_len;

    for (int i = 0; i < len; ++i) {
        memcpy(output+i*n, input_arr+input_order[i]*n, sizeof(int)*row_len);
    }
    memcpy(input_arr, output, row_len*len*sizeof(int));
}


int get_max_len(int* arr, int len)
{
    int x = 0;
    for (int i = 0; i < len; ++i)
        if (arr[i] > x)
            x = arr[i];
    return x;
}

int get_max_elem(int* arr, int len)
{
    int x = -4;
    for (int i = 0; i < len; ++i)
        if (arr[i] > x)
            x = arr[i];
    return x;
}

void get_reverse_sorted_positions(int* output, int* input, int len)
{
    /*
    input is array of ints in any order
    if rs = reverse(sorted(arr)), then output[n] = x, then rs[n] = input[x]
    in other words, output is the sequential positions in input of the corresponding values of rs
    */

    int new_input[PARTITIONS_NUM];
    memcpy(new_input, input, len*sizeof(int));


    int output_len = 0;

    while (1) {

        int x = get_max_elem(new_input, len);

        if (x == -1)
            break;

        for (int i = 0; i < len; ++i) {
            if (new_input[i] == x) {
                output[output_len] = i;
                ++output_len;
                new_input[i] = -1;
            }
        }
    }

}

int get_min_len(int* arr, int len)
{
    int x = arr[0];
    for (int i = 1; i < len; ++i) {
        if (!arr[i])
            break;
        if (arr[i] < x)
            x = arr[i];
    }
    return x;
}

int get_min_elem(int* arr, int len)
{
    int x = arr[0];
    for (int i = 1; i < len; ++i) {
        if (arr[i] < x)
            x = arr[i];
    }
    return x;
}

int contains(int* arr, int len, int elem)
{
    for (int i = 0; i < len; ++i) {
        if (!arr[i])
            break;
        if (arr[i] == elem)
            return 1;
    }
    return 0;
}

int check_unique_pts_allowed(int* arr_1, int* arr_2, int len)
{
    // returns 1 if arr_1 and arr_2 have 0's in different places

    for (int i = 0; i < len; ++i) {
        if (arr_1[i] == 0 && arr_2[i] == 0)
            return 0;
    }
    return 1;
}

int get_parts_num(int* arr, int len)
{
    for (int i = 0; i < len; ++i)
        if (!arr[i])
            return i;
    return len;
}

int get_partition_num(int x)
{
    int partition_nums[17] = {0, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231};
    return partition_nums[x];
}

void init_array(int* p, int val, int len)
{
    for (int i = 0; i < len; ++i)
        p[i] = val;
}

int input_partition_list(int partition_list[MAX_MOD * PARTITIONS_NUM], int n, int max_seg, int min_seg, int max_parts, int min_parts)
{
    /*
    returns num of partitions in list
    */

    for (int i = 0; i < MAX_MOD * 128; ++i)
        partition_list[i] = 0;
    int partition_list_temp[MAX_MOD * 128];
    char filename[64];
    sprintf(filename, "partition lists/%d.bin", n);
    FILE* fp = fopen(filename, "rb");
    fread((void*) partition_list_temp, sizeof(int), MAX_MOD * 128, fp);
    fclose(fp);
    int partition_num = get_partition_num(n);

    int j = 0;

    for (int i = 0; i < partition_num; ++i) {
        int* x = &partition_list_temp[i*(n+1)];
        if (x[0] < min_seg)
            continue;

        int seg, parts;

        for (int k = 0; k < n; ++k) {
            seg = x[n-1-k];
            if (seg) {
                parts = n - k;
                break;
            }
        }

        if (seg > max_seg || parts > max_parts || parts < min_parts)
            continue;

        memcpy(&partition_list[MAX_MOD*j], &partition_list_temp[i*(n+1)], (n+1)*sizeof(int));
        ++j;
    }
    return j;
}

int output_partition_list(int partition_list[MAX_MOD*PARTITIONS_NUM], int n,
                          int max_seg, int min_seg, int max_parts, int min_parts,
                          char* output_filename)
{
    /*
    returns num of partitions in list
    */
    int P = PARTITIONS_NUM;
    for (int i = 0; i < MAX_MOD * P; ++i)
        partition_list[i] = 0;
    int partition_list_temp[MAX_MOD * P];
    char filename[64];
    sprintf(filename, "partition lists/%d.bin", n);
    FILE* fp = fopen(filename, "rb");
    fread((void*) partition_list_temp, sizeof(int), MAX_MOD * P, fp);
    fclose(fp);

    printa(partition_list_temp, 16);
    int partition_num = get_partition_num(n);

    printf("pn %d\n", partition_num);

    int j = 0;

    for (int i = 0; i < partition_num; ++i) {
        int* x = &partition_list_temp[i*(n+1)];
        if (x[0] < min_seg)
            continue;

        int seg, parts;

        for (int k = 0; k < n; ++k) {
            seg = x[n-1-k];
            if (seg) {
                parts = n - k;
                break;
            }
        }

        if (seg > max_seg || parts > max_parts || parts < min_parts)
            continue;

        memcpy(&partition_list[MAX_MOD*j], &partition_list_temp[i*(n+1)], (n+1)*sizeof(int));
        ++j;
    }

    fp = fopen(output_filename, "w");
    for (int i = 0; i < j; ++i) {
        for (int k = 0; k < 16; ++k) {
            fprintf(fp, "%d ", partition_list[MAX_MOD*i+k]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return j;
}

int input_partition_list_txt(int partition_list[MAX_MOD*PARTITIONS_NUM], char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (fp == 0) {
        printf("cannot find file: %s", filename);
        fclose(fp);
        return -1;
    }
    int p = 0;
    while (1) {
        int x;
        fscanf(fp, "%d", &x);
        if (x == 17)
            break;
        partition_list[p] = x;
        ++p;
    }
    fclose(fp);
    return p / 16;
}

int input_partition_list_2(int partition_list[MAX_MOD * PARTITIONS_NUM], int n, int max_seg, int min_seg, int max_parts, int min_parts)
{
    /*
    returns num of partitions in list
    */

    for (int i = 0; i < MAX_MOD * 128; ++i)
        partition_list[i] = 0;
    int partition_list_temp[MAX_MOD * 128];
    char filename[64];
    sprintf(filename, "partition lists/%d.txt", n);
    FILE* fp = fopen(filename, "r");

    int len = 0;
    int x;
    while (fscanf(fp, "%d", &x) != EOF) {
        partition_list_temp[len] = x;
        ++len;
    }
    fclose(fp);
    int partition_num = len / 16; //get_partition_num(n);
    int j = 0;

    for (int i = 0; i < partition_num; ++i) {
        int* x = &partition_list_temp[i*(16)];
        if (x[0] < min_seg)
            continue;

        int seg, parts;

        for (int k = 0; k < n; ++k) {
            seg = x[n-1-k];
            if (seg) {
                parts = n - k;
                break;
            }
        }

        if (seg > max_seg || parts > max_parts || parts < min_parts)
            continue;

        memcpy(&partition_list[MAX_MOD*j], &partition_list_temp[i*(16)], (16)*sizeof(int));
        ++j;
    }
    return j;
}


void output_all_partition_array(char* name,
                                int matrix[MAX_MOD][PARTITIONS_NUM * MAX_MOD],
                                int result[PARTITIONS_NUM][MAX_MOD],
                                int rep_tried[PARTITIONS_NUM][MAX_MOD],
                                int partitions_used,
                                int rows)

{
    char* vals[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D"};
    FILE* fp = fopen(name, "w");

    int* pointer_array[MAX_MOD];
    for (int i = 0; i < MAX_MOD; ++i)
        pointer_array[i] = &matrix[i][1];

    int max_len_list[PARTITIONS_NUM];
    init_array(max_len_list, 0, PARTITIONS_NUM);

    for (int i = 0; i < partitions_used; ++i)
        max_len_list[i] = get_max_len(result[i], MAX_MOD);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < partitions_used; ++j) {
            int len = result[j][i];
            int max_len = max_len_list[j];
            if (len && rep_tried[j][i])
                --pointer_array[i];
            for (int k = 0; k < len; ++k) {
                fprintf(fp, "%s", vals[*pointer_array[i]]);
                ++pointer_array[i];
            }
            for (int k = len; k < max_len + 1; ++k)
                fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


void get_matrix_with_rep(int new_matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                int result[PARTITIONS_NUM][MAX_MOD],
                                int rep_tried[PARTITIONS_NUM][MAX_MOD],
                                int partitions_used,
                                int rows, int* offset_arr)

{
    for (int i = 0; i < MAX_MOD; ++i)
        for (int j = offset_arr[i]; j < PARTITIONS_NUM*MAX_MOD; ++j)
        new_matrix[i][j] = -1;
    int* pointer_array[MAX_MOD];
    for (int i = 0; i < MAX_MOD; ++i)
        pointer_array[i] = &matrix[i][1];

    int h;

    for (int i = 0; i < rows; ++i) {
            h = 0;
        for (int j = 0; j < partitions_used; ++j) {
            int len = result[j][i];
            if (len && rep_tried[j][i])
                --pointer_array[i];
            for (int k = 0; k < len; ++k) {
                new_matrix[i][h+offset_arr[i]] = *pointer_array[i];
                ++h;
                ++pointer_array[i];
            }
        }
        offset_arr[i] += h;
    }
}

void get_position_sums(int sums[MAX_MOD], int concatenated_result[PARTITIONS_NUM][MAX_MOD], int rows, int len)
{
    for (int i = 0; i < MAX_MOD; ++i)
        sums[i] = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < len; ++j) {
            sums[i] += concatenated_result[j][i];
        }
    }
}

int get_last_lens(int concatenated_result[PARTITIONS_NUM][MAX_MOD], int concatenated_result_len, int rows,
                   int matrix_with_rep[MAX_MOD][PARTITIONS_NUM*MAX_MOD])
{

    int found = 0;
    int sums[MAX_MOD];

    get_position_sums(sums, concatenated_result, rows, concatenated_result_len);
    printa(sums, 4);

    for (int i = 0; i < rows; ++i){
        int x = sums[i];
        int len = 0;
        while (1) {
            int n = matrix_with_rep[i][x];
//            printf("n = %d\n", n);
//            getchar();
            if (n == -1)
                break;
            found = 1;
            ++len;
            ++x;
        }
        concatenated_result[concatenated_result_len][i] = len;
    }
    return found;
}

int get_concatenated_result(int result[PARTITIONS_NUM][MAX_MOD], int concatenated_result[PARTITIONS_NUM][MAX_MOD],
                            int result_len, int inc, int rows,
                            int matrix_with_rep[MAX_MOD][PARTITIONS_NUM*MAX_MOD])
{
    /*
    assumes that

    */

    int P = PARTITIONS_NUM;
    int M = MAX_MOD;
    int concatenated_result_len = 0;

    for (int i = 0; i < PARTITIONS_NUM; ++i)
        for (int j = 0; j < MAX_MOD; ++j)
            concatenated_result[i][j] = 0;

    for (int i = 0; i < rows; ++i) {
        int len = 0;
        concatenated_result_len = 0;
        for (int j = 0; j < result_len; ++j) {
            int x = result[j][i];
            len += x;
            if ((j+1) % inc == 0 || j == result_len - 1) {
                concatenated_result[concatenated_result_len][i] = len;
                len = 0;
                concatenated_result_len++;
            }
        }
    }
    return concatenated_result_len + get_last_lens(concatenated_result, concatenated_result_len, rows, matrix_with_rep);
}


void get_result_rows_cumsum(struct V* block)
{
    for (int i = 0; i < block -> rows; ++i) {
        int* ptr = block -> result_rows[i];
        int x = 0;
        for (int j = 0; j < PARTITIONS_NUM; ++j) {
            block -> result_rows_cumsum[i][j] = x;
            x += block -> result_rows[i][j];

        }
    }
}

void get_result_rows(struct V* block)
{
    for (int i = 0; i < PARTITIONS_NUM; ++i) {
        for (int j = 0; j < MAX_MOD; ++j) {
            block -> result_rows[j][i] = block -> result[i][j];
        }
    }

    get_result_rows_cumsum(block);
}


void shuffle_pt_list(struct V* block)
{

    int N = PARTITIONS_NUM;
    int arr[N];
    for (int i = 0; i < N; ++i) {
        arr[i] = i;
    }
    shuffle(arr, block -> pt_num);

    for (int i = 0; i < block -> pt_num; ++i) {
        int k = arr[i];
        block -> pt_order[i] = k;

        int* pt = block -> pt_list_original + MAX_MOD * k;

        for (int j = 0; j < MAX_MOD; ++j) {
            block -> pt_list[i*MAX_MOD + j] = pt[j];
        }
    }
}

void shuffle_pt_list_thread(struct V* block, int* arr)
{

    int N = PARTITIONS_NUM;

    for (int i = 0; i < block -> pt_num; ++i) {
        int k = arr[i];
        block -> pt_order[i] = k;


        int* pt = block -> pt_list_original + MAX_MOD * k;

        for (int j = 0; j < MAX_MOD; ++j) {
            block -> pt_list[i*MAX_MOD + j] = pt[j];
        }
    }
}

void shuffle_pt_list_partial(struct V* block, int len, int thread_id)
{
    int N = PARTITIONS_NUM;
    int arr[N];
    for (int i = 0; i < N; ++i) {
        arr[i] = i;
    }
    for (int i = 0; i < thread_id; ++i)
        shuffle(arr, len);

    EnterCriticalSection(&critical);

    for (int i = 0; i < len; ++i) {
        int k = arr[i];
        block -> pt_order[i] = k;
        int* pt = block -> pt_list_original + MAX_MOD * k;

        for (int j = 0; j < MAX_MOD; ++j) {
            block -> pt_list[i*MAX_MOD + j] = pt[j];
        }
    }
    LeaveCriticalSection(&critical);
}

void set_variables_global(struct V* global, char* output_filename)
{
//    global -> sequences = malloc(12*12*sizeof(int));
//    global -> sequences_possible = malloc(12*PARTITIONS_NUM*sizeof(int));

    int i;
    for (i = 0; i < MAX_MOD; ++i)
        for (int j = 0; j < PARTITIONS_NUM * MAX_MOD; ++j)
            global -> matrix[i][j] = -1;

    char filename[64];
    sprintf(filename, "input files/");
    FILE* fp;

    open_file:
        printf("matrix input file name: ");
        scanf("%s", &filename[12]);
        fp = fopen(filename, "r");
        if (fp == 0) {
            printf("cannot find file: %s\n", filename);
            goto open_file;
        }

    sprintf(output_filename, "output files/%s", &filename[12]);
    fscanf(fp, "%d%d%d%d%d%d%d%d", &global -> m, &global -> max_seg_len, &global -> min_seg_len,
           &global -> max_parts_num, &global -> min_parts_num,
           &global -> rows, &global -> cols, &global -> pts_used);

    int n = global -> m;

    for (i = 0; i < global -> rows; ++i){
        global -> matrix[i][0] = n;
        int j = 0;
        while (1) {
            ++j;
            int x;
            fscanf(fp, "%d", &x);
            if (x == 17) {
                global -> matrix[i][j] = n;
                break;
            }
            else
                global -> matrix[i][j] = x;
        }
        global -> matrix[i][global -> cols + 1] = n;
    }
    fclose(fp);


    while (1) {

    printf("partition list file name: ");

    sprintf(filename, "partition lists txt/");

    scanf("%s", filename+20);
    puts(filename);
    global -> pt_num = input_partition_list_txt(global -> pt_list, filename);//"12a.txt");
//    printf("pt num= %d\n", global -> pt_num);
    if (global -> pt_num != -1)
        break;

    }

    memcpy(global -> pt_list_original, global -> pt_list, sizeof(global -> pt_list_original));

    init_array(global -> rev_pt_list, 0, MAX_MOD*128);

    for (i = 0; i < global -> pt_num; ++i) {
        int k = i*MAX_MOD;
        int h = 0;
        for (int j = global -> m - 1; j >= 0; --j) {
            int x = global -> pt_list[k+j];
            if (x) {
                global -> rev_pt_list[k+h] = x;
                ++h;
            }
        }
    }

    if (global -> pts_used == 0)
        global -> pts_used = global -> pt_num;

    for (i = 0; i < MAX_MOD; ++i)
        global -> pointer_array[i] = &global -> matrix[i][1];

    init_array(global -> v_pos, -1, PARTITIONS_NUM);

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> allowed[i], 0, MAX_MOD);

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> allowed[i], 1, n);

    /*
    Partitions
    */

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> result[i], 0, MAX_MOD);

    init_array(global -> pts_allowed, 1, PARTITIONS_NUM);
    init_array(global -> p_pos, -1, PARTITIONS_NUM);
    init_array(global -> which_pt, -1, PARTITIONS_NUM);


    /*
    REP
    */

    init_array(global -> rep_count, 0, MAX_MOD);


    for (i = 0; i < MAX_MOD; ++i)
        global -> pc_count[i] = 0;
    for (i = 0; i < global -> rows; ++i) {
        for (int j = 1; j < global -> cols + 1; ++j) {
            ++global -> pc_count[global -> matrix[i][j]];
        }
    }

    for (i = 0; i < global -> m; ++i) {
        if (global -> pc_count[i] >= global -> pts_used) {
            global -> max_rep_count[i] = 0;
        }
        else {
            global -> max_rep_count[i] = global -> pts_used - global -> pc_count[i];
        }
    }

    init_array(global -> rep_tried[0], 0, PARTITIONS_NUM*MAX_MOD);

    global -> c = -1;
    global -> tries = 0;
    global -> name_int = 0;
    global -> n = global -> m;
}

void set_variables(struct V* global, int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD], int rows, int cols,
                   int matrix_start, int pt_num, int pts_used, int rep_amount, int m,
                   int* pt_list)
{
    char filename[64];
    sprintf(filename, "sequences/");

    printf("sequence file: ");
    scanf("%s",filename+10);
    puts(filename);

    FILE* fp = fopen(filename, "r");



    fscanf(fp, "%d%d%d%d%d", &global -> group_size,
                            &global -> seq_num,
                            &global -> seq_len,
                            &global -> valid_selection_string_num,
                            &global -> valid_selection_string_len);

    global -> sequences = malloc(global -> seq_num * global -> seq_len * sizeof(int));
    global -> sequences_possible = malloc(global -> seq_num *PARTITIONS_NUM*sizeof(int));
    global -> valid_selection_strings = malloc(global -> valid_selection_string_num
                                               * global -> valid_selection_string_len
                                               * sizeof(int)); // size is num of different valid selections strins
                                                            // times len of each string. contains the values of the
                                                          // of the selection strings
    int k = 0;

    for (int i = 0; i < global -> seq_num; ++i) {
        for (int j = 0; j < global -> seq_len; ++j) {
            fscanf(fp, "%d", &global -> sequences[k]);// = (i+j)%12;
            ++k;
        }
    }

    k = 0;



    for (int i = 0; i < global -> valid_selection_string_num; ++i) {
        for (int j = 0; j < global -> valid_selection_string_len; ++j) {
            fscanf(fp, "%d", &global -> valid_selection_strings[k]);// = (j*6+i)%12;//(i+j)%12;
            ++k;
        }
    }

    fclose(fp);

    k = 0;

    int i;
    for (i = 0; i < MAX_MOD; ++i)
        for (int j = 0; j < PARTITIONS_NUM * MAX_MOD; ++j)
            global -> matrix[i][j] = -1;

    global -> rows = rows;
    global -> cols = cols;
    global -> pts_used = pts_used;
    global -> m = m;
    global -> pt_num = pt_num;
    global -> n = m;
    int n = m;

    for (i = 0; i < global -> rows; ++i){
        for (int j = 1; j < global -> cols + 1; ++j){
            global -> matrix[i][j] = matrix[i][j+matrix_start];
        }
        global -> matrix[i][global -> cols + 1] = n;
        global -> matrix[i][0] = n;
    }

    for (i = 0; i < pt_num*MAX_MOD; ++i) {
        global -> pt_list[i] = pt_list[i];
        global -> pt_list_original[i] = pt_list[i];
    }

    init_array(global -> rev_pt_list, 0, MAX_MOD*128);

    for (i = 0; i < global -> pt_num; ++i) {
        int k = i*MAX_MOD;
        int h = 0;
        for (int j = global -> m - 1; j >= 0; --j) {
            int x = global -> pt_list[k+j];
            if (x) {
                global -> rev_pt_list[k+h] = x;
                ++h;
            }
        }
    }

    for (i = 0; i < MAX_MOD; ++i)
        global -> pointer_array[i] = &global -> matrix[i][1];

    init_array(global -> v_pos, -1, PARTITIONS_NUM);

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> allowed[i], 0, MAX_MOD);

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> allowed[i], 1, n);

    /*
    Partitions
    */

    for (i = 0; i < PARTITIONS_NUM; ++i)
        init_array(global -> result[i], 0, MAX_MOD);

    init_array(global -> pts_allowed, 1, PARTITIONS_NUM);
    init_array(global -> p_pos, -1, PARTITIONS_NUM);
    init_array(global -> which_pt, -1, PARTITIONS_NUM);


    /*
    REP
    */

    init_array(global -> rep_count, 0, MAX_MOD);


    for (i = 0; i < MAX_MOD; ++i)
        global -> pc_count[i] = 0;
    for (i = 0; i < global -> rows; ++i) {
        for (int j = 1; j < global -> cols + 1; ++j) {
            ++global -> pc_count[global -> matrix[i][j]];
        }
    }

    for (i = 0; i < global -> m; ++i) {
        global -> max_rep_count[i] = rep_amount;
    }

    init_array(global -> rep_tried[0], 0, PARTITIONS_NUM*MAX_MOD);

    global -> c = -1;
    global -> tries = 0;
    global -> name_int = 0;
}

int get_next_partition(struct V* global, char* output_filename, int pp)

{

    int i;
    int j;
    int tt = 0;
    if (pp > 0)
        goto main_backward;
    main_forward: // start getting new aggregate
        ++global -> c;
        if (global -> c == global -> pts_used) {
            return 1;
            ++pp;
            goto main_backward;
        }
        do {
            ++global -> which_pt[global -> c];
        } while (!global -> pts_allowed[global -> which_pt[global -> c]]);

        if (global -> which_pt[global -> c] == global -> pt_num) {
            global -> which_pt[global -> c] = -1;
            goto main_backward;
        }
        global -> _pt = &global -> pt_list[global -> which_pt[global -> c]*MAX_MOD];
        goto _forward;

    main_backward:
        --global -> c;
        if (global -> c == -1) {
            return 0;
        }
        global -> pts_allowed[global -> which_pt[global -> c]] = 1;
        global -> _pt = &global -> pt_list[global -> which_pt[global -> c]*MAX_MOD];
        for (i = 0; i < global -> rows; ++i) {
            global -> pointer_array[i] -= global -> result[global -> c][i];
        }
        goto backward;

    _forward:
        ++global -> p_pos[global -> c];
        ++global -> v_pos[global -> c];
        global -> len = global -> _pt[global -> p_pos[global -> c]];
        if (!global -> len) {
            for (i = 0; i < global -> m; ++i) {
                global -> pointer_array[i] += global -> result[global -> c][i];
            }
            global -> pts_allowed[global -> which_pt[global -> c]] = 0;
            goto main_forward;
        }
        if (global -> v_pos[global -> c] == global -> rows) {
            goto backward;
        }
        if (global -> result[global -> c][global -> v_pos[global -> c]]) {
            --global -> p_pos[global -> c];
            goto _forward;
        }
        global -> p = global -> pointer_array[global -> v_pos[global -> c]];

        for (i = 0; i < global -> len; ++i) {
            if (!global -> allowed[global -> c][global -> p[i]]) {
                for (j = 0; j < i; ++j)
                    global -> allowed[global -> c][global -> p[j]] = 1;
                --global -> p_pos[global -> c];
                goto _forward_rep;
            }
            global -> allowed[global -> c][global -> p[i]] = 0;
        }
        global -> result[global -> c][global -> v_pos[global -> c]] = global -> len;
        if (global -> _pt[global -> p_pos[global -> c] + 1] && global -> _pt[global -> p_pos[global -> c] + 1] != global -> len)
            global -> v_pos[global -> c] = -1;
        goto _forward;

    _forward_rep:
        ++global -> p_pos[global -> c];
        global -> p = global -> pointer_array[global -> v_pos[global -> c]] - 1;

        if (global -> rep_count[*(global -> p)] == global -> max_rep_count[*(global -> p)]) {
            --global -> p_pos[global -> c];
            goto _forward;
        }
        for (i = 0; i < global -> len; ++i) {
            if (!global -> allowed[global -> c][global -> p[i]]) {
                for (j = 0; j < i; ++j)
                    global -> allowed[global -> c][global -> p[j]] = 1;
                --global -> p_pos[global -> c];
                goto _forward;
            }
            global -> allowed[global -> c][global -> p[i]] = 0;
        }

        global -> rep_tried[global -> c][global -> v_pos[global -> c]] = 1;
        global -> result[global -> c][global -> v_pos[global -> c]] = global -> len;
        --global -> pointer_array[global -> v_pos[global -> c]];
        ++global -> rep_count[*(global -> p)];

        if (global -> _pt[global -> p_pos[global -> c] + 1] && global -> _pt[global -> p_pos[global -> c] + 1] != global -> len)
            global -> v_pos[global -> c] = -1;
        goto _forward;

    backward:
        --global -> p_pos[global -> c];
        if (global -> p_pos[global -> c] == -1) {
            global -> v_pos[global -> c] = -1;
            --global -> c;
            goto main_forward;
        }
        --global -> v_pos[global -> c];
        global -> len = global -> _pt[global -> p_pos[global -> c]];
        if (global -> result[global -> c][global -> v_pos[global -> c]] != global -> len) {
            ++global -> p_pos[global -> c];
            goto backward;
        }
        global -> result[global -> c][global -> v_pos[global -> c]] = 0;
        global -> p = global -> pointer_array[global -> v_pos[global -> c]];
        for (i = 0; i < global -> len; ++i)
            global -> allowed[global -> c][global -> p[i]] = 1;
        --global -> p_pos[global -> c];
        if (global -> rep_tried[global -> c][global -> v_pos[global -> c]]) {
            global -> rep_tried[global -> c][global -> v_pos[global -> c]] = 0;
            --global -> rep_count[*(global -> p)];
            ++global -> pointer_array[global -> v_pos[global -> c]];
            goto _forward;
        }
        else
            goto _forward_rep;
}
int func(struct V* global)
{

    if (global -> group_size == 17)
        return 1;

    int* r = global -> result[global -> c];
    int c = global -> c;

    int M = MAX_MOD;
    int P = PARTITIONS_NUM;
    int rows = global -> rows;

    int matrix[16][16];
    for (int i = 0; i < rows; ++i){
        int* ptr = global -> pointer_array[i];
        ptr -= r[i];
        for (int j = 0; j < r[i]; ++j) {
            matrix[i][j] = *ptr;
            ++ptr;
        }
        matrix[i][r[i]] = global -> n;
    }

    int seq_len = global -> seq_len;
    int pos[M];
    init_array(pos, 0, M);
    int seq_num = global -> seq_num;
    int* seq = global -> sequences;
    int* seq_possible = global -> sequences_possible + c*seq_num;
    init_array(seq_possible, 0, seq_num);
    int total_found = 0;

    if (global -> group_size < 2) {
        for (int i = 0; i < seq_num; ++i) {
            int* s = seq + i * seq_len;
            init_array(pos, 0, M);
            total_found = 0;

            for (int j = 0; j < seq_len; ++j){
                int elem = s[j];
                int found = 0;
                for (int k = 0; k < rows; ++k) {
                    int val = matrix[k][pos[k]];
                    if (val == elem) {
                        ++total_found;
                        if (total_found == seq_len)
                            seq_possible[i] = 1;

                        ++pos[k];
                        found = 1;
                        //global -> valid_choice = i;
                        break;
                    }
                }
                if (found == 0 || total_found == seq_len)
                    break;
            }
        }
    }

    else {
        int pc_counts[M*M];
        init_array(pc_counts, 0, M*M);
        int size = global -> group_size;
        int matrix_positions[M];

        init_array(matrix_positions, 1, M);

        for (int i = 0; i < rows; ++i)  {
            for (int j = 0; j < c; ++j) {
                int x = global -> result[j][i];
                int start = matrix_positions[i];
                for (int k = start; k < start + x; ++k) {
                    int elem = global -> matrix[i][k];
                    ++pc_counts[(i/global -> group_size)*M+elem];

                }
                matrix_positions[i] += x;
            }
        }

        int pc_counts_mins[M];
        init_array(pc_counts_mins, 0, M);

        for (int i = 0; i < rows; ++i) {
            pc_counts_mins[i] = get_min_elem(pc_counts+M*i, global -> m);
        }

        for (int i = 0; i < seq_num; ++i) {
            int pc_counts_new[M*M];
            memcpy(pc_counts_new, pc_counts, M*M*sizeof(int));
            int* s = seq + i * seq_len;
            init_array(pos, 0, M);
            total_found = 0;

            for (int j = 0; j < seq_len; ++j){
                int elem = s[j];
                int found = 0;
                for (int k = 0; k < rows; ++k) {
                    int val = matrix[k][pos[k]];
                    if (val == elem) {
                        int z = k / global -> group_size;

                        ++pc_counts_new[z*M+elem];

                        if (pc_counts_new[z*M+elem] - pc_counts_mins[z] > 1) {// && z > 0) {
                            found = 0;
                            break;
                        }

                        ++total_found;

                        if (total_found == seq_len)
                            seq_possible[i] = 1;

                        ++pos[k];
                        found = 1;
                        //global -> valid_choice = i;
                        break;
                    }
                }
                if (found == 0 || total_found == seq_len) // found seq;
                    break;
            }
        }
    }

    int found = 0;
//
//    printa(global -> valid_selection_strings, 50);
//    printa(global -> valid_selection_strings+72, 50);
//    printa(global -> valid_selection_strings+144,50);

    int len = global -> valid_selection_string_len;
//    printf("len = %d\n", len);
//    printa(global -> valid_selection_strings + 72, 72);
//    puts("R");
//    printf("num %d\n",  global -> valid_selection_string_num);
//    getchar();
//    getchar();

    for (int i = 0; i < global -> valid_selection_string_num; ++i) { // number of different valid strings
//        printf("i = %d\n", i);
//        printf("%d\n",  global -> valid_selection_string_num);
        int* v = global -> valid_selection_strings + i * len;

//        printa(v, 72);
//        getchar();
//        getchar();

//        printa(v, 16);
//        getchar();
        int found = 1;
        for (int j = 0; j < c+1; ++j) {
            int* s = global -> sequences_possible + j*seq_num;
            if (s[v[j]] == 0) {
                found = 0;
                break;
            }
        }
        if (found == 1) {
            global -> valid_choice = i;
//            puts("found");
//            printf("c i = %d %d\n", c);
            return 1;
        }
    }
//    puts("not found");
//    printf("c = %d\n", c);
    return 0;
}
int get_next_partition_with_reps(struct V* g, char* output_filename, int pp, int reps_allowed, int (*func)(struct V*))

{
    int i;
    int j;
    int tt = 0;
    if (pp > 0)
        goto main_backward;

    int valid;

    int rep_count = 0;

    int pt_num_original = g -> pt_num;
    g -> pt_num *= 2;
    int xx = g -> pt_num / 4;
    g -> pt_num -= xx;
    memcpy(g -> pt_list + pt_num_original*MAX_MOD, g -> pt_list, sizeof(int)*MAX_MOD*pt_num_original);

    main_forward: // start getting new aggregate
        ++g -> c;
        if (g -> c == g -> pts_used) {
            --g -> c;
            if (func(g) == 1)
                return 1;
            else
                goto main_backward;
            ++pp;
//            goto main_backward;
        }
        do {
            ++g -> which_pt[g -> c];
        } while (!g -> pts_allowed[g -> which_pt[g -> c]]);

        if (g -> which_pt[g -> c] == g -> pt_num) {
            g -> which_pt[g -> c] = -1;
            goto main_backward;
        }
        if (g -> which_pt[g -> c] >= pt_num_original && rep_count >= reps_allowed) {
            g -> which_pt[g -> c] = -1;
            goto main_backward;
        }

        g -> _pt = &g -> pt_list[g -> which_pt[g -> c]*MAX_MOD];
        goto _forward;

    main_backward:
        --g -> c;
        if (g -> c == -1) {
            //puts("result not possible");
            return 0;
        }
        g -> pts_allowed[g -> which_pt[g -> c]] = 1;
        g -> _pt = &g -> pt_list[g -> which_pt[g -> c]*MAX_MOD];
        if (g -> which_pt[g -> c] >= pt_num_original)
            --rep_count;
        for (i = 0; i < g -> rows; ++i) {
            g -> pointer_array[i] -= g -> result[g -> c][i];
        }
        goto backward;

    _forward:
        ++g -> p_pos[g -> c];
        ++g -> v_pos[g -> c];
        g -> len = g -> _pt[g -> p_pos[g -> c]];
        if (!g -> len) {
            for (i = 0; i < g -> rows; ++i) {
                g -> pointer_array[i] += g -> result[g -> c][i];
            }
            g -> pts_allowed[g -> which_pt[g -> c]] = 0;
            if (g -> which_pt[g -> c] >= pt_num_original)
                ++rep_count;
            valid = func(g);

            if (valid == 0) {
                ++g -> c;
                goto main_backward;
            }
            goto main_forward;
        }
        if (g -> v_pos[g -> c] == g -> rows) {
            goto backward;
        }
        if (g -> result[g -> c][g -> v_pos[g -> c]]) {
            --g -> p_pos[g -> c];
            goto _forward;
        }
        g -> p = g -> pointer_array[g -> v_pos[g -> c]];

        for (i = 0; i < g -> len; ++i) {
            if (!g -> allowed[g -> c][g -> p[i]]) {
                for (j = 0; j < i; ++j)
                    g -> allowed[g -> c][g -> p[j]] = 1;
                --g -> p_pos[g -> c];
                goto _forward_rep;
            }
            g -> allowed[g -> c][g -> p[i]] = 0;
        }
        g -> result[g -> c][g -> v_pos[g -> c]] = g -> len;
        if (g -> _pt[g -> p_pos[g -> c] + 1] && g -> _pt[g -> p_pos[g -> c] + 1] != g -> len)
            g -> v_pos[g -> c] = -1;
        goto _forward;

    _forward_rep:
        ++g -> p_pos[g -> c];
        g -> p = g -> pointer_array[g -> v_pos[g -> c]] - 1;

        if (g -> rep_count[*(g -> p)] == g -> max_rep_count[*(g -> p)]) {
            --g -> p_pos[g -> c];
            goto _forward;
        }
        for (i = 0; i < g -> len; ++i) {
            if (!g -> allowed[g -> c][g -> p[i]]) {
                for (j = 0; j < i; ++j)
                    g -> allowed[g -> c][g -> p[j]] = 1;
                --g -> p_pos[g -> c];
                goto _forward;
            }
            g -> allowed[g -> c][g -> p[i]] = 0;
        }

        g -> rep_tried[g -> c][g -> v_pos[g -> c]] = 1;
        g -> result[g -> c][g -> v_pos[g -> c]] = g -> len;
        --g -> pointer_array[g -> v_pos[g -> c]];
        ++g -> rep_count[*(g -> p)];

        if (g -> _pt[g -> p_pos[g -> c] + 1] && g -> _pt[g -> p_pos[g -> c] + 1] != g -> len)
            g -> v_pos[g -> c] = -1;
        goto _forward;

    backward:
        --g -> p_pos[g -> c];
        if (g -> p_pos[g -> c] == -1) {
            g -> v_pos[g -> c] = -1;
            --g -> c;
            goto main_forward;
        }
        --g -> v_pos[g -> c];
        g -> len = g -> _pt[g -> p_pos[g -> c]];
        if (g -> result[g -> c][g -> v_pos[g -> c]] != g -> len) {
            ++g -> p_pos[g -> c];
            goto backward;
        }
        g -> result[g -> c][g -> v_pos[g -> c]] = 0;
        g -> p = g -> pointer_array[g -> v_pos[g -> c]];
        for (i = 0; i < g -> len; ++i)
            g -> allowed[g -> c][g -> p[i]] = 1;
        --g -> p_pos[g -> c];
        if (g -> rep_tried[g -> c][g -> v_pos[g -> c]]) {
            g -> rep_tried[g -> c][g -> v_pos[g -> c]] = 0;
            --g -> rep_count[*(g -> p)];
            ++g -> pointer_array[g -> v_pos[g -> c]];
            goto _forward;
        }
        else
            goto _forward_rep;
}

void output_blocks(struct V blocks[], int n, int start_int)
{
    for (int i = 0; i < n; ++i) {

        char filename[64];
        sprintf(filename, "output files/%d.txt", i+start_int);
        output_all_partition_array(filename,
        blocks[i].matrix, blocks[i].result, blocks[i].rep_tried, blocks[i].pts_used, blocks[i].rows);
    }
}

int get_all_block_partitions(struct V block)
{
    int c = 0;
    char* output_filename;
    int first_pt = 1;
    while (1) {
        int x = get_next_partition(&block, output_filename, !first_pt);
        first_pt = 0;
        ++c;
        printf("c = %d\n", c);
        if (x == 0)
            break;
    }
    return c;
}

void revert_pts_allowed(struct V* block)
{
    /*
    in block, pts_allowed list is based on random order of pts
    converts pts_allowed to correspond to original ordering of pts (using block -> pt_order);
    */
    int pts_allowed[PARTITIONS_NUM];
    int n = block -> pt_num;
    init_array(pts_allowed, 1, PARTITIONS_NUM);

    for (int i = 0; i < n; ++i) {
        int x = block -> pt_order[i];
        int b = block -> pts_allowed[i];
        pts_allowed[x] = b;
    }
    memcpy(block -> pts_allowed, pts_allowed, sizeof(block -> pts_allowed));
}

void get_matrix_with_rep_from_blocks(struct V* blocks[], int block_num, int new_matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                     int extra_ints)
{
//    int new_matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD];
    // extra_ints is either 0 or 1

    /*

    assuming last block doesnt have any reps and completes the original input matrix
    */
    int offset_arr[16];
    for (int i = 0; i < 16; ++i)
        offset_arr[i] = 0;

    for (int i = 0; i < block_num; ++i) {
            struct V* block = blocks[i];

                EnterCriticalSection(&critical);
                printf("\npt list: ");
                printa(block -> pt_list, 16);


                LeaveCriticalSection(&critical);


        if (i == block_num - 1 && extra_ints){
            for (int j = 0; j < block -> rows; ++j) {
                for (int k = 1; k < block -> cols + 1; ++k) {
                    new_matrix[j][k-1+offset_arr[j]] = block -> matrix[j][k];
                }
            }
        }
        else {
//            struct V* block = blocks[i];

            get_matrix_with_rep(new_matrix,
                                    block -> matrix,
                                    block -> result,
                                    block -> rep_tried,
                                    block -> pts_used,
                                    block -> rows,
                                    offset_arr);
        }
    }
}

int get_result_from_blocks(int result[PARTITIONS_NUM][MAX_MOD], struct V* blocks[], int block_num)
{
    int n = 0;
    int k = 0;
    for (int i = 0; i < block_num; ++i) {
        int x = blocks[i] -> pts_used;
        for (int j = 0; j < x; ++j) {
            memcpy(result[k], blocks[i] -> result[j], MAX_MOD*sizeof(int));
            ++k;
        }
    }
    return k;
}

int partition_block(struct V* block, int* ptr, int output_start)
{

    /*
    convert pts_used to decimal number and store in result
    result is list of pts_allowed
    */

    int P = PARTITIONS_NUM;
    int first_pt = 1;
    char* output_filename;
    int x = 1;
    int c = 0;

    while (x) {
        int x = get_next_partition(block, output_filename, !first_pt);
        if (x == 0)
            break;
        first_pt = 0;

        struct V output_block;
        memcpy(&output_block, block, sizeof(output_block));
        struct V output_block_arr[] = {output_block};

        //output_blocks(output_block_arr, 1, output_start);

        memcpy(ptr+c*PARTITIONS_NUM, block -> pts_allowed, sizeof(block -> pts_allowed));
        ++output_start;
        ++c;
    }
    return output_start;
}

int get_valid_path(int* result, int* input[], int* input_lens, int input_len, int pt_num)
{
    /*
    input: 2d binary array, though effectively 3d, [P][1000][P]
    input_lens: values for 2nd dimension of input, example if input_lens[x] = y, then input[x][y][P]
    input_len len of input, if = x then , input = [x][1000][P]
    */
    int P = PARTITIONS_NUM;

    int pts_allowed[P];

    for (int i = 0; i < P; ++i){
        result[i] = 0;
        pts_allowed[i] = 1;
    }

    int c = 0; // position in input
    int d = 0; // position in input[c]

    int n = pt_num;

    while (1) {
        int* ptr = input[c] + d*P;

        if (input_lens[c] == 0)
            return 0;

        int valid = check_unique_pts_allowed(ptr, pts_allowed, n);
        if (valid == 1) {
            for (int i = 0; i < n; ++i) {
                if (ptr[i] == 0)
                    pts_allowed[i] = 0;
            }
            result[c] = d;
            ++c;
            if (c == input_len)
                return 1;
            else {
                d = 0;
                continue;
            }
        }
        else {
            ++d;
            while (d == input_lens[c]) {
                result[c] = 0;
                --c;
                if (c == -1)
                    return 0;
                ptr = input[c] + result[c]*P;
                for (int i = 0; i < n; ++i) {
                    if (ptr[i] == 0) {
                        pts_allowed[i] = 1;
                    }
                }
                ++result[c];
                d = result[c];
            }
        }
    }
    return -1;
}

int get_all_valid_paths(int* result_arr, int* pts_allowed_arr, int* input[], int* input_lens, int input_len, int pt_num)
{
    /*
    input: 2d binary array, though effectively 3d, [P][1000][P]
    input_lens: values for 2nd dimension of input, example if input_lens[x] = y, then input[x][y][P]
    input_len len of input, if = x then , input = [x][1000][P]
    */
    int P = PARTITIONS_NUM;

    int pts_allowed[P];

    int result[P];

    for (int i = 0; i < P; ++i){
        result[i] = 0;
        pts_allowed[i] = 1;
    }

    int c = 0; // position in input
    int d = 0; // position in input[c]

    int n = pt_num;

    int result_len = 0;
    int* ptr;
    int valid;
    while (1)
    {
        if (c < input_len) {
            ptr = input[c] + d*P;
            valid = check_unique_pts_allowed(ptr, pts_allowed, n);
        }
        else {
            valid = 0;
        }

        if (valid == 1) {
            for (int i = 0; i < n; ++i) {
                if (ptr[i] == 0)
                    pts_allowed[i] = 0;
            }
            result[c] = d;
            ++c;
            if (c == input_len) {
                memcpy(result_arr+result_len*P, result, sizeof(result));
                memcpy(pts_allowed_arr+result_len*P, pts_allowed, sizeof(result));
                ++result_len;

                //return 1;



            }
            else {
                d = 0;
                continue;
            }
        }
        else {
            ++d;
            while (d == input_lens[c]) {
                result[c] = 0;
                --c;
                if (c == -1)
                    return result_len;
                ptr = input[c] + result[c]*P;
                for (int i = 0; i < n; ++i) {
                    if (ptr[i] == 0) {
                        pts_allowed[i] = 1;
                    }
                }
                ++result[c];
                d = result[c];
            }
        }
    }
    return result_len;
}

void get_pts_used_arr(int* pts_used_arr, struct V blocks[], int block_num, int total_pts_used, int inc)
{

    int t = 0;
    for (int i = 0; i < block_num; ++i)
        t += blocks[i].pts_used;

    int t1 = total_pts_used - t;

    int x = t / inc;

    for (int i = 0; i < x; ++i)
        pts_used_arr[i] = inc;
    if (x * inc < t) {
        pts_used_arr[x] = t - x * inc;
        pts_used_arr[x+1] = t1;
    }
    else
        pts_used_arr[x] = t1;
}

int get_all_pt_array_with_reps(struct V* global_input, struct V* output_block)
{
    struct V global;
    struct V blocks[3];
    struct V block_0;
    char output_filename[64];

    memcpy(&global, global_input, sizeof(global));

    int n = global.pt_num;

    set_variables(&blocks[0], global.matrix,
                  global.rows, global.cols,
                   0, global.pt_num,global.pt_num,global.max_rep_count[0], global.m, global.pt_list);

    memcpy(&block_0, &blocks[0], sizeof(block_0));

    int block_num = 1;

    int pp = 0;
    int first_pt = 1;

    int tries = 0;

    /*
    PART 2
    */
    int P = PARTITIONS_NUM;
    pp = 0;
    int pts_allowed_counts[P];
    init_array(pts_allowed_counts, 0, P);
    puts("ordering partitions...");

    int first_run = 100;
    int pts_missing = global.pt_num * .25;

    for (int i = 0; i < first_run; ++i){
        pp = 0;
        printf("i %d\n", i);
        if (pp == 0) {
            ++tries;
            memcpy(&blocks[0], &block_0, sizeof(blocks[0]));
            shuffle_pt_list(&blocks[0]);
            first_pt = 1;
        }
        blocks[0].pts_used -= pts_missing;
        int x = get_next_partition(&blocks[0], output_filename, 0);
        revert_pts_allowed(&blocks[0]);
        for (int j = 0; j < blocks[0].pt_num; ++j) {
            pts_allowed_counts[j] += blocks[0].pts_allowed[j];
        }
    }

    int positions[P];
    get_reverse_sorted_positions(positions, pts_allowed_counts, blocks[0].pt_num);


    reorder_arr(block_0.pt_list, positions, MAX_MOD, blocks[0].pt_num);

    pp = 0;
    int M = MAX_MOD;
    reverse_matrix(&block_0);

    int reps_allowed = 0;

    while (pp < block_num){
        if (pp == 0) {
            ++tries;
            printf("tries = %d\n", tries);
            memcpy(&blocks[0], &block_0, sizeof(blocks[0]));
            first_pt = 1;
        }

        int x = get_next_partition_with_reps(&blocks[0], output_filename, 0, reps_allowed, func);
        if (x == 1) {
            ++pp;
            if (pp == block_num) {
                break;
            }
            if (pp == 1)  {
                revert_pts_allowed(&blocks[0]);
            }
            first_pt = 1;
            for (int i = 0; i < PARTITIONS_NUM; ++i)
                blocks[pp].pts_allowed[i] = blocks[pp-1].pts_allowed[i];
        }
        else {
            --pp;
            first_pt = 0;
        }
    }
    memcpy(output_block, &blocks[0], sizeof(blocks[0]));
    return 1;
    //output_blocks(blocks, block_num, start);
}

void get_pt_list_ordering(struct V* global_input, struct V* output_block)
{

    /*
    block_0 is initially set then saved as reference, blocks[0] is partitioned, at start of each sample
    reset blocks[0] to block_0
    */
    struct V global;
    struct V blocks[3];
    struct V block_0;
    char output_filename[64];

    memcpy(&global, global_input, sizeof(global));

    int n = global.pt_num;

    set_variables(&blocks[0], global.matrix,
                  global.rows, global.cols,
                   0, global.pt_num,global.pt_num,global.max_rep_count[0], global.m, global.pt_list);

    memcpy(&block_0, &blocks[0], sizeof(block_0));

    int block_num = 1;

    int pp = 0;
    int first_pt = 1;

    /*
    PART 2
    */

    int P = PARTITIONS_NUM;
    pp = 0;
    int pts_allowed_counts[P];
    init_array(pts_allowed_counts, 0, P);
    puts("ordering partitions...");

    int samples = 100;
//    int pts_missing = global.pt_num * .25;
    int pts_missing = global.pts_used * .25;

    printf("missing %d\n", pts_missing);

    for (int i = 0; i < samples; ++i){
        pp = 0;
        printf("i = %d\n", i);
        if (pp == 0) {
            memcpy(&blocks[0], &block_0, sizeof(blocks[0]));
            shuffle_pt_list(&blocks[0]);
            first_pt = 1;
        }
        blocks[0].pts_used -= pts_missing;
        int x = get_next_partition(&blocks[0], output_filename, 0);
        revert_pts_allowed(&blocks[0]);
        for (int j = 0; j < blocks[0].pt_num; ++j) {
            pts_allowed_counts[j] += blocks[0].pts_allowed[j];
        }
    }

    int positions[P];
    get_reverse_sorted_positions(positions, pts_allowed_counts, blocks[0].pt_num);

    reorder_arr(block_0.pt_list, positions, MAX_MOD, blocks[0].pt_num);

    memcpy(output_block, &block_0, sizeof(block_0));
//    memcpy(output_block -> pt_list, block_0.pt_list, block_0.pt_list);

}

DWORD WINAPI one_thread(void* params)
{
    struct Params* p = (struct Params*) params;
    struct V global_input; memcpy(&global_input, &p -> global_input, sizeof(global_input));

    struct V output_block; memcpy(&output_block, &p -> output_block, sizeof(global_input));

    int* output_int = p -> output_int;
    int thread_id;
    EnterCriticalSection(&critical);
    int* t_id = p -> thread_id;
    ++(*t_id);
    thread_id = *t_id;
    printf("thread id %d\n", thread_id);
    LeaveCriticalSection(&critical);

    output_block.pt_num = global_input.pt_num;

    int len = output_block.pts_used / 3;

    memcpy(output_block.pt_list_original, output_block.pt_list, sizeof(output_block.pt_list));

    shuffle_pt_list_partial(&output_block, 4, thread_id);

    char output_filename[64];

    int x = get_next_partition(&output_block, output_filename, 0);

    struct V blocks[] = {output_block};

    output_blocks(blocks, 1, thread_id);

    return 0;
}


DWORD WINAPI get_all_pt_array_thread(void* params)
{
    int P = PARTITIONS_NUM;

    struct V global;
    struct V blocks[3];
    struct V block_0;
    char output_filename[64];

    struct Params* p = (struct Params*) params;

    int thread_id;
    EnterCriticalSection(&critical);
    int* t_id = p -> thread_id;
    ++(*t_id);
    thread_id = *t_id;
    printf("thread id %d\n", thread_id);
    LeaveCriticalSection(&critical);
    int* shuffled_arr = p -> shuffled_arr + thread_id*P*100;

    memcpy(&global, &p -> global_input, sizeof(global));

    int n = global.pt_num;

    set_variables(&blocks[0], global.matrix,
                  global.rows, global.cols,
                   0, global.pt_num,global.pts_used, global.max_rep_count[0], global.m, global.pt_list);

    memcpy(&block_0, &blocks[0], sizeof(block_0));

    int block_num = 1;

    int pp = 0;
    int first_pt = 1;

    int tries = 0;

    /*
    PART 2
    */

    int pts_allowed_counts[P];
    init_array(pts_allowed_counts, 0, P);
    puts("ordering partitions...");

    int first_run = 100;
    int pts_missing = global.pts_used * .25;

    for (int i = 0; i < first_run; ++i){
        pp = 0;
        if (pp == 0) {
            ++tries;
            memcpy(&blocks[0], &block_0, sizeof(blocks[0]));
            shuffle_pt_list_thread(&blocks[0], shuffled_arr+P*i);
            first_pt = 1;
        }
        blocks[0].pts_used -= pts_missing;
        int x = get_next_partition(&blocks[0], output_filename, 0);
        revert_pts_allowed(&blocks[0]);
        for (int j = 0; j < blocks[0].pt_num; ++j) {
            pts_allowed_counts[j] += blocks[0].pts_allowed[j];
        }
    }
    puts("done");
    int positions[P];
    get_reverse_sorted_positions(positions, pts_allowed_counts, blocks[0].pt_num);

    reorder_arr(block_0.pt_list, positions, MAX_MOD, blocks[0].pt_num);

    pp = 0;
    int M = MAX_MOD;
    if (thread_id % 2 == 1) {
        reverse_matrix(&block_0);
    }

    int reps_allowed = 0;
    int found = 0;

    while (pp < block_num){
        if (pp == 0) {
            ++tries;
            printf("tries = %d\n", tries);
            memcpy(&blocks[0], &block_0, sizeof(blocks[0]));
            first_pt = 1;
        }

        int x = get_next_partition_with_reps(&blocks[0], output_filename, 0, reps_allowed, func);
        printf("choice = %d\n", blocks[0].valid_choice);
        if (x == 0 && pp == 0) {
            puts("result not possible");
            break;
        }
        if (x == 1) {
            found = 1;
            ++pp;
            if (pp == block_num) {
                break;
            }
            if (pp == 1)  {
                revert_pts_allowed(&blocks[0]);
            }
            first_pt = 1;
            for (int i = 0; i < PARTITIONS_NUM; ++i)
                blocks[pp].pts_allowed[i] = blocks[pp-1].pts_allowed[i];
        }
        else {
            --pp;
            first_pt = 0;
        }
    }
    if (found == 1)
        output_blocks(blocks, block_num, thread_id);
    free(blocks[0].sequences);
    free(blocks[0].sequences_possible);
    free(blocks[0].valid_selection_strings);

    return 0;
}

int shuffle_once()
{
    InitializeCriticalSection(&critical);
    srand(time(NULL));
    int P = PARTITIONS_NUM;
    int M = MAX_MOD;
    struct V global_input;
    struct V output_block;
    char output_filename[64];
    set_variables_global(&global_input, output_filename);

    get_pt_list_ordering(&global_input, &output_block);

    struct Params p;
    int output_int = 0;
    int thread_id = 0;
    memcpy(&p.global_input, &global_input, sizeof(global_input));
    memcpy(&p.output_block, &output_block, sizeof(output_block));
    p.output_int = &output_int;
    p.thread_id = &thread_id;

    HANDLE threads[100];

    printf("thread num = ");
    int thread_num;
    scanf("%d", &thread_num);

    for (int i = 0; i < thread_num; ++i) {
        threads[i] = CreateThread(0, 0, one_thread, &p, 0, 0);
    }

    for (int i = 0; i  < thread_num; ++i) {
            WaitForSingleObject(threads[i], INFINITE);
    }
    return 0;

}

//
//int func(struct V* global)
//{
//
////    return 1;
//
////    puts("YO");
//    int* r = global -> result[global -> c];
//    int c = global -> c;
//    int M = MAX_MOD;
//    int P = PARTITIONS_NUM;
//    int rows = global -> rows;
//
//    int matrix[16][16];
//    for (int i = 0; i < rows; ++i){
//        int* ptr = global -> pointer_array[i];
//        ptr -= r[i];
//        for (int j = 0; j < r[i]; ++j) {
//            matrix[i][j] = *ptr;
//            ++ptr;
//        }
//        matrix[i][r[i]] = global -> n;
//    }
//
////    printf("c = %d\n", global -> c);
////    getchar();
////    printa(global -> result[global -> c], 16);
////    for (int i = 0; i < global -> rows; ++i) {
////        printa(matrix[i], 12);
////    }
//
//    int seq_len = global -> seq_len;
//    int pos[M];
//    init_array(pos, 0, M);
//    int seq_num = global -> seq_num;
//    int* seq = global -> sequences;
//    int* seq_possible = global -> sequences_possible + c*seq_num;
//    init_array(seq_possible, 0, seq_num);
//    int total_found = 0;
//
//
//
//
//    int group_ints_allowed[M][P];
//
//    for (int ii = 0; ii < M; ++ii)
//        memcpy(group_ints_allowed[ii], global -> group_ints_allowed[ii], P*sizeof(int));
//
//    int total_int_count[M];
//    memcpy(total_int_count, global -> total_int_count, M*sizeof(int));
//
//    for (int i = 0; i < seq_num; ++i) {
//        int* s = seq + i * seq_len;
//        init_array(pos, 0, M);
//        total_found = 0;
//
//        for (int j = 0; j < seq_len; ++j){
//            int elem = s[j];
//            int found = 0;
//            for (int k = 0; k < rows; ++k) {
//                int val = matrix[k][pos[k]];
//                if (val == elem) {
//                    ++total_found;
//                    if (total_found == seq_len)
//                        seq_possible[i] = 1;
//
//                    ++pos[k];
//                    found = 1;
//                    //global -> valid_choice = i;
//                    break;
//                }
//            }
//            if (found == 0 || total_found == seq_len)
//                break;
//        }
//    }
//    int found = 0;
////
////    printa(global -> valid_selection_strings, 50);
////    printa(global -> valid_selection_strings+72, 50);
////    printa(global -> valid_selection_strings+144,50);
//
//    int len = global -> valid_selection_string_len;
////    printf("len = %d\n", len);
////    printa(global -> valid_selection_strings + 72, 72);
////    puts("R");
////    printf("num %d\n",  global -> valid_selection_string_num);
////    getchar();
////    getchar();
//
//    for (int i = 0; i < global -> valid_selection_string_num; ++i) { // number of different valid strings
////        printf("i = %d\n", i);
////        printf("%d\n",  global -> valid_selection_string_num);
//        int* v = global -> valid_selection_strings + i * len;
//
////        printa(v, 72);
////        getchar();
////        getchar();
//
////        printa(v, 16);
////        getchar();
//        int found = 1;
//        for (int j = 0; j < c+1; ++j) {
//            int* s = global -> sequences_possible + j*seq_num;
//            if (s[v[j]] == 0) {
//                found = 0;
//                break;
//            }
//        }
//        if (found == 1) {
//            global -> valid_choice = i;
////            puts("found");
////            printf("c i = %d %d\n", c);
//            return 1;
//        }
//
//    }
////    puts("not found");
////    printf("c = %d\n", c);
//    return 0;
//}
//void set_variables_old(struct V* global, int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD], int rows, int cols,
//                   int matrix_start, int pt_num, int pts_used, int rep_amount, int m,
//                   int* pt_list)
//{
////    global -> group_size = 2;
////    global -> group_combinatoriality = 0;
////    for (int i = 0; i < MAX_MOD; ++i)
////        for (int j = 0; j < PARTITIONS_NUM; ++j)
////            global -> group_ints_allowed[i][j] = 1;
////    for (int i = 0; i < MAX_MOD; ++i)
////        global -> total_int_count[i] = 0;
//
//
//    global -> sequences = malloc(12*12*sizeof(int));
//    global -> sequences_possible = malloc(12*PARTITIONS_NUM*sizeof(int));
//    global -> seq_len = 12;
//    global -> seq_num = 12;
//
//    global -> valid_selection_string_len = 72; // the len of each valid string of valid selections
//                                            // typically = pts_used, sinch agg has it's own sequence
//    global -> valid_selection_string_num = 12;
//
//    global -> valid_selection_strings = malloc(global -> valid_selection_string_num
//                                               * global -> valid_selection_string_len
//                                               * sizeof(int)); // size is num of different valid selections strins
//                                                            // times len of each string. contains the values of the
//                                                          // of the selection strings
//    int k = 0;
//    for (int i = 0; i < 12; ++i) {
//        for (int j = 0; j < 72; ++j) {
//            global -> valid_selection_strings[k] = (j*6+i)%12;//(i+j)%12;
//            ++k;
//        }
//    }
//    for (int i = 0; i < 12; ++i) {
//
//    printa(global -> valid_selection_strings+i*72, 72);
//    puts("////////////");
//    }
////    printa(global -> valid_selection_strings+72, 72);
//
//
//
//    getchar();
//    getchar();
//
//
//
//    k = 0;
//    for (int i = 0; i < 12; ++i) {
//        for (int j = 0; j < 12; ++j) {
////            printf("%d\n", (i+j)%12);
//            global -> sequences[k] = (i+j)%12;
//            ++k;
//        }
//    }
//
////    for (int i = 0 ; i < 12; ++i)
////        printa(global -> sequences+i*12, 12);
//
//
//
//    int i;
//    for (i = 0; i < MAX_MOD; ++i)
//        for (int j = 0; j < PARTITIONS_NUM * MAX_MOD; ++j)
//            global -> matrix[i][j] = -1;
//
//    global -> rows = rows;
//    global -> cols = cols;
//    global -> pts_used = pts_used;
//    global -> m = m;
//    global -> pt_num = pt_num;
//    global -> n = m;
//    int n = m;
//
//    for (i = 0; i < global -> rows; ++i){
//        for (int j = 1; j < global -> cols + 1; ++j){
//            global -> matrix[i][j] = matrix[i][j+matrix_start];
//        }
//        global -> matrix[i][global -> cols + 1] = n;
//        global -> matrix[i][0] = n;
//    }
//
//    for (i = 0; i < pt_num*MAX_MOD; ++i) {
//        global -> pt_list[i] = pt_list[i];
//        global -> pt_list_original[i] = pt_list[i];
//    }
//
//    init_array(global -> rev_pt_list, 0, MAX_MOD*128);
//
//    for (i = 0; i < global -> pt_num; ++i) {
//        int k = i*MAX_MOD;
//        int h = 0;
//        for (int j = global -> m - 1; j >= 0; --j) {
//            int x = global -> pt_list[k+j];
//            if (x) {
//                global -> rev_pt_list[k+h] = x;
//                ++h;
//            }
//        }
//    }
//
//    for (i = 0; i < MAX_MOD; ++i)
//        global -> pointer_array[i] = &global -> matrix[i][1];
//
//    init_array(global -> v_pos, -1, PARTITIONS_NUM);
//
//    for (i = 0; i < PARTITIONS_NUM; ++i)
//        init_array(global -> allowed[i], 0, MAX_MOD);
//
//    for (i = 0; i < PARTITIONS_NUM; ++i)
//        init_array(global -> allowed[i], 1, n);
//
//    /*
//    Partitions
//    */
//
//    for (i = 0; i < PARTITIONS_NUM; ++i)
//        init_array(global -> result[i], 0, MAX_MOD);
//
//    init_array(global -> pts_allowed, 1, PARTITIONS_NUM);
//    init_array(global -> p_pos, -1, PARTITIONS_NUM);
//    init_array(global -> which_pt, -1, PARTITIONS_NUM);
//
//
//    /*
//    REP
//    */
//
//    init_array(global -> rep_count, 0, MAX_MOD);
//
//
//    for (i = 0; i < MAX_MOD; ++i)
//        global -> pc_count[i] = 0;
//    for (i = 0; i < global -> rows; ++i) {
//        for (int j = 1; j < global -> cols + 1; ++j) {
//            ++global -> pc_count[global -> matrix[i][j]];
//        }
//    }
//
//    for (i = 0; i < global -> m; ++i) {
//        global -> max_rep_count[i] = rep_amount;
//    }
//
//    init_array(global -> rep_tried[0], 0, PARTITIONS_NUM*MAX_MOD);
//
//    global -> c = -1;
//    global -> tries = 0;
//    global -> name_int = 0;
//}
