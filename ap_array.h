#ifndef _AP_ARRAY
#define _AP_ARRAY
#define MAX_MOD 16
#define PARTITIONS_NUM 128

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>

CRITICAL_SECTION critical;

struct V {

    int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD];
    int n;
    int m;
    int max_seg_len, min_seg_len, max_parts_num, min_parts_num;

    int pt_list[MAX_MOD * PARTITIONS_NUM];
    int pt_list_original[MAX_MOD * PARTITIONS_NUM];

    int pt_order[PARTITIONS_NUM];

    int rows, cols;
    int pt_num;
    int rev_pt_list[MAX_MOD * PARTITIONS_NUM];
    int pts_used;

    int* pointer_array[MAX_MOD];
    int v_pos[PARTITIONS_NUM];

    int allowed[PARTITIONS_NUM][MAX_MOD];
    int result[PARTITIONS_NUM][MAX_MOD];
    int result_rows[MAX_MOD][PARTITIONS_NUM];
    int result_rows_cumsum[MAX_MOD][PARTITIONS_NUM];
    int pts_allowed[PARTITIONS_NUM];
    int p_pos[PARTITIONS_NUM];
    int which_pt[PARTITIONS_NUM];
    int* _pt;

    /*
    REP
    */

    int rep_count[MAX_MOD];
    int max_rep_count[MAX_MOD];
    int pc_count[MAX_MOD];
    int rep_tried[PARTITIONS_NUM][MAX_MOD];
    int c;
    int len, j;
    int* p;
    int tries;
    int name_int;

    int* sequences;  //{0,1,2,...,9,10,11,1,2,3,...,10,11,0,2,3,...}
    int* sequences_possible; // size P * len(sequences)//{{0,1,0,...,}
    int seq_num;
    int seq_len;

    int* valid_selection_strings; // below two leangth of array multiple is len
    int valid_selection_string_len; // length each each individual valid str (= pts_used)
    int valid_selection_string_num; // different valid strings

    int valid_choice ;

//    int group_combinatoriality;
//    int group_ints_allowed[MAX_MOD][PARTITIONS_NUM];
//    int total_int_count[MAX_MOD];

    int group_size;
};

struct Params {
    struct V global_input;
    struct V output_block;
    int* output_int;
    int* thread_id;
    int* shuffled_arr;
};

void printa(int* arr, int n);

void get_cumsum_arr(int* output, int* input, int len);

void get_cumsum_arr_2(int* output, int* input, int len, int constant);

void reverse_matrix(struct V* block);

void shuffle(int* arr, int n);

int enumerate_combinations(int* ouput, int n, int k);

int enumerate_compositions(int* output, int pt_num);

void reorder_arr(int* input_arr, int* input_order, int row_len, int len);

int get_max_len(int* arr, int len);

void get_reverse_sorted_positions(int* output, int* input, int len);

int get_min_len(int* arr, int len);
int get_min_elem(int* arr, int len);

int contains(int* arr, int len, int elem);

int check_unique_pts_allowed(int* arr_1, int* arr_2, int len);

int get_parts_num(int* arr, int len);

int get_partition_num(int x);

//int get_partition_list(int partition_list[MAX_MOD * PARTITIONS_NUM], int partition_num, int max_seg, int min_seg, int seg_required, int max_parts, int min_parts);

void init_array(int* p, int val, int len);

int input_partition_list(int partition_list[MAX_MOD * PARTITIONS_NUM], int n, int max_seg, int min_seg, int max_parts, int min_parts);

int output_partition_list(int partition_list[MAX_MOD*PARTITIONS_NUM], int n, int max_seg, int min_seg, int max_parts, int min_parts, char* output_filename);

int input_partition_list_txt(int partition_list[MAX_MOD*PARTITIONS_NUM], char* filename);


int input_partition_list_2(int partition_list[MAX_MOD * PARTITIONS_NUM], int n, int max_seg, int min_seg, int max_parts, int min_parts);

void output_all_partition_array(char* name,
                                int matrix[MAX_MOD][PARTITIONS_NUM * MAX_MOD],
                                int result[PARTITIONS_NUM][MAX_MOD],
                                int rep_tried[PARTITIONS_NUM][MAX_MOD],
                                int partitions_used,
                                int rows);

void get_matrix_with_rep(int new_matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                int result[PARTITIONS_NUM][MAX_MOD],
                                int rep_tried[PARTITIONS_NUM][MAX_MOD],
                                int partitions_used,
                                int rows, int* offset_arr);

void get_position_sums(int sums[MAX_MOD], int concatenated_result[PARTITIONS_NUM][MAX_MOD], int rows, int len);

int get_last_lens(int concatenated_result[PARTITIONS_NUM][MAX_MOD], int concatenated_result_len, int rows,
                   int matrix_with_rep[MAX_MOD][PARTITIONS_NUM*MAX_MOD]);

int get_concatenated_result(int result[PARTITIONS_NUM][MAX_MOD], int concatenated_result[PARTITIONS_NUM][MAX_MOD],
                            int result_len, int inc, int rows,
                            int matrix_with_rep[MAX_MOD][PARTITIONS_NUM*MAX_MOD]);

void get_result_rows_cumsum(struct V* block);

void get_result_rows(struct V* block);

void shuffle_pt_list(struct V* block);

void shuffle_pt_list_thread(struct V* block, int* arr);

void shuffle_pt_list_partial(struct V* block, int len, int thread_id);

void set_variables_global(struct V* global, char* output_filename);
void set_variables(struct V* global, int matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD], int rows, int cols,
                   int matrix_start, int pt_num, int pts_used, int rep_amount, int m,
                   int* pt_list);
int get_next_partition(struct V* global, char* output_filename, int pp);

int get_next_partition_with_reps(struct V* global, char* output_filename, int pp, int reps_allowed, int (*func)(struct V*));
void output_blocks(struct V blocks[], int n, int start_int);
int get_all_block_partitions(struct V block);

void revert_pts_allowed(struct V* block);

void get_matrix_with_rep_from_blocks(struct V* blocks[], int block_num, int new_matrix[MAX_MOD][PARTITIONS_NUM*MAX_MOD],
                                     int extra_ints);
int get_result_from_blocks(int result[PARTITIONS_NUM][MAX_MOD], struct V* blocks[], int block_num);

int get_valid_path(int* result, int* input[], int* input_lens, int input_len, int pt_num);

int get_all_valid_paths(int* result_arr, int* pts_allowed_arr, int* input[], int* input_lens, int input_len, int pt_num);
void get_pts_used_arr(int* pts_used_arr, struct V blocks[], int block_num, int total_pts_used, int inc);

int get_all_pt_array_with_reps(struct V* global_input, struct V* output_block);

void get_pt_list_ordering(struct V* global_input, struct V* output_block);

DWORD WINAPI one_thread(void* params);

DWORD WINAPI get_all_pt_array_thread(void* params);

int shuffle_once();

#endif
