/*
** author: Nathan Ezra Schmidt
** email: nathanezraschmidt@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "ap_array.h"

int main()
{

    /*
    int output_partition_list(int partition_list[MAX_MOD*PARTITIONS_NUM], int n,
                          int max_seg, int min_seg, int max_parts, int min_parts,
                          char* output_filename)
    */



//    int pp[MAX_MOD*PARTITIONS_NUM];
//
//    char yo[64];
//    sprintf(yo, "partition lists txt/4.txt");
//
//    output_partition_list(pp, 4,
//                        8, 1, 8, 1,
//                          yo);
//    return 0;

    srand(time(NULL));
    InitializeCriticalSection(&critical);

    struct V global_input;
    char output_filename[64];
    set_variables_global(&global_input, output_filename);

    int thread_num = 1;
//    printf("thread num = ");
//    scanf("%d", &thread_num);

    int P = PARTITIONS_NUM;

    int* shuffled_arr = malloc(P*100*thread_num*sizeof(int));
    for (int i = 0; i < 100*thread_num; ++i) {
        for (int j = 0; j < P; ++j) {
            shuffled_arr[i*P+j] = j;
        }
        shuffle(shuffled_arr+i*P, global_input.pt_num);
    }

    struct Params p;
    int output_int = 0;
    int thread_id = -1;
    memcpy(&p.global_input, &global_input, sizeof(global_input));
    p.output_int = &output_int;
    p.thread_id = &thread_id;
    p.shuffled_arr = shuffled_arr;

    HANDLE threads[100];

    for (int i = 0; i < thread_num; ++i) {
        threads[i] = CreateThread(0, 0, get_all_pt_array_thread, &p, 0, 0);
    }
    for (int i = 0; i < thread_num; ++i) {
        WaitForSingleObject(threads[i], INFINITE);
    }

    //
    /*
     global -> sequences = malloc(12*12*sizeof(int));
    global -> sequences_possible = malloc(12*PARTITIONS_NUM*sizeof(int));
    global -> seq_len = 12;
    global -> seq_num = 12;

    global -> valid_selection_string_len = 72; // the len of each valid string of valid selections
                                            // typically = pts_used, sinch agg has it's own sequence
    global -> valid_selection_string_num = 12;

    global -> valid_selection_strings = malloc(global -> valid_selection_string_num
                                               * global -> valid_selection_string_len
                                               * sizeof(int)); // siz

    */
//    free(glo
    free(shuffled_arr);
    return 0;

}
