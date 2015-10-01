#include <assert.h>
#include <stdio.h>
#include "max_weighted_match.h"


void test1(){
    uint32_t shape_1 = 4;
    uint32_t shape_2 = 4;
    const val_type mat[] = {80, 40, 50, 46,
                            40, 70, 20, 25,
                            30, 10, 20, 30,
                            35, 20, 25, 30};
    uint32_t shape_mapping[4];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    uint32_t  map_size = 4;
    uint32_t correct_ret[] = {
        0,1,3,2
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
    assert(ret == 205); 
}

void test2(){
    uint32_t shape_1 = 5;
    uint32_t shape_2 = 5;
    const val_type mat[] = {12, 8, 11, 18, 11,
                            14, 22, 8, 12, 14,
                            14, 14, 16, 14, 15,
                            19, 11, 14, 17, 15,
                            13, 9, 17, 20, 11};
    uint32_t shape_mapping[5];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    uint32_t  map_size = 5;
    uint32_t correct_ret[] = {
        3,1,4,0,2
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
    assert(ret == 91);
}


void test3(){
    uint32_t shape_1 = 5;
    uint32_t shape_2 = 5;
    const val_type mat[] = {20, 14, 6, 10, 22,
                            16, 8, 22, 20, 10,
                            8, 6, 24, 14, 12,
                            20, 22, 2, 8, 6,
                            4, 16, 22, 6, 24};
    uint32_t shape_mapping[5];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 110);
    uint32_t  map_size = 5;
    uint32_t correct_ret[] = {
        0,3,2,1,4
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}


void test4(){
    uint32_t shape_1 = 3;
    uint32_t shape_2 = 3;
    const val_type mat[] = {1, 4, 5,
                            5, 7, 6,
                            5, 8, 8};
    uint32_t shape_mapping[3];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 18);
    uint32_t  map_size = 3;
    uint32_t correct_ret[] = {
        2,0,1
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}

void test5(){
    uint32_t shape_1 = 3;
    uint32_t shape_2 = 3;
    const val_type mat[] = {5, 9, 1,
                            10, 3, 2,
                            8, 7, 4};
    uint32_t shape_mapping[3];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 23);
    uint32_t  map_size = 3;
    uint32_t correct_ret[] = {
        1,0,2
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}

void test6(){
    uint32_t shape_1 = 4;
    uint32_t shape_2 = 4;
    const val_type mat[] = {30, 50, 63, 45,
                            42, 36, 54, 48,
                            31, 41, 67, 21,
                            49, 54, 37, 28};
    uint32_t shape_mapping[4];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 214);
    uint32_t  map_size = 4;
    uint32_t correct_ret[] = {
        1,3,2,0
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}

void test7(){
    uint32_t shape_1 = 5;
    uint32_t shape_2 = 5;
    const val_type mat[] = {3, 5, 6, 7, 1,
                            2, 2, 0, 2, 4,
                            2, 4, 4, 1, 0,
                            0, 1, 1, 0, 0,
                            1, 2, 1, 6, 3};
    uint32_t shape_mapping[5];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 20);
    uint32_t  map_size = 5;
    uint32_t correct_ret[] = {
        2,4,1,0,3
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}

void test8(){
    uint32_t shape_1 = 5;
    uint32_t shape_2 = 5;
    const val_type mat[] = {5, 2, 6, 3, 3,
                            4, 7, 6, 5, 4,
                            2, 6, 3, 1, 6,
                            4, 5, 2, 7, 4,
                            3, 2, 5, 4, 3};
    uint32_t shape_mapping[5];
    uint32_t ret = kuhn_munkres(mat, shape_1, shape_2, shape_mapping);
    assert(ret == 30);
    uint32_t  map_size = 5;
    uint32_t correct_ret[] = {
        0,1,4,3,2
    };
    for (int i = 0; i < map_size;i++) {
        //printf("shape_mapping[%d] = %d, correct_ret[%d] = %d\n", i, shape_mapping[i], i ,correct_ret[i]);
        assert(shape_mapping[i] == correct_ret[i]);
    }
}


int main(int argc, char *argv[]) {
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
}

