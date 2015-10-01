#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "max_weighted_match.h"
#define false 0
#define true 1
#define STACK_VERSION
//#undef STACK_VERSION
static inline val_type max(val_type a,val_type b){
    return a>b?a:b;
}
#define cell(mat, shape_2, c1, c2) ((mat)[(shape_2) * (c1) + (c2)])
//assert len(ret) == shape_1
static inline void max_every_row(const val_type* mat, uint32_t m, uint32_t n, val_type* ret){
    for (uint32_t i = 0; i < m; i++) {
        for (uint32_t j = 0; j < n;j++) {
            ret[i] = max(ret[i], cell(mat, n,i,j));
        }
    }
}

static inline void transpose(val_type* dst,const val_type* src, uint32_t m, uint32_t n){
    for (uint32_t i = 0; i< m; i++) {
        for (uint32_t j = 0; j < n; j++) {
            uint32_t dst_pos = j * m  + i;
            dst[dst_pos] = cell(src,n,i,j);
        }
    }
}

val_type kuhn_munkres(const val_type* mat, uint32_t m, uint32_t n, uint32_t *ret_mapping){
    assert(m>0 && n>0 && mat);
    int transposed = false;
    val_type* dst = NULL;
#ifdef STACK_VERSION
    val_type dst_stack[m*n];
#endif
    if (m>n) {
        transposed = true;
        //一定要m<=n
        //如果不对,那么转秩
#ifndef STACK_VERSION
        dst = malloc(sizeof(val_type)*m*n);
#else
        dst = dst_stack;
#endif
        transpose(dst, mat,m,n);
        uint32_t tmp = m;
        m = n;
        n = tmp;
        mat = dst;
    }
#ifndef STACK_VERSION
    val_type* us = malloc(sizeof(val_type)*m); //label for every row
    val_type* vs = malloc(sizeof(val_type)*n);//label for every column

    int32_t* match_u = malloc(sizeof(int32_t)*m);//match for every row
    int32_t* match_v = malloc(sizeof(int32_t)*n);//match for every column
    int32_t* tmp_match = malloc(sizeof(int32_t)*n);

    int32_t* stack = malloc(sizeof(int32_t)*(m+n));
#else
    val_type us[m]; //label for every row
    val_type vs[n]; //label for every column

    int32_t match_u[m];//match for every row
    int32_t match_v[n];//match for every column
    int32_t tmp_match[n];

    int32_t stack[m+n];
#endif
    int32_t p,q;
    //dual就是最小化点的权重总和
    //因为一条边两个端点的总和肯定大于等于边的权重
    memset(us, 0, sizeof(val_type)*m);
    memset(vs, 0, sizeof(val_type)*n);
    max_every_row(mat,m,n,us);
    //初始化结束

    memset(match_u, -1, sizeof(uint32_t)*m);
    memset(match_v, -1, sizeof(uint32_t)*n);
    for (uint32_t i = 0; i < m; i++) {//每次增广路增加一条边,最多增加m次(m<=n)

        memset(tmp_match, -1, sizeof(uint32_t)*n);//初始时, 没有任何column有match过

        //匈牙利算法是对于每个点尝试加入最大匹配,dfs进行踢人.
        //因此对于每个没有被匹配上的点,尝试加入匹配,并踢开原来匹配的人.
        //match_u[i]保证当前访问的点是没有被匹配上的.
        //stack保存了当前被踢开的人的栈,首先i加入
        //踢开一个点,然后这个点加入,踢开另一个点

        for (stack[p=q=0] = i; p <= q && match_u[i] < 0; p++) {

            uint32_t k = stack[p];

            //对于集合X中的点k,遍历集合Y中的点j.
            //寻找点:如果两个点的标号和正好等于边的权重,而且这个点没有被访问过.
            //如果找到了:
            //    将这个点j原来匹配的点压入栈中,将当前点j的临时匹配设为k.
            //    如果这个点原来没有匹配. 那么这次寻找增广路已经成功了.
            //        那么
            
            for (uint32_t j = 0; j < n && match_u[i] < 0; j++) {
                if (us[k] + vs[j] == cell(mat,n,k,j) && tmp_match[j] < 0) {
                    stack[++q] = match_v[j];
                    tmp_match[j] = k;
                    if (stack[q] < 0) {
                        for (p = j; p>=0;j = p) {//将临时匹配设为正式匹配.
                            match_v[j] = k = tmp_match[j];
                            p = match_u[k];
                            match_u[k] = j;
                        }
                    }
                }
            }
        }

        if (match_u[i] < 0) {//如果当前没有完美匹配.
            //我们求当前相等子图的完备匹配失败了，是因为对于某个X顶点，我们找不到一条从它出发的交错路。
            //这时我们获得了一棵交错树，它的叶子结点全部是X顶点。
            i--;
            p=max_val;
            for (uint32_t k=0;k<=q;k++) {
                for (uint32_t j = 0; j< n;j++) {
                    if (tmp_match[j] < 0 && us[stack[k]] + vs[j] - cell(mat,n,stack[k],j) < p) {
                        p = us[stack[k]] + vs[j] - cell(mat,n,stack[k],j);
                    }
                }
            }
            
		//但是朴素的实现方法， 时间复杂度 为O(n^4)——需要找O(n)次增广路，
		//每次增广最多需要修改O(n)次顶标，每次修改顶标时由于要枚举边来求d值，复杂度为O(n^2)。
		//实际上KM算法的复杂度是可以做到O(n^3)的。
		//我们给每个Y顶点一个“松弛量” 函数 slack，每次开始找增广路时初始化为无穷大。在寻找增广路的过程中，检查边(i,j)时，如果它不在相等子图中，则让slack[j]变成原值与A[ i ]+B[j]-w[i,j]的较小值。
		
        }
        //现在我们把交错树中X顶点的顶标全都减小某个值d，Y顶点的顶标全都增加同一个值d
        for (uint32_t j = 0; j < n;j++){
            vs[j] += tmp_match[j] < 0?0:p;
        }
        for (uint32_t k = 0; k<= q;){
            //printf("stack[%d]=%d\n", k, stack[k]);
            if (stack[k] < 0) {
                break;
            }
            us[stack[k++]] -= p; 
        }
        //因为每次都保证了l1[i] + l2[j] >= w[i,j]的
        //1）两端都在交错树中的边(i,j)，l1[ i ]+l2[j]的值没有变化。也就是说，它原来属于相等子图，现在仍属于相等子图。
        //2）两端都不在交错树中的边(i,j)，l1[ i ]和l2[j]都没有变化。也就是说，它原来属于（或不属于）相等子图，现在仍属于（或不属于）相等子图。
        //3）X端不在交错树中，Y端在交错树中的边(i,j)，它的l1[ i ]+l2[j]的值有所增大。它原来不属于相等子图，现在仍不属于相等子图。
        //4）X端在交错树中，Y端不在交错树中的边(i,j)，它的l1[ i ]+l2[j]的值有所减小。也就说，它原来不属于相等子图，现在可能进入了相等子图，因而使相等子图得到了扩大。

    }
    uint32_t ret = 0;
    for (uint32_t i = 0; i<m;i++) {
        ret += cell(mat,n,i,match_u[i]);
    }
    if (ret_mapping) {
        if (transposed) {
            for (uint32_t i = 0; i<m;i++) {
                ret_mapping[match_u[i]] = i;
            }
        } else {
            for (uint32_t i = 0;i<m;i++) {
                ret_mapping[i] = match_u[i]; 
            }
        }
    }
#ifndef STACK_VERSION
    if (transposed) {
        free(dst);
    }
    free(us);
    free(vs);
    free(match_u);
    free(match_v);
    free(tmp_match);
    free(stack);
#endif
    return ret;
}
