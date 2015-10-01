#ifndef __MAX_WEIGHTED_MATCH_H__
#define __MAX_WEIGHTED_MATCH_H__ 1
#include <stdint.h>
#include <limits.h>

#define val_type int32_t
#define max_val INT_MAX
/**
 * 权重值都是非负数
 * 暂时不引入float的原因就是整数操作快一点, 
 * 如果实际应用中有浮点数,也只需要乘以一个大整数,转换成整数.
 * 
 * @author chenkovsky (9/30/2015)
 * 
 * @param mat 矩阵
 * @param shape_1 矩阵dim1
 * @param shape_2 矩阵dim2
 * @param ret_mapping 
 *                      如果传入数组,那么将返回的映射关系写在这个数组里面,
 *                      ret_mapping[i]=j说明dim1为i映射到dim2的j
 *  
 * @return val_type 最大权重和
 */
val_type kuhn_munkres(const val_type* mat, uint32_t shape_1, uint32_t shape_2, uint32_t *ret_mapping);
#endif

