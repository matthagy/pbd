/* -*- Mode: c -*-
 * array.c - Array Data Structures
 *--------------------------------------------------------------------------
 * Copyright (C) 2009, Matthew Hagy (hagy@gatech.edu)
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <math.h>
#include <assert.h>

#include "array.h"
#include "mem.h"
#include "debug.h"

array_t * 
CEX_make_array(size_t el_size, size_t init_alloc)
{
        array_t *arr;

        if (unlikely(el_size<=0)) {
                Fatal("bad el_size %lu", (unsigned long)el_size);
        }
        arr = XNEW(array_t, 1);
        ARR_LENGTH(arr) = 0;
        ARR_ALLOCED(arr) = init_alloc;
        ARR_EL_SIZE(arr) = el_size;
        ARR_DATA(arr) = CEX_malloc(el_size * init_alloc);
        ARR_ALIGNMENT(arr) = 0;
        return arr;
}

static inline void
free_array_data(array_t *arr)
{
        if (ARR_ALIGNMENT(arr)==0) {
                CEX_free(ARR_DATA(arr));
        } else {
                CEX_aligned_free(ARR_DATA(arr));
        }
}

void 
CEX_free_array(array_t *arr)
{
        if (arr) {
                free_array_data(arr);
                CEX_free(arr);
        }
}

void 
CEX_prealloc_array(array_t *arr, size_t sz)
{
        if (ARR_ALLOCED(arr) >= sz || sz==0)
                return;
        ARR_DATA(arr) = ARR_ALIGNMENT(arr)==0 ? 
                          CEX_realloc(ARR_DATA(arr), ARR_EL_SIZE(arr) * sz) :
                          CEX_aligned_realloc(ARR_DATA(arr), ARR_EL_SIZE(arr) * sz,
                                    ARR_ALIGNMENT(arr));
        ARR_ALLOCED(arr) = sz;
}

void 
CEX_grow_array(array_t *arr)
{
        if (ARR_ALLOCED(arr)==0) {
                CEX_prealloc_array(arr, 16);
        } else {
                CEX_prealloc_array(arr, ARR_ALLOCED(arr) << 1);
        }
}

void 
CEX_align_array(array_t *arr, size_t alignment)
{
        if (alignment<=ARR_ALIGNMENT(arr)) 
                return;
        void *new_data = CEX_aligned_malloc(ARR_ALLOCED(arr)*ARR_EL_SIZE(arr),
                                            alignment);
        memcpy(new_data, ARR_DATA(arr), ARR_LENGTH(arr)*ARR_EL_SIZE(arr));
        free_array_data(arr);        
        ARR_DATA(arr) = new_data;
        ARR_ALIGNMENT(arr) = alignment;
}


array_t *
CEX_copy_array(array_t *arr)
{
        array_t *cp;

        cp = CEX_make_array(ARR_EL_SIZE(arr), ARR_LENGTH(arr));
        ARR_LENGTH(cp) = ARR_LENGTH(arr);
        memcpy(ARR_DATA(cp), ARR_DATA(arr), ARR_EL_SIZE(arr)*ARR_LENGTH(arr));
        return cp;
}

void 
CEX_append_ex(array_t *arr, void *data)
{
        CEX_prealloc_array(arr, ARR_LENGTH(arr)+1);
        memcpy(ARR_ADDRESS_ELEMENT(arr, ARR_LENGTH(arr)),
               data, ARR_EL_SIZE(arr));
        ARR_LENGTH(arr) ++;
}

void 
CEX_extend_array(array_t *dst, array_t *src)
{
        if (unlikely(ARR_EL_SIZE(dst) != ARR_EL_SIZE(src))) {
                Fatal("extending from incompatible element sizes; dst=%lu src=%lu",
                      (unsigned long)ARR_EL_SIZE(dst), 
                      (unsigned long)ARR_EL_SIZE(src));
        }
        CEX_prealloc_array(dst, ARR_LENGTH(dst) + ARR_LENGTH(src));
        memcpy(ARR_ADDRESS_ELEMENT(dst, ARR_LENGTH(dst)),
               ARR_ADDRESS_ELEMENT(src, 0),
               ARR_LENGTH(src) * ARR_EL_SIZE(src));
        ARR_LENGTH(dst) += ARR_LENGTH(src);
}

array_t *
CEX_concat_arrays(array_t *a, array_t *b)
{
        array_t *cp;

        cp = CEX_copy_array(a);
        CEX_extend_array(cp, b);
        return cp;
}


/* this could be optimized for langage size elements (ie. chars and integers)
 * currently not used anywhere where this would really help
 */
void
CEX_reverse_array(array_t *arr)
{
        int size;
        char *start, *end, *tmp;

        size = ARR_EL_SIZE(arr);
        start = ARR_DATA_AS(char, arr);
        end = (char *)ARR_ADDRESS_ELEMENT(arr, ARR_LENGTH(arr)-1);
        tmp = (char *)alloca(size);
        while (start < end) {
                memcpy(tmp, start, size);
                memcpy(start, end, size);
                memcpy(end, tmp, size);
                start += size;
                end -= size;
        }
}

array_t *
CEX_reversed_array(array_t *arr)
{
        array_t *cp;

        cp = CEX_copy_array(arr);
        CEX_reverse_array(cp);
        return cp;
}

static inline int normalize_index(array_t *arr, int i)
        GCC_ATTRIBUTE((always_inline));

static inline int
normalize_index(array_t *arr, int i)
{
        if (i<0) {
                i += ARR_LENGTH(arr);
        }
        if (unlikely(i<0 || i>=ARR_LENGTH(arr))) {
                Fatal("bad array index %d for array of length %lu",
                      i, (unsigned long)ARR_LENGTH(arr));
        }
        return i;
}

void 
CEX_remove_array_index(array_t *arr, int i)
{
        i = normalize_index(arr, i);
        memcpy(ARR_ADDRESS_ELEMENT(arr, i),
               ARR_ADDRESS_ELEMENT(arr, i+1),
               ARR_EL_SIZE(arr) * (ARR_LENGTH(arr) - i - 1));
        ARR_LENGTH(arr) --;
}

void 
CEX_sort_array(array_t *arr, array_el_comparer cmp)
{
        qsort(ARR_DATA(arr), ARR_LENGTH(arr), ARR_EL_SIZE(arr), cmp);
}

array_t *
CEX_sorted_array(array_t *arr, array_el_comparer cmp)
{
        array_t *cp;
        
        cp = CEX_copy_array(arr);
        CEX_sort_array(cp, cmp);
        return cp;
}

void 
CEX_zero_array_elements(array_t *arr)
{
        bzero(ARR_DATA(arr), ARR_LENGTH(arr) * ARR_EL_SIZE(arr));
}

/* Character Arrays */
array_t *
CEX_make_char_array(size_t init_alloc)
{
        return CEX_make_array(sizeof(char), init_alloc);
}

array_t *
CEX_make_char_array_from_string(const char *str)
{
        return CEX_make_char_array_from_string_ex(str, strlen(str));
}

array_t *
CEX_make_char_array_from_string_ex(const char *str, int length)
{
        array_t *arr = CEX_make_char_array(length);
        memcpy(ARR_DATA(arr), (void *)str, length);
        ARR_LENGTH(arr) = length;
        return arr;
}

void
CEX_char_array_as_string(array_t *arr, char *buffer)
{
        int length = ARR_LENGTH(arr);
        memcpy((void *)buffer, ARR_DATA(arr), length);
        buffer[length] = '\0';
}

void
CEX_print_char_array(FILE *fp, array_t *arr)
{
        REQ_CARR(arr);
        int first = 1;
        fprintf(fp, "[");
        char *ptr; int counter;
        CARR_FOREACH(arr, ptr, counter) {
                if (first) {
                        first = 0;
                } else {
                        fprintf(fp, " ");
                }
                fprintf(fp, "%d", (int)*ptr);
        }
        fprintf(fp, "]\n");
}


/* Integer Arrays */
array_t * 
CEX_make_int_array(size_t init_alloc)
{
        return CEX_make_array(sizeof(int), init_alloc);
}

array_t *
CEX_make_arange(int start, int end, int step)
{
        if (unlikely(step == 0 && start!=end)  ||
            unlikely((end-start) * step < 0)) {
                Fatal("%d step for range [%d:%d) is undefined",
                      step, start, end);
        }
        int i;
        array_t *arr;

        arr = CEX_make_int_array((int)ceil(fabs((double)(end - start) / (double)step)));
        if (step>0) {
                for (i=start; i<end; i+=step) {
                        IARR_APPEND(arr, i);
                }
        } else {
                for (i=start; i>end; i+=step) {
                        IARR_APPEND(arr, i);
                }
        }
        return arr;

}

void 
CEX_print_int_array(FILE *fp, array_t *arr)
{
        int first;
        int *ptr, counter;

        REQ_IARR(arr);
        first = 1;
        fprintf(fp, "[");
        IARR_FOREACH(arr, ptr, counter) {
                if (first) {
                        first = 0;
                } else {
                        fprintf(fp, " ");
                }
                fprintf(fp, "%d", *ptr);
        }
        fprintf(fp, "]\n");
}

static int iarr_cmp(int *a, int *b)
{
        return *a - *b;
}

void 
CEX_sort_int_array(array_t *arr)
{
        REQ_IARR(arr);
        CEX_sort_array(arr, (array_el_comparer)&iarr_cmp);
}

array_t *
CEX_sorted_int_array(array_t *arr)
{
        REQ_IARR(arr);
        return CEX_sorted_array(arr, (array_el_comparer)&iarr_cmp);
}

void 
CEX_remove_array_of_indices(array_t *data, array_t *indices)
{
        array_t *sorted_indices;
        int *iptr, icnt;

        REQ_IARR(indices);
        sorted_indices = CEX_copy_array(indices);
        IARR_FOREACH(sorted_indices, iptr, icnt) {
                *iptr = normalize_index(data, *iptr);
        }
        CEX_sort_int_array(sorted_indices);
        IARR_FOREACH_BACKWARDS(sorted_indices, iptr, icnt) {
                CEX_remove_array_index(data, *iptr);
        }
        CEX_free_array(sorted_indices);
}

array_t * 
CEX_splice_array_of_indices(array_t *data, array_t *indices)
{
        array_t *cp;
        int *iptr, icnt;

        REQ_IARR(indices);
        cp = CEX_make_array(ARR_EL_SIZE(data), ARR_LENGTH(indices));
        IARR_FOREACH(indices, iptr, icnt) {
                CEX_append_ex(cp, ARR_ADDRESS_ELEMENT(data, *iptr));
        }
        return cp;
        
}

int
CEX_find_int_array(array_t *arr, int key)
{
       int *iptr, counter;

        IARR_FOREACH(arr, iptr, counter) {
                if (unlikely(*iptr == key)) {
                        return (int)(iptr - ARR_DATA_AS(int, arr));
                }
        }
        return -1;
}

int 
CEX_in_int_array(array_t *arr, int key)
{
        return CEX_find_int_array(arr, key) != -1;
}
 
/* Vector Arrays */
array_t * 
CEX_make_vec_array(size_t init_alloc)
{
        return CEX_make_array(sizeof(vec_t), init_alloc);
}

void 
CEX_print_vec_array(FILE *fp, array_t *arr)
{
        int first;
        int counter;
        vec_t *ptr;
        
        REQ_VARR(arr);
        first = 1;
        fprintf(fp, "[");
        VARR_FOREACH(arr, ptr, counter) {
                if (first) {
                        first = 0;
                } else {
                        fprintf(fp, " ");
                }
                fprintf(fp, Vec3G_FRMT, Vec3_ARGS(*ptr));
        }
        fprintf(fp, "]\n");
}


void 
CEX_assign_array_length(array_t *arr, int len)
{
        ARR_LENGTH(arr) = len;
}
