/* -*- Mode: c -*-
 * array.h - Generic Array Data Structure
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


#ifndef _ARRAY_H
#define _ARRAY_H

#include <stdio.h>
#include "vector.h"
#include "opt.h"

typedef struct {
        void * CEX_RESTRICT data;
        size_t length;
        size_t alloced; 
        size_t el_size;
        size_t alignment;
} array_t;

/* Generic Operations */
array_t * CEX_make_array(size_t el_size, size_t init_alloc);
void CEX_free_array(array_t *);
void CEX_prealloc_array(array_t *, size_t ensure_alloc);
void CEX_grow_array(array_t *);
void CEX_align_array(array_t *, size_t);
array_t * CEX_copy_array(array_t *);
void CEX_append_ex(array_t *, void *);
void CEX_extend_array(array_t *, array_t *);
array_t * CEX_concat_arrays(array_t *, array_t *);
void CEX_reverse_array(array_t *);
array_t *CEX_reversed_array(array_t *);
void CEX_remove_array_index(array_t *, int);
typedef int (*array_el_comparer)(const void *, const void *);
void CEX_sort_array(array_t *, array_el_comparer);
array_t *CEX_sorted_array(array_t *, array_el_comparer);
void CEX_zero_array_elements(array_t *);

/* helper function for python interface */
void CEX_assign_array_length(array_t *arr, int len);

#define ARR_EL_SIZE(arr) ((arr)->el_size)
#define ARR_LENGTH(arr) ((arr)->length)
#define ARR_ALLOCED(arr) ((arr)->alloced)
#define ARR_DATA(arr) ((arr)->data)
#define ARR_ALIGNMENT(arr) ((arr)->alignment)
#define ARR_DATA_AS(tp, arr) ((tp*)ARR_DATA(arr))
#define ARR_INDEX_AS(tp, arr, inx) ARR_DATA_AS(tp, arr)[inx]
#define ARR_ADDRESS_ELEMENT(arr, inx) ((void *)(ARR_DATA_AS(char, arr) + ((inx) * ARR_EL_SIZE(arr))))

/* #define ARR_FLAGS(arr) ((arr)->flags) */
/* #define ARR_FL_GET(arr, fl) (ARR_FL_GET(arr) & (fl)) */
/* #define ARR_FL_SET(arr, fl) (ARR_FLAGS(arr) |= (fl)) */
/* #define ARR_FL_UNSET(arr, fl) (ARR_FLAGS(arr) &= (~(fl))) */

/* #define ARR_FL_ALIGNED (1<<0) */

/* #define ARR_ALIGNED(arr) ARR_FL_GET(arr, ARR_FL_ALIGNED) */

#define ARR_APPEND(tp, arr, value) do {                         \
        if (unlikely(ARR_LENGTH(arr)==ARR_ALLOCED(arr))) {      \
                CEX_grow_array(arr);                            \
        }                                                       \
        ARR_DATA_AS(tp, arr)[ARR_LENGTH(arr)++] = (value);      \
} while(0);

#define ARR_FOREACH(tp, arr, ptr, counter)                      \
        for (({array_t *_arr = (arr);                           \
               ptr=ARR_DATA_AS(tp, _arr);                       \
               counter=ARR_LENGTH(_arr); });                    \
             counter-- > 0; ptr++)

#define ARR_FOREACH_BACKWARDS(tp, arr, ptr, counter)                    \
        for (({array_t *_arr = (arr);                                   \
                counter=ARR_LENGTH(_arr);                               \
                ptr=ARR_DATA_AS(tp, arr) + counter - 1;});              \
             counter-- > 0; ptr--)

/* use gcc's typeof extension to figureout array type */
#define XARR_FOREACH(arr, ptr, counter)                                 \
        for (({array_t *_arr = (arr);                                   \
               ptr=(__typeof__(ptr))((_arr)->data);                     \
               counter=ARR_LENGTH(_arr);});                             \
             counter-- > 0; ptr++)
/* misserable hack due to limitation of loop scope variables and
 * statement expression */
#define XARR_FOREACH_SLICE(arr, start, end, step, ptr, counter)         \
        for (int _step = step, _nothing GCC_ATTRIBUTE((unused)) =       \
             ({array_t *_arr = (arr);                                   \
               int _start = (start);                                    \
               int _end = (end);                                        \
               ptr=((__typeof__(ptr))((_arr)->data)) + start;           \
                     counter=_end - _start ;});                         \
             counter > 0; ptr+=_step, counter -= _step)

static inline void truncate_array(array_t *, size_t length) GCC_ATTRIBUTE((always_inline));
static inline void clear_array(array_t *) GCC_ATTRIBUTE((always_inline));

/* Character Arrays 
 * Not NULL terminated
 */
array_t *CEX_make_char_array(size_t init_alloc);
array_t *CEX_make_char_array_from_string(const char *str);
array_t *CEX_make_char_array_from_string_ex(const char *str, int length);
void CEX_char_array_as_string(array_t *, char *);
void CEX_print_char_array(FILE *, array_t *);

#define IS_CARR(arr) (ARR_EL_SIZE(arr)==sizeof(char))
#define REQ_CARR(ar) do {                               \
        if (unlikely(!IS_CARR(ar))) {                   \
                Fatal("non char array (el_size=%lu)",    \
                      (unsigned long)ARR_EL_SIZE(ar));   \
        }                                               \
} while (0)                

#define CARR_APPEND(arr, value) ARR_APPEND(char, arr, value)
#define CARR_FOREACH(arr, ptr, counter) ARR_FOREACH(char, arr, ptr, counter)
#define CARR_FOREACH_BACKWARDS(arr, ptr, counter)       \
        ARR_FOREACH_BACKWARDS(char, arr, ptr, counter)


/* Integer Arrays */
array_t * CEX_make_int_array(size_t init_alloc);
array_t * CEX_make_arange(int start, int end, int step);
void CEX_print_int_array(FILE *, array_t *);

#define IS_IARR(arr) (ARR_EL_SIZE(arr)==sizeof(int))
#define REQ_IARR(ar) do {                               \
        if (unlikely(!IS_IARR(ar))) {                   \
                Fatal("non int array (el_size=%lu)",     \
                      (unsigned long)ARR_EL_SIZE(ar));   \
        }                                               \
} while (0)                
#define IARR_APPEND(arr, value) ARR_APPEND(int, arr, value)
#define IARR_FOREACH(arr, ptr, counter) ARR_FOREACH(int, arr, ptr, counter)
#define IARR_FOREACH_BACKWARDS(arr, ptr, counter)       \
        ARR_FOREACH_BACKWARDS(int, arr, ptr, counter)
void CEX_sort_int_array(array_t *);
array_t* CEX_sorted_int_array(array_t *);
/* second array is an Integer array_t */
void CEX_remove_array_of_indices(array_t *, array_t *);
array_t* CEX_splice_array_of_indices(array_t *, array_t *);
int CEX_find_int_array(array_t *arr, int i);
int CEX_in_int_array(array_t *arr, int i);


/* Vector Arrays */
array_t * CEX_make_vec_array(size_t init_alloc);
void CEX_print_vec_array(FILE *fp, array_t *ar);

#define IS_VARR(arr) (ARR_EL_SIZE(arr)==sizeof(vec_t))
#define REQ_VARR(ar) do {                                  \
    if (unlikely(!IS_VARR(ar))) {                          \
            Fatal("non vector array (el_size=%lu)",         \
                  (unsigned long)ARR_EL_SIZE(ar));          \
    }                                                      \
} while (0)                
#define VARR_APPEND(arr, value) ARR_APPEND(vec_t, arr, value)
#define VARR_FOREACH(arr, ptr, counter) ARR_FOREACH(vec_t, arr, ptr, counter)
#define VARR_FOREACH_BACKWARDS(arr, ptr, counter)       \
        ARR_FOREACH_BACKWARDS(vec_t, arr, ptr, counter)

#include "array-inline.h"

#endif /* _ARRAY_H */

