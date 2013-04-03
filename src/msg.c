/* -*- Mode: c -*-
 * msg.c  - Composite Messages
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

#include "msg.h"
#include "mem.h"

static inline msg_t *
make_msg(int mode, array_t *buffer)
{
        if (!(mode==MSG_R || mode==MSG_W)) {
                Fatal("bad message mode %d", mode);
        }
        msg_t *msg = XNEW(msg_t, 1);
        MSG_MODE(msg) = mode;
        MSG_BUFFER(msg) = buffer;
        MSG_START(msg) = MSG_PTR(msg) = ARR_DATA_AS(char, buffer);
        MSG_END(msg) = MSG_PTR(msg) + ARR_ALLOCED(buffer);
        MSG_OWN_BUFFER(msg) = 1;
        return msg;
}


void
CEX_free_msg(msg_t *msg)
{
        if (msg) {
                if (!MSG_OWN_BUFFER(msg)) {
                        Fatal("freeing message view");
                }
                CEX_free_array(MSG_BUFFER(msg));
                CEX_free(msg);
        }
}
void
CEX_msg_grow_writer(msg_t *msg)
{
        REQ_MSG_OWN_BUFFER(msg);
        REQ_WMSG(msg);
        int length = CEX_msg_len(msg);
        CEX_prealloc_msg(msg, length==0 ? 64 : length<<1);
        assert(!MSG_EOFP(msg));
}

void
CEX_prealloc_msg(msg_t *msg, int req_len)
{
        REQ_MSG_OWN_BUFFER(msg);
        array_t *buffer = MSG_BUFFER(msg);
        int start_offset = (int)(MSG_START(msg) - ARR_DATA_AS(char, buffer));
        int ptr_offset = CEX_msg_tell(msg);
        CEX_prealloc_array(buffer, req_len+start_offset);
        MSG_START(msg) = ARR_DATA_AS(char, buffer) + start_offset;
        MSG_END(msg) = ARR_DATA_AS(char, buffer) + ARR_ALLOCED(buffer);
        MSG_PTR(msg) = MSG_START(msg) + ptr_offset;
        assert(CEX_msg_tell(msg)==ptr_offset);
}

int
CEX_msg_len(msg_t *msg)
{
        return (int)(MSG_END(msg) - MSG_START(msg));        
}

int
CEX_msg_tell(msg_t *msg)
{
        return (int)(MSG_PTR(msg) - MSG_START(msg));
}

void
CEX_msg_seek(msg_t *msg, int index)
{
        if (index<0 || index >= CEX_msg_len(msg)) {
                Fatal("seek %d out of range of message of length %d",
                      index, CEX_msg_len(msg));
        }
        MSG_PTR(msg) = MSG_START(msg) + index;
}

/* Writing Messages */
msg_t *
CEX_make_write_msg(int init_alloc)
{
        return make_msg(MSG_W, CEX_make_array(sizeof(char), init_alloc));
}

void 
CEX_finalize_write_msg(msg_t *msg)
{
        REQ_WMSG(msg);
        REQ_MSG_OWN_BUFFER(msg);
        MSG_END(msg) = MSG_PTR(msg);
}

static void
char_writer(msg_t *msg, array_t *arr, int i)
{
        CEX_msg_write_char(msg, ARR_INDEX_AS(char, arr, i));
}

void
CEX_msg_write_char_array(msg_t *msg, array_t *arr)
{
        REQ_WMSG(msg);
        REQ_CARR(arr);
        CEX_msg_write_array(msg, char_writer, arr);
}

static void
int_writer(msg_t *msg, array_t *arr, int i)
{
        CEX_msg_write_int(msg, ARR_INDEX_AS(int, arr, i));
}

void
CEX_msg_write_int_array(msg_t *msg, array_t *arr)
{
        REQ_WMSG(msg);
        REQ_IARR(arr);
        CEX_msg_write_array(msg, int_writer, arr);
}

static inline void
trim_zeros(char *ptr)
{
        int len = strlen(ptr) + 1;
        while (likely(*ptr=='0' && ptr[-1] >= '0' && ptr[-1] <= '9')) {
                memmove(ptr, ptr+1, len);
                ptr--;
        }
}

void 
CEX_msg_write_double(msg_t *msg, double value)
{
        REQ_WMSG(msg);
        char buffer[256];
        int written = snprintf(buffer, sizeof(buffer), "%.10le", value);
        if(written > sizeof(buffer)) {
                Fatal("overflow writing %le to buffer length %lu", value, sizeof(buffer));
        }
#if 1
        //trim_zeros(buffer + strlen(buffer) - 1);
        char *ptr = buffer;
        while (*(ptr++) != 'e');
        ptr--; ptr--;
        trim_zeros(ptr);
#endif
        //printf("%le %.256s\n", value, buffer);
        array_t *arr = CEX_make_char_array_from_string(buffer);
        CEX_msg_write_char_array(msg, arr);
        CEX_free_array(arr);
}

void
CEX_msg_write_vec(msg_t *msg, vec_t vec)
{
        REQ_WMSG(msg);
        CEX_msg_write_double(msg, vec.x);
        CEX_msg_write_double(msg, vec.y);
        CEX_msg_write_double(msg, vec.z);
}

static void
vec_writer(msg_t *msg, array_t *arr, int i)
{
        REQ_WMSG(msg);
        REQ_VARR(arr);
        CEX_msg_write_vec(msg, ARR_INDEX_AS(vec_t, arr, i));
}

void
CEX_msg_write_vec_array(msg_t *msg, array_t *arr)
{
        REQ_WMSG(msg);
        REQ_VARR(arr);
        CEX_msg_write_array(msg, &vec_writer, arr);
}

/* Reading Message */
msg_t *
CEX_make_read_msg(array_t *buffer)
{
        return make_msg(MSG_R, buffer);
}

static void
char_reader(msg_t *msg, array_t *arr)
{
        REQ_CARR(arr);
        ARR_APPEND(char, arr, CEX_msg_read_char(msg));
}

array_t *
CEX_msg_read_char_array(msg_t *msg)
{
        return CEX_msg_read_array(msg, sizeof(char), &char_reader);
}

static void 
uint_reader(msg_t *msg, array_t *arr)
{
        ARR_APPEND(unsigned int, arr, CEX_msg_read_uint(msg));
}

array_t *
CEX_msg_read_uint_array(msg_t *msg)
{
        return CEX_msg_read_array(msg, sizeof(unsigned int), &uint_reader);
}

static void 
int_reader(msg_t *msg, array_t *arr)
{
        REQ_IARR(arr);
        ARR_APPEND(int, arr, CEX_msg_read_int(msg));
}

array_t *
CEX_msg_read_int_array(msg_t *msg)
{
        return CEX_msg_read_array(msg, sizeof(int), &int_reader);
}

double
CEX_msg_read_double(msg_t *msg)
{
        array_t *arr = CEX_msg_read_char_array(msg);
        char buffer[ARR_LENGTH(arr)+1];
        CEX_char_array_as_string(arr, buffer);
        double value;
        sscanf(buffer, "%le", &value);
        CEX_free_array(arr);
        return value;
}

static void
double_reader(msg_t *msg, array_t *arr)
{
        ARR_APPEND(double, arr, CEX_msg_read_double(msg));
}

array_t *
CEX_msg_read_double_array(msg_t *msg)
{
        return CEX_msg_read_array(msg, sizeof(double), &double_reader);
}

vec_t
CEX_msg_read_vec(msg_t *msg)
{
        vec_t vec;
        vec.x = CEX_msg_read_double(msg);
        vec.y = CEX_msg_read_double(msg);
        vec.z = CEX_msg_read_double(msg);
        return vec;
}

static void
vec_reader(msg_t *msg, array_t *arr)
{
        REQ_VARR(arr);
        ARR_APPEND(vec_t, arr, CEX_msg_read_vec(msg));
}

array_t *
CEX_msg_read_vec_array(msg_t *msg)
{
        return CEX_msg_read_array(msg, sizeof(vec_t), &vec_reader);
}
                
