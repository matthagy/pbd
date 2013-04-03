/* -*- Mode: c -*-
 * msg-inline.h - Composite Messaging Inlines
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

#include <assert.h>

/* Writing Messages */
void CEX_msg_grow_writer(msg_t *msg);

static inline void
CEX_msg_write_char(msg_t *msg, char c)
{
        REQ_WMSG(msg);
        if (unlikely(MSG_EOFP(msg))) {
                CEX_msg_grow_writer(msg);
        }
        *(MSG_PTR(msg)++) = c;
}

static inline void 
CEX_msg_write_uint(msg_t *msg, unsigned int ui)
{
        if ((ui >> 24) & ~0xffU) {
                Fatal("unsigned int %ui too large", ui);
        }
        CEX_msg_write_char(msg, (char)((ui>>24) & 0xffU));
        CEX_msg_write_char(msg, (char)((ui>>16) & 0xffU));
        CEX_msg_write_char(msg, (char)((ui>>8) & 0xffU));
        CEX_msg_write_char(msg, (char)((ui>>0) & 0xffU));
}

static inline void
CEX_msg_write_int(msg_t *msg, int i)
{
        int sign = i >= 0 ? 1 : 0;
        unsigned int base = (unsigned int)(sign ? i : -i);
        CEX_msg_write_uint(msg, base | (sign << 31));
}

static inline void
CEX_msg_write_array(msg_t *msg, msg_element_writer writer, array_t *arr)
{
        unsigned int length = (unsigned int)ARR_LENGTH(arr);
        if (unlikely((size_t)length != ARR_LENGTH(arr))) {
                Fatal("cannot coere array length to unsigned int");
        }
        CEX_msg_write_uint(msg, length);
        for (int i=0; i < (int)length; i++) {
                writer(msg, arr, i);
        }
}

/* Reading Message */
static inline char 
CEX_msg_read_char(msg_t *msg)
{
        REQ_RMSG(msg);
        if (unlikely(MSG_EOFP(msg))) {
                Fatal("underflow in message reading");
        }
        return *(MSG_PTR(msg)++);
}

static inline unsigned int
CEX_msg_read_uint(msg_t *msg)
{
        unsigned int a,b,c,d;
        a = (unsigned int)CEX_msg_read_char(msg) & 0xffUL;
        b = (unsigned int)CEX_msg_read_char(msg) & 0xffUL;
        c = (unsigned int)CEX_msg_read_char(msg) & 0xffUL;
        d = (unsigned int)CEX_msg_read_char(msg) & 0xffUL;
        return (a<<24) | (b<<16) | (c<<8) | d;
}

static inline int
CEX_msg_read_int(msg_t *msg)
{
        unsigned int base = CEX_msg_read_uint(msg);
        int sign = (int)(base >> 31);
        int ibase = (int)(base & ((1UL<<31UL)-1UL));
        return likely(sign==1) ? ibase : -ibase;
}

static inline array_t *
CEX_msg_read_array(msg_t *msg, size_t el_size, msg_element_reader reader)
{
        size_t length = (size_t)CEX_msg_read_uint(msg);
        array_t *arr = CEX_make_array(el_size, length);
        for (int i=length; i-- > 0;) {
                reader(msg, arr);
        }
        if (unlikely(ARR_LENGTH(arr)!=length)) {
                Fatal("expected array of length %lu; reader gave %lu elements",
                      (unsigned long)length, (unsigned long)ARR_LENGTH(arr));
        }
        return arr;
}
