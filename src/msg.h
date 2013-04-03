/* -*- Mode: c -*-
 * msg.h - Composite Messages
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


#ifndef _MSG_H
#define _MSG_H

#include "opt.h"
#include "debug.h"
#include "array.h"


#define MSG_R ((int)'r')
#define MSG_W ((int)'w')

typedef struct {
        array_t *buffer;
        char *start, *end, *ptr;
        int mode, own_buffer;
} msg_t;

#define MSG_BUFFER(msg) ((msg)->buffer)
#define MSG_START(msg) ((msg)->start)
#define MSG_END(msg) ((msg)->end)
#define MSG_PTR(msg) ((msg)->ptr)
#define MSG_MODE(msg) ((msg)->mode)
#define MSG_OWN_BUFFER(msg) ((msg)->own_buffer)
#define MSG_EOFP(msg) (MSG_PTR(msg)==MSG_END(msg))

void CEX_free_msg(msg_t *);
int CEX_msg_len(msg_t *); //GCC_ATTRIBUTE((pure));
int CEX_msg_tell(msg_t *); //GCC_ATTRIBUTE((pure));
void CEX_msg_seek(msg_t *, int);
void CEX_prealloc_msg(msg_t *msg, int req_len);

#define REQ_MSG_OWN_BUFFER(msg) do {                     \
      if (unlikely(!MSG_OWN_BUFFER(msg))) {              \
           Fatal("cannot accept message view");          \
      }                                                  \
} while (0)                

/* Writing Messages */
msg_t *CEX_make_write_msg(int init_alloc);
#define IS_WMSG(msg) (MSG_MODE(msg)==MSG_W)
#define REQ_WMSG(msg) do {                                \
      if (unlikely(!IS_WMSG(msg))) {                      \
              Fatal("non writing message (mode=%d)",     \
                    MSG_MODE(msg));                      \
      }                                                  \
} while (0)                

void CEX_finalize_write_msg(msg_t *);
static inline void CEX_msg_write_char(msg_t *, char) GCC_ATTRIBUTE((always_inline));
/* 32-bit integers */
static inline void CEX_msg_write_uint(msg_t *, unsigned int) GCC_ATTRIBUTE((always_inline));
static inline void CEX_msg_write_int(msg_t *, int) GCC_ATTRIBUTE((always_inline));
typedef void (*msg_element_writer)(msg_t *, array_t *, int i);
static inline void CEX_msg_write_array(msg_t *, msg_element_writer, array_t *arr) GCC_ATTRIBUTE((always_inline));

void CEX_msg_write_char_array(msg_t *, array_t *);
void CEX_msg_write_double(msg_t *, double);
void CEX_msg_write_vec(msg_t *, vec_t);
void CEX_msg_write_int_array(msg_t *, array_t *);
void CEX_msg_write_vec_array(msg_t *, array_t *);
typedef void (*submsg_writer)(msg_t *, void *);
void CEX_msg_write_submsg(msg_t *, submsg_writer, void *);

/* Reading Message */
msg_t *CEX_make_read_msg(array_t *buffer);
#define IS_RMSG(msg) (MSG_MODE(msg)==MSG_R)
#define REQ_RMSG(msg) do {                                \
      if (unlikely(!IS_RMSG(msg))) {                      \
              Fatal("non reading message (mode=%d)",     \
                    MSG_MODE(msg));                      \
      }                                                  \
} while (0)                
static inline char CEX_msg_read_char(msg_t *) GCC_ATTRIBUTE((always_inline));
/* 32-bit integers */
static inline unsigned int CEX_msg_read_uint(msg_t *) GCC_ATTRIBUTE((always_inline));
static inline int CEX_msg_read_int(msg_t *) GCC_ATTRIBUTE((always_inline));
typedef void (*msg_element_reader)(msg_t *, array_t *);
static inline array_t *CEX_msg_read_array(msg_t *, size_t el_size, msg_element_reader) GCC_ATTRIBUTE((always_inline));

array_t *CEX_msg_read_char_array(msg_t *);
array_t *CEX_msg_read_uint_array(msg_t *);
array_t *CEX_msg_read_int_array(msg_t *);

double CEX_msg_read_double(msg_t *);
array_t *CEX_msg_read_double_array(msg_t *);

vec_t CEX_msg_read_vec(msg_t *);
array_t *CEX_msg_read_vec_array(msg_t *);

#define REQ_MSG_EOFP(msg_form) do {                      \
      msg_t *_tmp_msg = (msg_form);                      \
      REQ_RMSG(_tmp_msg);                                     \
      if (unlikely(!MSG_EOFP(_tmp_msg))) {                    \
              Fatal("expected EOFP with %d of %d "            \
                    "character remaining",                         \
                    CEX_msg_len(_tmp_msg) - CEX_msg_tell(_tmp_msg), \
                    CEX_msg_len(_tmp_msg));                        \
      }                                                            \
} while (0)                

#include "msg-inline.h"


#endif /* _MSG_H */
