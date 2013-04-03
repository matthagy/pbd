/* -*- Mode: c -*-
 * main.c - Entry point and command loop for simulator
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "debug.h"
#include "compat.h"
#include "mem.h"
#include "comm.h"
#include "periodic.h"
#include "bd.h"
#include "init.h"

/* entry point */
static void main_master(int argc , char **argv);
static void main_slave(void);

int
main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &CEX_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &CEX_size);
        if (IS_MASTER()) {
                main_master(argc, argv);
        } else {
                main_slave();
        }
        MPI_Finalize();
        return 0;
}

/* master/slave common routine */
static int exit_master_command_loop=0;
static int exit_slave_command_loop=0;

static void
exit_main_loop(void)
{
        if (IS_MASTER()) {
                exit_master_command_loop = 1;
        } else {
                exit_slave_command_loop = 1;
        }
}

static void 
msg_send(int rank, msg_t *msg)
{
        //xprintf("sending len=%d", CEX_msg_len(msg));
        int res = MPI_Send(MSG_START(msg), CEX_msg_len(msg),
                           MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        if (unlikely(res!=0)) {
                Fatal("MPI_Send returned %d", res);
        }
}

static void
msg_recv(int rank, msg_t *msg)
{
        MPI_Status status;
        MPI_Probe(rank, 0, MPI_COMM_WORLD, &status);
        int nbytes = GET_MPI_STATUS_BYTES(status);
        //xprintf("recieving len=%d", nbytes);
        CEX_prealloc_msg(msg, nbytes);
        int res = MPI_Recv(MSG_START(msg), nbytes, MPI_CHAR,
                           rank, 0, MPI_COMM_WORLD, &status);
        if (unlikely(res!=0)) {
                Fatal("MPI_Recv returned %d", res);
        }
        MSG_PTR(msg) = MSG_START(msg);
        MSG_END(msg) = MSG_START(msg) + nbytes;
}

static void
prepare_read_msg(msg_t *msg)
{
        MSG_PTR(msg) = MSG_START(msg);
}

static void
prepare_write_msg(msg_t *msg)
{
        MSG_START(msg) = ARR_DATA_AS(char, MSG_BUFFER(msg));
        MSG_PTR(msg) = MSG_START(msg);
        MSG_END(msg) = MSG_START(msg) + ARR_LENGTH(MSG_BUFFER(msg));
}

static void
finalize_write_msg(msg_t *msg)
{
        MSG_END(msg) = MSG_PTR(msg);
}


/* command declarations */
typedef void (*command_func)(msg_t *recv, msg_t *send);
typedef struct {
        const char *name;
        command_func func;
} command_t;
/* master/slave common commands */
static void exit_command(msg_t *recv, msg_t *send);
static void poll_size_command(msg_t *recv, msg_t *send);
static void poll_thread_name_command(msg_t *recv, msg_t *send);
static void set_thread_name_command(msg_t *recv, msg_t *send);
/* slave specific commands */
/* master specific commands */
static void send_msg_command(msg_t *recv, msg_t *send);
static void recv_msg_command(msg_t *recv, msg_t *send);
/* init commands */
static void initialize_system_command(msg_t *recv, msg_t *send);
static void initialize_random_command(msg_t *recv, msg_t *send);
static void initialize_cell_state_command(msg_t *recv, msg_t *send);
static void initialize_cell_comm_command(msg_t *recv, msg_t *send);
static void initialize_cell_junctions_command(msg_t *recv, msg_t *send);
static void thread_update_neighbors_command(msg_t *recv, msg_t *send);
static void thread_update_forces_command(msg_t *recv, msg_t *send);
static void slave_simulation_loop_command(msg_t *recv, msg_t *send);
static void master_simulate_cycles_command(msg_t *recv, msg_t *send);
static void collect_thread_positions_and_tags_command(msg_t *recv, msg_t *send);
static void collect_thread_state_command(msg_t *recv, msg_t *send);

static command_t commands[] = {
        {"exit", &exit_command},
        {"poll_size", &poll_size_command},
        {"poll_thread_name", &poll_thread_name_command},
        {"set_thread_name", &set_thread_name_command},
        {"send_msg", &send_msg_command},
        {"recv_msg", &recv_msg_command},
        {"initialize_system", &initialize_system_command},
        {"initialize_random", &initialize_random_command},
        {"initialize_cell_state", &initialize_cell_state_command},
        {"initialize_cell_comm", &initialize_cell_comm_command},
        {"initialize_cell_junctions", &initialize_cell_junctions_command},
        {"thread_update_neighbors", &thread_update_neighbors_command},
        {"thread_update_forces", &thread_update_forces_command},
        {"slave_simulation_loop", &slave_simulation_loop_command},
        {"master_simulate_cycles", &master_simulate_cycles_command},
        {"collect_thread_positions_and_tags", &collect_thread_positions_and_tags_command},
        {"collect_thread_state", &collect_thread_state_command},
        {NULL, NULL} /* setinel */
};

static void 
perform_command(msg_t *recv, msg_t *send)
{
        prepare_read_msg(recv);
        prepare_write_msg(send);
        array_t *command_name_arr = CEX_msg_read_char_array(recv);
        char command_name_buffer[ARR_LENGTH(command_name_arr)+1];
        CEX_char_array_as_string(command_name_arr, command_name_buffer);
        CEX_free_array(command_name_arr);
        for(command_t *cptr=commands; cptr->name!=NULL; cptr++) {
                if (cptr->name[0]==command_name_buffer[0] &&
                    strcmp(cptr->name, command_name_buffer)==0) {
                        //xprintf("performing");
                        cptr->func(recv, send);
                        //xprintf("finalize");
                        finalize_write_msg(send);
                        //xprintf("performed");
                        return;
                }
        }
        Fatal("unkown command '%.200s'", command_name_buffer);
}

/* slave MPI reciever loop */
static void
main_slave(void)
{
        REQ_SLAVE();
        msg_t *recv = CEX_make_read_msg(CEX_make_char_array(2048));
        msg_t *send = CEX_make_write_msg(2048);
        exit_slave_command_loop = 0;
        for (;!exit_slave_command_loop;) {
                msg_recv(0, recv);
                //xprintf("recieved %d", CEX_msg_len(recv));
                perform_command(recv, send);
                //xprintf("will send %d", CEX_msg_len(send));
                msg_send(0, send);
                MSG_END(send) = MSG_START(send) + ARR_LENGTH(MSG_BUFFER(send));
        }
        CEX_free_msg(recv);
        CEX_free_msg(send);
}


/* master IO loop with fifos to control process */
/* fifo IO declarations */
static void setup_fifos(int argc, char **argv);
static void xread(void *buf, size_t size);
static void xwrite(void *buf, size_t size);
static void xflush(void);
static unsigned int read_uint(void);
static void write_uint(unsigned int value);

static void perform_remote_command(int rank, msg_t *recv, msg_t *send);

static void
main_master(int argc, char **argv)
{
        REQ_MASTER();
        setup_fifos(argc, argv);
        msg_t *recv = CEX_make_read_msg(CEX_make_char_array(2048));
        msg_t *send = CEX_make_write_msg(2048);
        exit_master_command_loop = 0;
        while (!exit_master_command_loop) {
                unsigned int msg_rank = read_uint();
                /* read a message from reading fifo */
                unsigned int msg_len = read_uint();
                //xprintf("msglen %lu for rank %d", msg_len, msg_rank);
                MSG_END(recv) = MSG_START(recv) + ARR_LENGTH(MSG_BUFFER(recv));
                CEX_prealloc_msg(recv, msg_len);
                xread(MSG_START(recv), msg_len);
                MSG_END(recv) = MSG_START(recv) + msg_len;
                if (msg_rank==CEX_rank) {
                        perform_command(recv, send);
                } else {
                        perform_remote_command(msg_rank, recv, send);
                }
                /* write send message to writing fifo */
                //xprintf("fifo write len %d", CEX_msg_len(send));
                write_uint(CEX_msg_len(send));
                xwrite(MSG_START(send), CEX_msg_len(send));
                xflush();
        }
        CEX_free_msg(recv);
        CEX_free_msg(send);
}

static void
validate_slave_rank(int rank)
{
        REQ_MASTER();
        if (rank<=0 || rank>=CEX_size) {
                Fatal("bad slave rank %d; for world size %d",
                      rank, CEX_size);
        }
}

static void
send_remote_command(int rank, msg_t *recv)
{
        REQ_MASTER();
        validate_slave_rank(rank);
        msg_send(rank, recv);
}

static void
recv_remote_command(int rank, msg_t *send)
{
        REQ_MASTER();
        validate_slave_rank(rank);
        prepare_write_msg(send);
        msg_recv(rank, send);
}

static void
perform_remote_command(int rank, msg_t *recv, msg_t *send)
{
        REQ_MASTER();
        send_remote_command(rank, recv);
        recv_remote_command(rank, send);
}


/* fifo IO definitions */
static FILE* reading_fifo;
static FILE* writing_fifo;
static FILE* open_fifo(const char *path, const char *mode);

static void 
setup_fifos(int argc, char **argv)
{
        REQ_MASTER();
        if (argc!=3) {
                Fatal("takes exactly 2 command line arguments; given %d", argc-1);
        }
        xprintf("reading commands from fifo %.200s", argv[1]);
        xprintf("writing results to fifo %.200s", argv[2]);
        reading_fifo = open_fifo(argv[1], "r");
        writing_fifo = open_fifo(argv[2], "w");
}

static FILE*
open_fifo(const char *path, const char *mode)
{
        REQ_MASTER();
        FILE *fp = fopen(path, mode);
        if (fp!=NULL)
                return fp;
        Fatal("failed to open %.200s for mode %.50s; %.200s (errno=%d)",
              path, mode, strerror(errno), errno);
}

static void
xread(void *buf, size_t size)
{
        REQ_MASTER();
        if (fread(buf, 1, size, reading_fifo) != size) {
                if (feof(reading_fifo)) {
                        Fatal("EOF on reading fifo");
                } else {
                        Fatal("io error on reading fifo; %.200s (errno=%d)",
                              strerror(errno), errno);
                }
        }
}

static void
xwrite(void *buf, size_t size)
{
        REQ_MASTER();
        if (fwrite(buf, 1, size, writing_fifo)!=size) {
                Fatal("io error on writing fifo; %.200s (errno=%d)",
                      strerror(errno), errno);
        }
}

static void
xflush(void)
{
        REQ_MASTER();
        if (fflush(writing_fifo)!=0) {
                Fatal("io error flushing fifo; %.200s (errno=%d)",
                      strerror(errno), errno);
        }
}

static unsigned int
read_uint(void)
{
        unsigned char buff[4];
        xread(buff, sizeof(buff));
        return (((unsigned int)buff[0]) << 24) | 
               (((unsigned int)buff[1]) << 16) | 
               (((unsigned int)buff[2]) << 8)  |
                ((unsigned int)buff[3]);
}

static void
write_uint(unsigned int value)
{
        char buff[4];

        buff[0] = (char)((value>>24) & 0xffU);
        buff[1] = (char)((value>>16) & 0xffU);
        buff[2] = (char)((value>>8) & 0xffU);
        buff[3] = (char)((value>>0) & 0xffU);
        xwrite(buff, sizeof(buff));
}

#if __INTEL_COMPILER
/* disable parameter was never refrenced warnining */
#  pragma warning (disable : 869)
#endif

/* common command definitions */
static void
exit_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        exit_main_loop();
}

static void 
poll_size_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        CEX_msg_write_uint(send, (unsigned int)CEX_size);
}

static void
poll_thread_name_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        array_t *name = CEX_make_char_array_from_string(CEX_thread_name);
        CEX_msg_write_char_array(send, name);
        CEX_free_array(name);
}

static void
set_thread_name_command(msg_t *recv, msg_t *send)
{
        array_t *name = CEX_msg_read_char_array(recv);
        REQ_MSG_EOFP(recv);
        static int have_set_thread_name = 0;
        if (have_set_thread_name) {
                CEX_free(CEX_thread_name);
        } else {
                have_set_thread_name = 1;
        }
        CEX_thread_name = XNEW(char, ARR_LENGTH(name)+1);
        CEX_char_array_as_string(name, CEX_thread_name);
        CEX_free_array(name);
}

/* master specific commands to allow asynchronous messaging */
static void
send_msg_command(msg_t *recv, msg_t *send)
{
        REQ_MASTER();
        int send_rank = CEX_msg_read_int(recv);
        int sub_len = CEX_msg_read_uint(recv);
        /* narrow to submsg */
        char *old_start = MSG_START(recv);
        char *old_end = MSG_END(recv);
        MSG_START(recv) = MSG_PTR(recv);
        MSG_END(recv) = MSG_START(recv) + sub_len;
        send_remote_command(send_rank, recv);
        MSG_PTR(recv) += sub_len;
        MSG_START(recv) = old_start;
        MSG_END(recv) = old_end;
        REQ_MSG_EOFP(recv);
}

static void
recv_msg_command(msg_t *recv, msg_t *send)
{
        REQ_MASTER();
        int send_rank = CEX_msg_read_int(recv);
        REQ_MSG_EOFP(recv);
        recv_remote_command(send_rank, send);
        MSG_PTR(send) = MSG_END(send);
}

/* init commands */
static void
initialize_system_command(msg_t *recv, msg_t *send)
{
        CEX_initialize_system(recv);
        REQ_MSG_EOFP(recv);
}

static void
initialize_random_command(msg_t *recv, msg_t *send)
{
        CEX_initialize_random(recv);
        REQ_MSG_EOFP(recv);
}

static void
initialize_cell_state_command(msg_t *recv, msg_t *send)
{
        CEX_initialize_cell_state(recv);
        REQ_MSG_EOFP(recv);
}

static void
initialize_cell_comm_command(msg_t *recv, msg_t *send)
{
        CEX_initialize_cell_comm(recv);
        REQ_MSG_EOFP(recv);
}

static void
initialize_cell_junctions_command(msg_t *recv, msg_t *send)
{
        CEX_initialize_cell_junctions(recv);
        REQ_MSG_EOFP(recv);
}

/* simulation commands */
static void
thread_update_neighbors_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        CEX_thread_update_neighbors();
}

static void
thread_update_forces_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        CEX_thread_update_forces();
}

static void
slave_simulation_loop_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        CEX_slave_simulation_loop();
}

static void
master_simulate_cycles_command(msg_t *recv, msg_t *send)
{
        int cycles = CEX_msg_read_int(recv);
        REQ_MSG_EOFP(recv);
        CEX_master_simulate_cycles(cycles);
}

static void
collect_thread_positions_and_tags_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        REQ_INIT();
#if 0        
        for (int i=0; i<ARR_LENGTH(CEX_positions); i++) {
                vec_t position = ARR_INDEX_AS(vec_t, CEX_positions, i);
                if (unlikely(position.x > CEX_box_size.x ||
                             position.y > CEX_box_size.y ||
                             position.z > CEX_box_size.z)) {
                        Fatal("Bad position " Vec3_FRMT("%.3f") " (nm)",
                               Vec3_ARGS_SCALED(1e9, position));
                }
        }
#endif
        int n_positions = ARR_LENGTH(CEX_positions);
        ARR_LENGTH(CEX_positions) = CEX_N_internal_particles;
        CEX_msg_write_vec_array(send, CEX_positions);
        ARR_LENGTH(CEX_positions) = n_positions;
        CEX_msg_write_int_array(send, CEX_tags);
}

static void
collect_thread_state_command(msg_t *recv, msg_t *send)
{
        REQ_MSG_EOFP(recv);
        REQ_INIT();
        CEX_msg_write_vec_array(send, CEX_positions);
        CEX_msg_write_int_array(send, CEX_tags);
        CEX_msg_write_int_array(send, CEX_internal_neighbors);
        CEX_msg_write_int_array(send, CEX_external_neighbors);
}
