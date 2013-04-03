
#include <assert.h>
#include "debug.h"
#include "compat.h"

/* communcation routines  */
static inline void
comm_send(comm_t *comm, void *data, int len, MPI_Datatype dt)
{
        int res;
        res = MPI_Send(data, len, dt,
                        comm->comm_rank, comm->current_rule->tag,
                        MPI_COMM_WORLD);
        if (unlikely(res!=0)) {
                Fatal("MPI_Send returned %d", res);
        }

}

static inline void
comm_recv(comm_t *comm, void *data, int len, MPI_Datatype dt)
{
        int res;
        MPI_Status status;
        res = MPI_Recv(data, len, dt,
                       comm->comm_rank, comm->current_rule->tag,
                       MPI_COMM_WORLD, &status);
        if (unlikely(res!=0)) {
                Fatal("MPI_Recv returned %d", res);
        }
}

static inline void
comm_send_int(comm_t *comm, int op)
{
        comm_send(comm, &op, 1, MPI_INT);
}

static inline int
comm_recv_int(comm_t *comm)
{
        int data=0;
        comm_recv(comm, &data, 1, MPI_INT);
        return data;
}

static inline int
msg_bytes(comm_t *comm)
{
        int res;
        MPI_Status status;

        res = MPI_Probe(comm->comm_rank, comm->current_rule->tag,
                        MPI_COMM_WORLD, &status);
        if (unlikely(res!=0)) {
                Fatal("MPI_Probe returned %d", res);
        }
        return GET_MPI_STATUS_BYTES(status);
}

static inline void 
comm_send_ints(comm_t *comm, array_t *ints) 
{
        REQ_IARR(ints);        
        comm_send(comm, ARR_DATA_AS(void, ints),
                  ARR_LENGTH(ints), MPI_INT);
}

static inline void 
comm_send_vecs(comm_t *comm, array_t *vecs) 
{
        REQ_VARR(vecs);        
        comm_send(comm, ARR_DATA_AS(void, vecs),
                  3*ARR_LENGTH(vecs), MPI_DOUBLE);
}

static inline void 
comm_send_ints_by_index(comm_t *comm, array_t *ints, array_t *indices)
{
        array_t *continuous;

        REQ_IARR(ints);
        REQ_IARR(indices);
        continuous = CEX_splice_array_of_indices(ints, indices);
        assert(ARR_LENGTH(continuous)==ARR_LENGTH(indices));
        comm_send_ints(comm, continuous);
        CEX_free_array(continuous);
}

static inline void 
comm_send_vecs_by_index(comm_t *comm, array_t *vecs, array_t *indices)
{
        array_t *continuous;

        REQ_VARR(vecs);
        REQ_IARR(indices);
        continuous = CEX_splice_array_of_indices(vecs, indices);
        assert(ARR_LENGTH(continuous)==ARR_LENGTH(indices));
        comm_send_vecs(comm, continuous);
        CEX_free_array(continuous);
}

static inline void 
comm_recv_ints(comm_t *comm, array_t *dst)
{
        int len;

        REQ_IARR(dst);
        len = MSG_SIZE(comm, int);
        CEX_prealloc_array(dst, len);
        comm_recv(comm, ARR_DATA_AS(void, dst), len, MPI_INT);
        ARR_LENGTH(dst) = len;
}

static inline void 
comm_recv_vecs(comm_t *comm, array_t *dst)
{
        int len;

        REQ_VARR(dst);
        len = MSG_SIZE(comm, vec_t);
        CEX_prealloc_array(dst, len);
        comm_recv(comm, ARR_DATA_AS(void, dst), 3*len, MPI_DOUBLE);
        ARR_LENGTH(dst) = len;
}

static inline void 
comm_recv_extend_ints(comm_t *comm, array_t *dst)
{
        int len;

        REQ_IARR(dst);
        len = MSG_SIZE(comm, int);
        CEX_prealloc_array(dst, ARR_LENGTH(dst) + len);
        comm_recv(comm, ARR_ADDRESS_ELEMENT(dst, ARR_LENGTH(dst)),
                  len, MPI_INT);
        ARR_LENGTH(dst) += len;
}

static inline void 
comm_recv_extend_vecs(comm_t *comm, array_t *dst)
{
        int len;

        REQ_VARR(dst);
        len = MSG_SIZE(comm, vec_t);
        CEX_prealloc_array(dst, ARR_LENGTH(dst) + len);
        comm_recv(comm, ARR_ADDRESS_ELEMENT(dst, ARR_LENGTH(dst)),
                  3*len, MPI_DOUBLE);
        ARR_LENGTH(dst) += len;
}

