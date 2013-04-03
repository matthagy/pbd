
/* convert pair-wise forces to a force table for each particle to
 * evaluate forces paritcle-wise instead of pair-wise
 */
static int nlen(int part_i);
static void nappend(int part_i, int neighbor_i);

typedef struct neighbor_t neighbor_t;
struct neighbor_t
{
        int index;
        neighbor_t *next;
};

static array_t *neighbor_heads=NULL;
static array_t *neighbor_block=NULL;
static array_t *neighbor_offsets=NULL;

static void
setup_force_aux(void)
{
        int N_neighbors = ARR_LENGTH(CEX_internal_neighbors) + 
                          (ARR_LENGTH(CEX_external_neighbors) >> 1);
        //        xprintf("N_neighbors=%d particles=%d", N_neighbors, CEX_N_internal_particles);
        if (neighbor_block==NULL) {
                neighbor_block = CEX_make_array(sizeof(neighbor_t), 0);
                CEX_align_array(neighbor_block, sizeof(int));
                neighbor_heads = CEX_make_array(sizeof(neighbor_t*), 0);
                neighbor_offsets = CEX_make_array(sizeof(int*),0);
                CEX_align_array(neighbor_offsets, sizeof(int));
        }
        CEX_prealloc_array(neighbor_block, N_neighbors);
        CEX_prealloc_array(neighbor_heads, CEX_N_internal_particles);
        ARR_LENGTH(neighbor_block) = N_neighbors;
        ARR_LENGTH(neighbor_heads) = CEX_N_internal_particles;
        CEX_zero_array_elements(neighbor_heads);
        for (int n_counter=ARR_LENGTH(CEX_internal_neighbors) >> 1,
                *n_ptr=ARR_DATA_AS(int, CEX_internal_neighbors);
             n_counter -- > 0;) {
                int part_i = *(n_ptr++);
                int part_j = *(n_ptr++);
                nappend(part_i, part_j);
                nappend(part_j, part_i);
        }
        for (int n_counter=ARR_LENGTH(CEX_external_neighbors) >> 1,
                *n_ptr=ARR_DATA_AS(int, CEX_external_neighbors);
             n_counter -- > 0;) {
                int i_inner = *(n_ptr++);
                int i_ext = *(n_ptr++);
                nappend(i_inner, i_ext);
        }
        CEX_prealloc_array(neighbor_offsets, CEX_N_internal_particles);
        ARR_LENGTH(neighbor_offsets) = CEX_N_internal_particles;
        int totallen=0;
        for (int i=0; i<CEX_N_internal_particles; i++) {
                totallen += nlen(i);
        }
        int* tbl = CEX_aligned_malloc(sizeof(int)*(totallen + CEX_N_internal_particles),
                                      sizeof(int));
        for (int i=0; i<CEX_N_internal_particles; i++) {
                ARR_INDEX_AS(int*, neighbor_offsets, i) = tbl;
                *(tbl++) = nlen(i);
                for (neighbor_t *nb = ARR_INDEX_AS(neighbor_t*, neighbor_heads, i);
                     nb!=NULL; nb = nb->next) {
                        *(tbl++) = nb->index;
                }
        }
}

static void
nappend(int part_i, int neighbor_i)
{
        assert(ARR_LENGTH(neighbor_block));
        neighbor_t *nb = &ARR_INDEX_AS(neighbor_t, neighbor_block, --ARR_LENGTH(neighbor_block));
        nb->next = ARR_INDEX_AS(neighbor_t *, neighbor_heads, part_i);
        nb->index = neighbor_i;
        ARR_INDEX_AS(neighbor_t *, neighbor_heads, part_i) = nb;
}

static int
nlen(int part_i)
{
        int len=0;
        for (neighbor_t *nb = ARR_INDEX_AS(neighbor_t*, neighbor_heads, part_i);
             nb!=NULL; len++, nb = nb->next);
        return len;
}



static void
evaluate_forces(void)
{
        _SETUP_FORCE_LOCALS
        XBZERO(vec_t, _forces, CEX_N_internal_particles);
        /* loop over all neighbors */
        #pragma omp parallel for schedule(static) \
                firstprivate(r_pair_cutoff_sqr, _box_size, _box_half, _positions, _forces, \
                             _linterp_x_min, _inv_linterp_x_prec, _linterp_table) 
        for (int i=0; i<CEX_N_internal_particles; i++) {
                vec_t pos_i = _positions[i], sum_force={0,0,0};
                int *ptr = ARR_INDEX_AS(int*, neighbor_offsets, i);
                for (int len=*(ptr++); len-- > 0; ) {
                        //xprintf("%d %d", i, len);
                        int neighbor_index = *(ptr++);
                        vec_t r;
                        _PER_SEP(r, pos_i, _positions[neighbor_index]);
                        double rsqr = Vec3_SQR(r);
                        if (unlikely(rsqr<r_pair_cutoff_sqr)) {
                                double force_div_rlen = _INTERPOLATE_FORCE(sqrt(rsqr));
                                vec_t force;
                                Vec3_MUL(force, r, force_div_rlen);
                                Vec3_SUBTO(sum_force, force);
                        }
                }
                _forces[i] = sum_force;
        }
}
