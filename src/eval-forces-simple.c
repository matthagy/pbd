
static void
setup_force_aux(void)
{
}

static inline void
evaluate_internal_forces()
{
        _SETUP_FORCE_LOCALS
        XBZERO(vec_t, _forces, CEX_N_internal_particles);
        /* loop over internal neighbors */
        for (int n_counter=ARR_LENGTH(CEX_internal_neighbors) >> 1,
             *n_ptr=ARR_DATA_AS(int, CEX_internal_neighbors);
             n_counter -- > 0;) {
                int part_i = *(n_ptr++);
                int part_j = *(n_ptr++);
                vec_t r;
                _PER_SEP(r, _positions[part_i], _positions[part_j]);
                double rsqr = Vec3_SQR(r);
                assert(part_i < CEX_N_internal_particles);
                assert(part_j < CEX_N_internal_particles);
                if (rsqr<r_pair_cutoff_sqr) {
                        /* XXX Assumes CEX_pair_force has already
                         * been divied by vector length, s.t. this multiplication
                         * also normalizes the force vector to direction unit vector */
                        double force_div_rlen = _INTERPOLATE_FORCE(sqrt(rsqr));
                        vec_t force;
                        Vec3_MUL(force, r, force_div_rlen);
                        /* How can we tell the compiler that these two segments
                         * of memory will never overalp? i.e. part_i != part_j
                         */
                        Vec3_SUBTO(_forces[part_i], force);
                        Vec3_ADDTO(_forces[part_j], force);
                }
        }
}

static inline void
evaluate_external_forces()
{
        _SETUP_FORCE_LOCALS
        int N_positions GCC_ATTRIBUTE((unused)) = ARR_LENGTH(CEX_positions);
        for (int n_counter=ARR_LENGTH(CEX_external_neighbors) >> 1,
             *n_ptr=ARR_DATA_AS(int, CEX_external_neighbors);
             n_counter -- > 0;) {
                int part_i = *(n_ptr++); /* inner particle */
                int part_j = *(n_ptr++); /* extern particle */
                assert(part_i < CEX_N_internal_particles);
                assert(part_j >= CEX_N_internal_particles);
                assert(part_j < N_positions);
                vec_t r;
                _PER_SEP(r, _positions[part_i], _positions[part_j]);
                double rsqr = Vec3_SQR(r);
                if (rsqr<r_pair_cutoff_sqr) {
                        double force_div_rlen = _INTERPOLATE_FORCE(sqrt(rsqr));
                        vec_t force;
                        Vec3_MUL(force, r, force_div_rlen);
                        Vec3_SUBTO(_forces[part_i], force);
                }
        }
}

static void
evaluate_forces(void)
{
        evaluate_internal_forces();
        evaluate_external_forces();
}

