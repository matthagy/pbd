
#include "array.h"
#include <assert.h>

static void
dump_iarray(array_t *ar)
{
        CEX_print_int_array(stdout, ar);
}

static void
add_vec(array_t *ar, double x, double y, double z)
{
        vec_t v;

        Vec3_SET(v, x, y, z);
        VARR_APPEND(ar, v);
}

static void
dump_varray(array_t *ar)
{
        CEX_print_vec_array(stdout, ar);
}

int 
main(int argc, char **argv)
{
        array_t *arr, *rarr, *conc, *srt, *indices, *splice;

        arr = CEX_make_arange(32, -23, -7);
        printf("%lu\n", (unsigned long)ARR_LENGTH(arr));
        printf("%lu\n", (unsigned long)ARR_ALLOCED(arr));
        dump_iarray(arr);
        rarr = CEX_reversed_array(arr);
        dump_iarray(rarr);
        conc = CEX_concat_arrays(arr, rarr);
        dump_iarray(conc);
        srt = CEX_sorted_int_array(conc);
        dump_iarray(srt);

        indices = CEX_make_arange(2,4,1);
        IARR_APPEND(indices, 0);
        IARR_APPEND(indices, 8);
        dump_iarray(indices);
        splice = CEX_splice_array_of_indices(srt, indices);
        dump_iarray(splice);
        CEX_remove_array_of_indices(srt, indices);
        dump_iarray(srt);

        int *ip, ic;
        XARR_FOREACH_SLICE(srt, 1, ARR_LENGTH(srt), 2, ip, ic) {
                printf("%d ", *ip);
        }
        printf("\n");

        array_t *pos, *rpos, *rsplc;

        pos = CEX_make_vec_array(0);
        CEX_align_array(pos, sizeof(double));
        add_vec(pos, 0,0,0);
        add_vec(pos, 5,0,3);
        add_vec(pos, 23,-5,34.23);
        rpos = CEX_reversed_array(pos);
        CEX_extend_array(pos, rpos);
        CEX_extend_array(pos, pos);
        dump_varray(pos);
        rsplc = CEX_splice_array_of_indices(pos, indices);
        dump_varray(rsplc);
        CEX_remove_array_of_indices(pos, indices);
        dump_varray(pos);
        CEX_align_array(pos, sizeof(double));
        dump_varray(pos);

        CEX_free_array(pos);
        CEX_free_array(rpos);

        array_t *cr;

        cr = CEX_make_char_array_from_string("some random string");
        CEX_print_char_array(stdout, cr);
        char buffer[ARR_LENGTH(cr)+1];
        CEX_char_array_as_string(cr, buffer);
        printf("value %.200s\n", buffer);
        
        
        CEX_free_array(cr);

        CEX_free_array(arr);
        CEX_free_array(rarr);
        CEX_free_array(conc);
        CEX_free_array(srt);
        CEX_free_array(indices);
        CEX_free_array(splice);


        return 0;
}
