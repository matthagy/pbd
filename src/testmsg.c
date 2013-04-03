
#include "msg.h"
#include <assert.h>
#include <math.h>

#ifdef NDEBUG
#  undef NDEBUG
#endif

#ifdef NDEBUG
#  error "FUCK"
#endif


#ifndef M_PI
#  define M_PI 3.1415926535897932385
#endif

static void	
setup_read(msg_t *msg)
{
        CEX_finalize_write_msg(msg);
        int len = CEX_msg_tell(msg);
        MSG_PTR(msg) = MSG_START(msg);
        MSG_END(msg) = MSG_PTR(msg) + len;
        MSG_MODE(msg) = MSG_R;
}

static void
test_basic()
{
        msg_t *msg;

        msg = CEX_make_write_msg(0);
        REQ_WMSG(msg);
        CEX_msg_write_char(msg, 'a');
        CEX_msg_write_char(msg, 'b');
        CEX_msg_write_char(msg, 'c');
        setup_read(msg);
        assert(CEX_msg_read_char(msg)=='a');
        assert(CEX_msg_read_char(msg)=='b');
        assert(CEX_msg_read_char(msg)=='c');
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
}

#define EXPECT(tp, frmt, ftp, expectation, expr) do {   \
        tp _expectation = (expectation);                \
        tp _value = (expr);                             \
        if (_expectation != _value) {                   \
                Fatal("expected " frmt " got " frmt,    \
                      (ftp)_expectation, (ftp)_value);  \
        }                                               \
} while(0)

#define EXPECT_CHAR(expectation, expr)          \
        EXPECT(char, "%d", int, expectation, expr)

static void
test_many_chars()
{
        msg_t *msg;

        msg = CEX_make_write_msg(0);
        REQ_WMSG(msg);
        for (int cycle=0; cycle<6; cycle++) {
                for (int ci=1<<cycle; ci<256; ci++) {
                        CEX_msg_write_char(msg, (char)ci);
                }
        }
        setup_read(msg);
        for (int cycle=0; cycle<6; cycle++) {
                for (int ci=1<<cycle; ci<256; ci++) {
                        EXPECT_CHAR((char)ci, CEX_msg_read_char(msg));
                }
        }
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
}

#define EXPECT_UINT(expectation, expr)          \
        EXPECT(unsigned int, "%u", unsigned int, expectation, expr)


static void
test_uint()
{
        msg_t *msg;

        msg = CEX_make_write_msg(0);
        CEX_msg_write_uint(msg, 0UL);
        CEX_msg_write_uint(msg, 12UL);
        CEX_msg_write_uint(msg, 365UL);
        CEX_msg_write_uint(msg, 0xC0EDA55UL);
        CEX_msg_write_uint(msg, 0xffffffffUL);
        setup_read(msg);
        EXPECT_UINT(0UL, CEX_msg_read_uint(msg));
        EXPECT_UINT(12UL, CEX_msg_read_uint(msg));
        EXPECT_UINT(365UL, CEX_msg_read_uint(msg));
        EXPECT_UINT(0xC0EDA55UL, CEX_msg_read_uint(msg));
        EXPECT_UINT(0xffffffffUL, CEX_msg_read_uint(msg));
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
}

static void
test_char_array()
{
        const char *test = "this is a test string";
        array_t *arr = CEX_make_char_array_from_string(test);
        msg_t *msg = CEX_make_write_msg(0);
        CEX_msg_write_char_array(msg, arr);
        CEX_free_array(arr);
        setup_read(msg);
        array_t *rarr = CEX_msg_read_char_array(msg);
        char buffer[ARR_LENGTH(rarr)+1];
        CEX_char_array_as_string(rarr, buffer);
        assert(strcmp(test, buffer)==0);
        CEX_free_array(rarr);
}

#define EXPECT_DOUBLE(epsilon, expectation, expr) do {        \
        double _expectation = (expectation);            \
        double _value = (expr);                         \
        double _epsilon = (epsilon);                                \
        if (fabs(_expectation - _value)>_epsilon) {           \
                Fatal("expected %.20e within %.5e "     \
                      "got %.20e",                      \
                      _expectation, _epsilon, _value);        \
        }                                               \
} while(0)



static void
test_double()
{
        msg_t *msg = CEX_make_write_msg(0);
        CEX_msg_write_double(msg, 0);
        CEX_msg_write_double(msg, M_PI);
        CEX_msg_write_double(msg, 32.07);
        setup_read(msg);
        EXPECT_DOUBLE(1e-10, 0, CEX_msg_read_double(msg));
        EXPECT_DOUBLE(1e-10, M_PI, CEX_msg_read_double(msg));
        EXPECT_DOUBLE(1e-10, 32.07, CEX_msg_read_double(msg));
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
}

static void
test_vec()
{
        vec_t vec;
        Vec3_SET(vec, M_PI, 323.3e-10, 0);
        msg_t *msg = CEX_make_write_msg(0);
        CEX_msg_write_vec(msg, vec);
        setup_read(msg);
        vec_t v2 = CEX_msg_read_vec(msg);
        EXPECT_DOUBLE(1e-10, vec.x, v2.x);
        EXPECT_DOUBLE(1e-10, vec.y, v2.y);
        EXPECT_DOUBLE(1e-10, vec.z, v2.z);
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
}

static void
add_vec(array_t *ar, double x, double y, double z)
{
        vec_t v;

        Vec3_SET(v, x, y, z);
        VARR_APPEND(ar, v);
}

static void
test_vec_array()
{
        array_t *pos, *rpos;

        pos = CEX_make_vec_array(0);
        add_vec(pos, 0,0,0);
        add_vec(pos, 5,0,M_PI);
        add_vec(pos, 23e12,-5e-13,34.23);
        rpos = CEX_reversed_array(pos);
        CEX_extend_array(pos, rpos);
        CEX_extend_array(pos, pos);
        CEX_free_array(rpos);
        
        msg_t *msg = CEX_make_write_msg(0);
        CEX_msg_write_vec_array(msg, pos);
        setup_read(msg);
        array_t *pos2 = CEX_msg_read_vec_array(msg);
        REQ_MSG_EOFP(msg);
        CEX_free_msg(msg);
        REQ_VARR(pos2);
        assert(ARR_LENGTH(pos) == ARR_LENGTH(pos2));
        for (int i=0; i <ARR_LENGTH(pos); i++) {
                vec_t v1 = ARR_INDEX_AS(vec_t, pos, i);
                vec_t v2 = ARR_INDEX_AS(vec_t, pos2, i);
#if 0
                xprintf("%d " Vec3_FRMT("%.5e") " " Vec3_FRMT("%.5e"),
                        i, Vec3_ARGS(v1), Vec3_ARGS(v2));
#endif
                EXPECT_DOUBLE(1e-10, v1.x, v2.x);
                EXPECT_DOUBLE(1e-10, v2.y, v2.y);
                EXPECT_DOUBLE(1e-10, v1.z, v2.z);
        }
        CEX_free_array(pos);
        CEX_free_array(pos2);
}

int
main(int argc, char **argv)
{
        test_basic();
        test_many_chars();
        test_uint();
        test_char_array();
        test_double();
        test_vec();
        test_vec_array();
        return 0;
}
