
#include "debug.h"

inline void 
truncate_array(array_t *ar, size_t length)
{
        if (unlikely(ARR_LENGTH(ar) < length)) {
                Fatal("cannot truncate array of length %lu to %lu",
                      (unsigned long)ARR_LENGTH(ar), 
                      (unsigned long)length);
        }
        ARR_LENGTH(ar) = length;
}

inline void 
clear_array(array_t *ar)
{
        ARR_LENGTH(ar) = 0;
}

