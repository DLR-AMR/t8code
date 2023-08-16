/* file: memory_usage.h
 *   >> captures the number of bytes used by malloc
 */

// NOTE: LINUS ONLY

/* system header files */
#  include <math.h>
#  include <t8.h>
#  include <stdlib.h>
#  include <malloc.h>

static inline
double
memory_usage(int display){
    /* get malloc info structure */
    struct mallinfo my_mallinfo = mallinfo();

    /*total memory reserved by the system for malloc currently */
    double reserved_mem = (double) my_mallinfo.arena;

    /* get all the memory currently allocated to user by malloc, etc. */
    double used_mem = (double) my_mallinfo.hblkhd
                  + (double) my_mallinfo.usmblks
                  + (double) my_mallinfo.uordblks;

    /* get memory not currently allocated to user but malloc controls */
    double free_mem = (double) my_mallinfo.fsmblks
                  + (double) my_mallinfo.fordblks;

    /* Print out concise malloc info line */
    if (display) {
        t8_global_essentialf("MEMORY USAGE: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
              used_mem / (1024.0 * 1024.0),
              used_mem,
              reserved_mem / (1024.0 * 1024.0),
              free_mem);
    }
    return used_mem;
}