#include "../../caldgemm/caldgemm_config_load.h"
sprintf( output_buffer + strlen(output_buffer), "Tiling %d x %d", TILING_X, TILING_Y);
#ifdef CALDGEMM_TRANSPOSED_A
sprintf( output_buffer + strlen(output_buffer), " - At");
#endif
#ifdef CALDGEMM_TRANSPOSED_B
sprintf( output_buffer + strlen(output_buffer), " - Bt");
#endif

#ifdef CALDGEMM_USE_VEC_MEMCPY_PREFETCH
sprintf( output_buffer + strlen(output_buffer), " - VEC_MEMCPY_PREFETCH");
#endif

#ifdef CALDGEMM_STREAMING_STORES_DIVIDE
sprintf( output_buffer + strlen(output_buffer), " - STREAMING_STORES_DIVIDE");
#endif

#ifdef CALDGEMM_STREAMING_STORES_MERGE
sprintf( output_buffer + strlen(output_buffer), " - STREAMING_STORES_MERGE");
#endif

#ifdef CALDGEMM_PREFETCH_MERGE_STORES
sprintf( output_buffer + strlen(output_buffer), " - PREFETCH_MERGE_STORES");
#endif

#ifdef CALDGEMM_MERGE_NOPS
sprintf( output_buffer + strlen(output_buffer), " - MERGE_NOPS %d", CALDGEMM_MERGE_NOPS);
#endif

#ifdef CALDGEMM_LATE_EXIT_CONDITION
sprintf( output_buffer + strlen(output_buffer), " - LATE_EXIT_CONDITION");
#endif

#ifdef CALDGEMM_SHIFT_TEXTURE
sprintf( output_buffer + strlen(output_buffer), " - SHIFT_TEXTURE %d", CALDGEMM_SHIFT_TEXTURE);
#endif

#ifdef CALDGEMM_DIVIDE_STATIC_BUFFER
sprintf( output_buffer + strlen(output_buffer), " - DIVIDE_STATIC_BUFFER");
#endif

#ifdef CALDGEMM_DIVIDE_BLOCKING
sprintf( output_buffer + strlen(output_buffer), " - DIVIDE_BLOCKING %d", CALDGEMM_DIVIDE_BLOCKING);
#endif

#ifdef CALDGEMM_DIVIDE_TRANSPOSE_TWOPHASE
sprintf( output_buffer + strlen(output_buffer), " - TRANSPOSE_TWOPHASE %d", CALDGEMM_TRANSPOSE_BLOCKING);
#endif
