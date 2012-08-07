#include "../../caldgemm/caldgemm_config_load.h"
HPL_fprintf( TEST->outfp, "Tiling %d x %d", TILING_X, TILING_Y);
#ifdef CALDGEMM_TRANSPOSED_A
HPL_fprintf( TEST->outfp, " - At");
#endif
#ifdef CALDGEMM_TRANSPOSED_B
HPL_fprintf( TEST->outfp, " - Bt");
#endif

#ifdef CALDGEMM_USE_VEC_MEMCPY_PREFETCH
HPL_fprintf( TEST->outfp, " - VEC_MEMCPY_PREFETCH");
#endif

#ifdef CALDGEMM_STREAMING_STORES_DIVIDE
HPL_fprintf( TEST->outfp, " - STREAMING_STORES_DIVIDE");
#endif

#ifdef CALDGEMM_STREAMING_STORES_MERGE
HPL_fprintf( TEST->outfp, " - STREAMING_STORES_MERGE");
#endif

#ifdef CALDGEMM_PREFETCH_MERGE_STORES
HPL_fprintf( TEST->outfp, " - PREFETCH_MERGE_STORES");
#endif

#ifdef CALDGEMM_MERGE_NOPS
HPL_fprintf( TEST->outfp, " - MERGE_NOPS %d", CALDGEMM_MERGE_NOPS);
#endif

#ifdef CALDGEMM_LATE_EXIT_CONDITION
HPL_fprintf( TEST->outfp, " - LATE_EXIT_CONDITION");
#endif

#ifdef CALDGEMM_SHIFT_TEXTURE
HPL_fprintf( TEST->outfp, " - SHIFT_TEXTURE %d" CALDGEMM_SHIFT_TEXTURE);
#endif

#ifdef CALDGEMM_DIVIDE_STATIC_BUFFER
HPL_fprintf( TEST->outfp, " - DIVIDE_STATIC_BUFFER");
#endif

#ifdef CALDGEMM_DIVIDE_BLOCKING
HPL_fprintf( TEST->outfp, " - DIVIDE_BLOCKING %d", CALDGEMM_DIVIDE_BLOCKING);
#endif

#ifdef CALDGEMM_DIVIDE_TRANSPOSE_TWOPHASE
HPL_fprintf( TEST->outfp, " - TRANSPOSE_TWOPHASE %d", CALDGEMM_TRANSPOSE_BLOCKING);
#endif
