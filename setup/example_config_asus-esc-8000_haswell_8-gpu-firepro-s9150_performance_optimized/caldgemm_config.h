/**
 * Compile time configuration of the CALDGEMM library.
 *
 * Copyright 2010:
 *  - David Rohr (drohr@jwdt.org)
 *  - Matthias Bach (bach@compeng.uni-frankfurt.de)
 *  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 *
 * This file is part of CALDGEMM.
 *
 * CALDGEMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CALDGEMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CALDGEMM.  If not, see <http://www.gnu.org/licenses/>.
 */

//CAL DGEMM Kernel Settings
#define CALDGEMM_TRANSPOSED_A					//Use Kernel for transposed A Matrix
//#define CALDGEMM_TRANSPOSED_B					//Use Kernel for transposed B Matrix
//#define CALDGEMM_88							//8x8 tiling (implies memexport)
//#define CALDGEMM_84							//8x4 tiling (implies memexport)
//#define CALDGEMM_48							//4x8 tiling (implies memexport)
#define CALDGEMM_44								//4x4 tiling
//#define CALDGEMM_USE_MEMEXPORT				//Use Memexport for output instead of color buffers
//#define CALDGEMM_COMPUTE_SHADER 64			//Use compute shader, define compute group size
//#define CALDGEMM_DIAGONAL_TEXTURE				//Alternate storage format, only valid for 4x4 kernel, obsolete
#define CALDGEMM_DUAL_ENTRY						//Unroll factor of 2 for 4x4 tiling
//#define CALDGEMM_SINGLE_BUFFER				//Use a single buffer, 4x4 tiling a transposed, experimental
//#define CALDGEMM_SINGLE_BUFFER_IMPROVED		//Alternative access scheme for single buffer, experimental
//#define CALDGEMM_DUAL_BUFFERS					//Double number of buffers, 4x4 tiling a transposed, experimental
#define CALDGEMM_LATE_EXIT_CONDITION			//Put exit condition at end of while loop
#define CALDGEMM_SHIFT_TEXTURE 1				//Shift even numbered rows in texture by n pixels
//#define CALDGEMM_44_BT_64						//64 bit DMA transfers for 4x4 B transposed kernel
//#define CALDGEMM_44_BT_64_CONVERT				//Perform 64 bit DMA transfer but transform to 128 bit for kernel input

//Other Settings
//#define TESTMODE								//Activate Test Mode for debugging
//#define CALDGEMM_LOOP_DETECTION				//Enable loop detection
//#define TEST_KERNEL
//#define TEST_PARAMETERS
//#define CALDGEMM_UNALIGNED_ADDRESSES

#ifndef STD_OUT
#define STD_OUT stdout							//Output for all messages
#endif

#define CALDGEMM_OUTPUT_THREADS 1				//Number of Output threads
#define CALDGEMM_OUTPUT_THREADS_SLOW 2			//Number of output threads when KeepBuffersMapped = false
#define CALDGEMM_EXTRA_OUTPUT_THREADS_LINPACK 0	//Number of additional output threads when running in linpack mode
#define REUSE_BBUFFERS							//Allocate many BBuffers on the GPU so B is not necessarily retransferred, used for A as well
//#define WASTE_MEMORY							//Allocate extra memory before and after every memory segment allocated
//#define CALDGEMM_BENCHMARK_KERNEL 1

//#define DEBUG_MSG_ALLOCATION					//Debug Messages considering GPU buffer allocation when in Debug = true
//#define DEBUG_MSG_TIMED						//Add timestamps to all messages

//#define CALDGEMM_SGEMM						//Experimental SGEMM implementation (requires MemExport)
//#define CALDGEMM_IGEMM						//Experimental IGEMM implementation (Integer instead of single) (requires SGEMM)
//#define CALDGEMM_BGEMM						//Experimental

#define CALDGEMM_MIN_TILE_DIM 32                              //Tile Dimension must be multiple of this
#define CALDGEMM_MIN_TILE_DIM2 128                            //Min dimension of a tile
#define CALDGEMM_MIN_CORRECTION_SIZE 768                      //Min tile size used to calculate correction ratio for tile distribution

//#define CALDGEMM_FORCE_K 16					//Force K Parameter to simulate different kernel perfoemance

#define _NO_AMD_CPU								//Set to run on CPU without 3dnow (which nowadays also include AMD CPUs)
#define _NO_AVX									//Do not use AVX instructions (Only relevant for OpenCL code atm)
#define _NO_ADL									//Do not use ADL library to read GPU temps
//#define _NO_AFFINITY							//Disable affinity setting
//#define USE_OLD_HUGE_MALLOC					//Use old method to allocate huge tables
//#define VTRACE

#define CALDGEMM_USE_VEC_MEMCPY_PREFETCH		//Use prefetching in Divide / Merge Buffer
#define CALDGEMM_STREAMING_STORES_DIVIDE		//Use streaming stores in Divide Buffer
#define CALDGEMM_STREAMING_STORES_MERGE			//Use streaming stores in Merge buffer
#define CALDGEMM_PREFETCH_MERGE_STORES			//Use prefetching in Merge buffer even when using streaming stores
#define CALDGEMM_MERGE_NOPS 20					//Add nops to slow down merge process freeing resources for other tasks
//#define CALDGEMM_MERGE_FLUSH				

//#define CALDGEMM_LDAB_INC 1					//Inc for LDA and LDB to avoid bank conflics
//#define CALDGEMM_LDB_INC 0					//Override LDAB_INC for LDB
//#define CALDGEMM_LDC_INC 0					//see above

//#define CALDGEMM_DIVIDE_STATIC_BUFFER			//Allocate tmpBuffer for divide staticly once and for all
#define CALDGEMM_DIVIDE_BLOCKING 128			//Blocking size for divideBuffer with SHIFT_TEXTURE = 1 (larger multiple of two)
//#define CALDGEMM_DIVIDE_TRANSPOSE_TWOPHASE	//Perform dividebuffer transposition in two phases such that fewer write combining buffers are used, Only works for 2 input buffers per matrix with A transposed!
#define CALDGEMM_TRANSPOSE_BLOCKING 8			//Blocking factor for the transposition (multiple of 2)

//#define CALDGEMM_QUERY_ALL_EVENTS				//Query for all events, not only the last one in a queue
//#define CALDGEMM_USE_CAL_WAIT_FOR_EVENTS		//Use different method for queriying CAL events
//#define CALDGEMM_USE_CAL_WAIT_FOR_EVENTS_NO_POLL	//Do not use active wait to reduce CPU utilization

//Settings for integrated OpenCL kernels, 3rd party kernels from -Ol library must override this
#define OCL_TILING_X 4
#define OCL_TILING_Y 4
#define OCL_TILED_KERNEL
#define OCL_USE_SIMPLE_BUFFERS
#define OCL_GROUP_SIZE_X 8
#define OCL_GROUP_SIZE_Y 8

//Custom header files for optimized height parameters
//#define CALDGEMM_CUSTOM_AUTO_HEIGHT "auto_height.h"		//Can define a custom header file that is included in caldgemm, that handles autoheight feature
//#define CALDGEMM_CUSTOM_HEIGHT_MOD "height_mod.h"		//Same for posterior height adoption

//#define CALDGEMM_OPENCL_EMULATE_STRIDED		//Emulate strided transfers in OpenCL via linear transfers
#define CALDGEMM_OPENCL_USE_ORIGINAL_POINTERS	//Use the original pointers returned by clEnqueueMapBuffer for the DMA transfers and supply an origin parameter for the correct offset
#define CALDGEMM_OPENCL_PROFILED_PIPELINE 0     //Use a profiling command queue to get timing information in pipelined runs.
