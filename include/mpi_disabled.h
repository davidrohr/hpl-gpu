#ifndef MPI_DISABLED_H
#define MPI_DISABLED_H

#include <time.h>

typedef void* MPI_Comm;
typedef void* MPI_Datatype;
typedef void* MPI_Request;
typedef int MPI_Status;

#define MPI_SUCCESS                   0
#define MPI_ERR_BUFFER                1
#define MPI_ERR_COUNT                 2
#define MPI_ERR_TYPE                  3
#define MPI_ERR_TAG                   4
#define MPI_ERR_COMM                  5
#define MPI_ERR_RANK                  6
#define MPI_ERR_REQUEST               7
#define MPI_ERR_ROOT                  8
#define MPI_ERR_GROUP                 9
#define MPI_ERR_OP                    10
#define MPI_ERR_TOPOLOGY              11
#define MPI_ERR_DIMS                  12
#define MPI_ERR_ARG                   13
#define MPI_ERR_UNKNOWN               14
#define MPI_ERR_TRUNCATE              15
#define MPI_ERR_OTHER                 16
#define MPI_ERR_INTERN                17
#define MPI_ERR_IN_STATUS             18
#define MPI_ERR_PENDING               19
#define MPI_ERR_ACCESS                20
#define MPI_ERR_AMODE                 21
#define MPI_ERR_ASSERT                22
#define MPI_ERR_BAD_FILE              23
#define MPI_ERR_BASE                  24
#define MPI_ERR_CONVERSION            25
#define MPI_ERR_DISP                  26
#define MPI_ERR_DUP_DATAREP           27
#define MPI_ERR_FILE_EXISTS           28
#define MPI_ERR_FILE_IN_USE           29
#define MPI_ERR_FILE                  30
#define MPI_ERR_INFO_KEY              31
#define MPI_ERR_INFO_NOKEY            32
#define MPI_ERR_INFO_VALUE            33
#define MPI_ERR_INFO                  34
#define MPI_ERR_IO                    35
#define MPI_ERR_KEYVAL                36
#define MPI_ERR_LOCKTYPE              37
#define MPI_ERR_NAME                  38
#define MPI_ERR_NO_MEM                39
#define MPI_ERR_NOT_SAME              40
#define MPI_ERR_NO_SPACE              41
#define MPI_ERR_NO_SUCH_FILE          42
#define MPI_ERR_PORT                  43
#define MPI_ERR_QUOTA                 44
#define MPI_ERR_READ_ONLY             45
#define MPI_ERR_RMA_CONFLICT          46
#define MPI_ERR_RMA_SYNC              47
#define MPI_ERR_SERVICE               48
#define MPI_ERR_SIZE                  49
#define MPI_ERR_SPAWN                 50
#define MPI_ERR_UNSUPPORTED_DATAREP   51
#define MPI_ERR_UNSUPPORTED_OPERATION 52
#define MPI_ERR_WIN                   53
#define MPI_ERR_LASTCODE              54
#define MPI_ERR_SYSRESOURCE          -2

#define MPI_BYTE NULL
#define MPI_INT NULL
#define MPI_DOUBLE NULL
#define MPI_FLOAT NULL
#define MPI_UNDEFINED NULL
#define MPI_COMM_NULL NULL
#define MPI_COMM_WORLD NULL
#define MPI_THREAD_SERIALIZED NULL
#define MPI_STATUS_IGNORE NULL

static inline int MPI_Send_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){return(MPI_SUCCESS);}
static inline int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request){return(MPI_SUCCESS);}
static inline int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status){return(MPI_SUCCESS);}
static inline int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){return(MPI_SUCCESS);}
static inline int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status){return(MPI_SUCCESS);}
static inline int MPI_Comm_rank(MPI_Comm comm, int *rank){*rank = 0;return(MPI_SUCCESS);}
static inline int MPI_Comm_size(MPI_Comm comm, int *size){*size = 1;return(MPI_SUCCESS);}
static inline int MPI_Init(int *argc, char ***argv){return(MPI_SUCCESS);}
static inline int MPI_Init_thread(int *argc, char ***argv, int required, int *provided){*provided=required;return(MPI_SUCCESS);}
static inline int MPI_Issend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){return(MPI_SUCCESS);}
static inline int MPI_Wait(MPI_Request *request, MPI_Status *status){return(MPI_SUCCESS);}
static inline int MPI_Finalize(void){return(MPI_SUCCESS);}
static inline int MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){return(MPI_SUCCESS);}
static inline int MPI_Type_free(MPI_Datatype *type){return(MPI_SUCCESS);}
static inline int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype){return(MPI_SUCCESS);}
static inline int MPI_Type_commit(MPI_Datatype *type){return(MPI_SUCCESS);}
static inline int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request){return(MPI_SUCCESS);}
static inline int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){return(MPI_SUCCESS);}
static inline int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm){return(MPI_SUCCESS);}
static inline int MPI_Comm_free(MPI_Comm *comm){return(MPI_SUCCESS);}
static inline int MPI_Abort(MPI_Comm comm, int errorcode){return(MPI_SUCCESS);}
static inline int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){return(MPI_SUCCESS);}
static inline int MPI_Get_count(MPI_Status* status, MPI_Datatype datatype, int* count){*count = 0;return(MPI_SUCCESS);}
static inline double MPI_Wtime(void)
{
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return(1e-9 * (double) t.tv_nsec + (double) t.tv_sec);
}

#endif
