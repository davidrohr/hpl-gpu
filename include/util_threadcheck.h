
#ifdef HPL_CHECK_MPI_THREADS
extern void setMpiThread();
extern void checkMpiThread_impl( const char *, int );

#define checkMpiThread \
   checkMpiThread_impl( __FILE__, __LINE__ )
#else
#define checkMpiThread
#endif
