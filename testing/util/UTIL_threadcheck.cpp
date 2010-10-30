#include <unistd.h>
#include <pthread.h>
#include <cstdio>

#ifdef HPL_CHECK_MPI_THREADS
static pthread_t checkMpiThread_id = pthread_t();

extern "C" void checkMpiThread_impl( const char *file, int line )
{
   pthread_t id = pthread_self();
   if ( !pthread_equal( id, checkMpiThread_id ) ) {
      char hostname[ 255 ];
      gethostname(hostname, 255);
      fprintf( stderr, "%s MPI ERROR: incorrect calling thread to MPI function at %s:%d\n", hostname, file, line );
   }
}

extern "C" void setMpiThread()
{
   checkMpiThread_id = pthread_self();
}
#endif
