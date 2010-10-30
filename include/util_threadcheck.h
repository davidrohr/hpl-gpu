
extern void setMpiThread();
extern void checkMpiThread_impl( const char *, int );

#define checkMpiThread \
   checkMpiThread_impl( __FILE__, __LINE__ )
