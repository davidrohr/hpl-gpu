/**
 * Some hacks to work around features missing in ancient glibc versions.
 *
 * The source code is property of the Frankfurt Institute for Advanced Studies
 * (FIAS). None of the material may be copied, reproduced, distributed,
 * republished, downloaded, displayed, posted or transmitted in any form or by
 * any means, including, but not limited to, electronic, mechanical,
 * photocopying, recording, or otherwise, without the prior written permission
 * of FIAS.
 *
 * Authors:
 * David Rohr (drohr@jwdt.org)
 * Matthias Bach (bach@compeng.uni-frankfurt.de)
 * Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 */

#ifndef _GLIBC_HACKS_H_
#define _GLIBC_HACKS_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>

#ifndef CPU_COUNT
int CPU_COUNT( cpu_set_t *set )
{
   int i, count = 0;
   for( i = 0; i < CPU_SETSIZE; ++i )
      if( CPU_ISSET( i, set ) )
         ++count;
   return count;
}
#endif /* CPU_COUNT */

#ifndef CPU_XOR
void CPU_XOR( cpu_set_t *destset,
              cpu_set_t *srcset1, cpu_set_t *srcset2 )
{
   int i, tmp1, tmp2;
   for( i = 0; i < CPU_SETSIZE; ++i )
   {
      tmp1 = CPU_ISSET( i, srcset1 );
      tmp2 = CPU_ISSET( i, srcset2 );
      if( ( tmp1 && !tmp2 ) || ( !tmp1 && tmp2 ) ) // XOR( tmp1, tmp2 )
         CPU_SET( i, destset );
      else
         CPU_CLR( i, destset );
   }
}
#endif /* CPU_XOR */

#endif /*_GLIBC_HACKS_H_*/

