/**
 * Some hacks to work around features missing in ancient glibc versions.
 *
 * Copyright 2010:
 *  - David Rohr (drohr@jwdt.org)
 *  - Matthias Bach (bach@compeng.uni-frankfurt.de)
 *  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 *
 * This file is part of HPL-GPU.
 *
 * HPL-GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HPL-GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HPL-GPU.  If not, see <http://www.gnu.org/licenses/>.
 *
 * In addition to the rules layed out by the GNU General Public License
 * the following exception is granted:
 *
 * Use with the Original BSD License.
 *
 * Notwithstanding any other provision of the GNU General Public License
 * Version 3, you have permission to link or combine any covered work with
 * a work licensed under the 4-clause BSD license into a single combined
 * work, and to convey the resulting work.  The terms of this License will
 * continue to apply to the part which is the covered work, but the special
 * requirements of the 4-clause BSD license, clause 3, concerning the
 * requirement of acknowledgement in advertising materials will apply to
 * the combination as such.
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

