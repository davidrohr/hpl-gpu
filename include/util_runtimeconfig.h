/**
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

#ifndef UTIL_RUNTIMECONFIG_H
#define UTIL_RUNTIMECONFIG_H

#define HPL_NB_MULTIPLIER_MAX 8

#ifdef __cplusplus
extern "C"
{
#endif

struct runtime_config_options
{
    char* paramdefs;
    int warmup;
    int fastrand;
    int disable_lookahead;
    int lookahead2_turnoff;
    int duration_find_helper;
    int caldgemm_async_fact_dgemm;
    int caldgemm_async_fact_first;
    int caldgemm_async_dtrsm;
    int caldgemm_async_fact_dtrsm;
    int hpl_nb_multiplier_count;
    int hpl_nb_multiplier_threshold[HPL_NB_MULTIPLIER_MAX];
    int hpl_nb_multiplier_factor[HPL_NB_MULTIPLIER_MAX];
};

extern struct runtime_config_options global_runtime_config;

#ifdef __cplusplus
}
#endif

#endif


