/*BHEADER*********************************************************************
 *
 *  This file is part of Parflow. For details, see
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/


#ifndef CHEM_RESTART_H
#define CHEM_RESTART_H

/* chem_initialize.c*/
typedef void (*RestartChemistryInvoke) (ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, int *any_file_dumped, int dump_files, double t, int file_number, char* file_prefix);
typedef PFModule *(*RestartChemistryInitInstanceXtraType) (Problem *problem, Grid *grid);
void RestartChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, int *any_file_dumped, int dump_files, double t, int file_number, char* file_prefix);
PFModule *RestartChemistryInitInstanceXtra(Problem *problem, Grid *grid);
void RestartChemistryFreeInstanceXtra(void);
PFModule *RestartChemistryNewPublicXtra(void);
void RestartChemistryFreePublicXtra(void);
int RestartChemistrySizeOfTempData(void);

#endif

