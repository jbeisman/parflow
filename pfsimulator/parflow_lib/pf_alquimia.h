/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
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

/*****************************************************************************
* Header file to include all Alquimia related functions
*
*-----------------------------------------------------------------------------
*
*****************************************************************************/
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_interface.h"

typedef struct _AlquimiaDataPF {
  // Per-cell chemistry data.
  AlquimiaProperties* chem_properties;
  AlquimiaState* chem_state;
  AlquimiaAuxiliaryData* chem_aux_data;
  AlquimiaAuxiliaryOutputData* chem_aux_output;

  // Chemistry engine -- one of each of these per thread in general.
  AlquimiaInterface chem;
  void* chem_engine;
  AlquimiaEngineStatus chem_status;

  // Chemistry metadata.
  AlquimiaSizes chem_sizes;
  AlquimiaProblemMetaData chem_metadata;

  // Initial and boundary conditions.
  AlquimiaGeochemicalCondition chem_ic;
  AlquimiaGeochemicalCondition chem_left_bc, chem_right_bc;
  AlquimiaState chem_left_state, chem_right_state;
  AlquimiaAuxiliaryData chem_left_aux_data, chem_right_aux_data;

  // Bookkeeping.
  AlquimiaState advected_chem_state;
  AlquimiaAuxiliaryData advected_chem_aux_data;
  double *test_double_ptr;
  double test_double;
  int test_int;
} AlquimiaDataPF;






/* problem_geochem_cond.c */
typedef void (*GeochemCondInvoke) (ProblemData *problem_data, Vector *geochemcond );
void GeochemCond (ProblemData *problem_data, Vector *geochemcond );
PFModule *GeochemCondInitInstanceXtra (void );
void GeochemCondFreeInstanceXtra (void );
PFModule *GeochemCondNewPublicXtra (void );
void GeochemCondPublicXtra (void );
void GeochemCondFreePublicXtra(void);
int GeochemCondSizeOfTempData (void );

PFModule *ChemAdvanceInitInstanceXtra(Problem *problem, Grid *grid);
void ChemAdvanceFreeInstanceXtra(void);
PFModule *ChemAdvanceNewPublicXtra(void);
void ChemAdvanceFreePublicXtra(void);
int ChemAdvanceSizeOfTempData(void);


/* chem_initialize.c*/
typedef void (*InitializeChemistryInvoke) (ProblemData *problem_data, AlquimiaDataPF *alquimia_data);
typedef PFModule *(*InitializeChemistryInitInstanceXtraType) (Problem *problem, Grid *grid);
void InitializeChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data);
PFModule *InitializeChemistryInitInstanceXtra(Problem *problem, Grid *grid);
void InitializeChemistryFreeInstanceXtra(void);
PFModule *InitializeChemistryNewPublicXtra(void);
void InitializeChemistryFreePublicXtra(void);
int InitializeChemistrySizeOfTempData(void);


/* set_chem_data.c */
typedef void (*SetChemDataInvoke) (ProblemData *problem_data);
typedef PFModule *(*SetChemDataInitInstanceXtraInvoke) (Problem *problem, Grid *grid);
void SetChemData(ProblemData *problem_data);
PFModule *SetChemDataInitInstanceXtra(Problem *problem, Grid *grid);
void SetChemDataFreeInstanceXtra(void);
PFModule *SetChemDataNewPublicXtra(void);
void SetChemDataFreePublicXtra(void);
int SetChemDataSizeOfTempData(void);

void Chem2PF_Single(Vector *pf_vector, double *chem_var);
void PF2Chem_Single(Vector *pf_vector, double *chem_var);
void Chem2PF_Multi(Vector *pf_vector, double *chem_var, int num_var);
void PF2Chem_Multi(Vector *pf_vector, double *chem_var, int num_var);


