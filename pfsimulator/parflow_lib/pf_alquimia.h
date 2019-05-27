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
*/
#ifndef PF_ALQUIMIA_H
#define PF_ALQUIMIA_H
 /**********************************************************************EHEADER

*****************************************************************************
* Header file to include all Alquimia related functions
*
*-----------------------------------------------------------------------------
*
*****************************************************************************/
#include "alquimia/alquimia_containers.h"
#include "alquimia/alquimia_interface.h"




typedef struct _ChemPrintFlags {
  int print_primary_mobile;
  int silo_primary_mobile;
  int print_mineral_rate;
  int silo_mineral_rate;
  int print_mineral_volfx;
  int silo_mineral_volfx;
  int print_mineral_surfarea;
  int silo_mineral_surfarea;
  int print_surf_dens;
  int silo_surf_dens;
  int print_CEC;
  int silo_CEC;
  int print_pH;
  int silo_pH;
  int print_aqueous_rate;
  int silo_aqueous_rate;
  int print_mineral_SI;
  int silo_mineral_SI;
  int print_primary_freeion;
  int silo_primary_freeion;
  int print_primary_activity;
  int silo_primary_activity;
  int print_secondary_freeion;
  int silo_secondary_freeion;
  int print_secondary_activity;
  int silo_secondary_activity;
  int print_sorbed;
  int silo_sorbed;
} ChemPrintFlags;




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
  AlquimiaState* chem_bc_state;
  AlquimiaAuxiliaryData* chem_bc_aux_data;
  AlquimiaProperties* chem_bc_properties;
  AlquimiaGeochemicalConditionVector ic_condition_list, bc_condition_list;
  
  // engine functionality
  AlquimiaEngineFunctionality chem_engine_functionality;


  // PF Vectors, needed for printing
  //state vars - also concen vector, but keep that outside
  Vector **total_immobilePF; //immobile portion of primary, num_primary or 0
  Vector **mineral_volume_fractionsPF; // num_minerals
  Vector **mineral_specific_surfacePF; // num_minerals
  Vector **surface_site_densityPF; // num_surface_sites
  Vector **cation_exchange_capacityPF; //num_ion_exchange_sites
  
  //aux output data
  Vector *pH; // ony one value per cell
  Vector **aqueous_kinetic_ratePF; // num_aqueous_kinetics 
  Vector **mineral_saturation_indexPF; // num_minerals
  Vector **mineral_reaction_ratePF; // num_minerals
  Vector **primary_free_ion_concentrationPF; // num_primary
  Vector **primary_activity_coeffPF; // num_primary
  Vector **secondary_free_ion_concentrationPF; // num_aqueous_complexes
  Vector **secondary_activity_coeffPF; // num_aqueous_complexes

  //printing flags
  ChemPrintFlags *print_flags;

} AlquimiaDataPF;



void FindIndexFromNameCaseInsensitive(const char* const name,
                                 const AlquimiaVectorString* const names,
                                 int* index);


/* problem_geochem_cond.c */
typedef void (*GeochemCondInvoke) (ProblemData *problem_data, Vector *geochemcond );
void GeochemCond (ProblemData *problem_data, Vector *geochemcond );
PFModule *GeochemCondInitInstanceXtra (void );
void GeochemCondFreeInstanceXtra (void );
PFModule *GeochemCondNewPublicXtra (void );
void GeochemCondPublicXtra (void );
void GeochemCondFreePublicXtra(void);
int GeochemCondSizeOfTempData (void );


/* chem_advance.c*/
typedef void (*AdvanceChemistryInvoke) (ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *saturation, double dt, double t, int *any_file_dumped, int dump_files, int file_number, char* file_prefix);
typedef PFModule *(*AdvanceChemistryInitInstanceXtraType) (Problem *problem, Grid *grid);
void AdvanceChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *saturation, double dt, double t, int *any_file_dumped, int dump_files, int file_number, char* file_prefix);
PFModule *AdvanceChemistryInitInstanceXtra(Problem *problem, Grid *grid);
void AdvanceChemistryFreeInstanceXtra(void);
PFModule *AdvanceChemistryNewPublicXtra(void);
void AdvanceChemistryFreePublicXtra(void);
int AdvanceChemistrySizeOfTempData(void);


/* chem_initialize.c*/
typedef void (*InitializeChemistryInvoke) (ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *saturation, int *any_file_dumped, int dump_files, double t, int file_number, char* file_prefix);
typedef PFModule *(*InitializeChemistryInitInstanceXtraType) (Problem *problem, Grid *grid);
void InitializeChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *saturation, int *any_file_dumped, int dump_files, double t, int file_number, char* file_prefix);
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



int  SubgridNumCells(Grid *grid, ProblemData *problem_data);
void AllocateChemCells(AlquimiaDataPF *alquimia_data, Grid *grid, ProblemData *problem_data);
void ChemDataToPFVectors(AlquimiaDataPF *alquimia_data, Vector **concentrations, ProblemData *problem_data);
void AdvectedPrimaryToChem(AlquimiaState* chem_state, AlquimiaSizes* chem_sizes, Vector **concentrations, ProblemData *problem_data);
void AllocatePFChemData(AlquimiaDataPF *alquimia_data, Grid *grid);
void ProcessGeochemICs(AlquimiaDataPF *alquimia_data, Grid *grid, ProblemData *problem_data, int num_ic_conds, NameArray ic_cond_na, Vector * saturation);
void ProcessGeochemBCs(AlquimiaDataPF *alquimia_data, int num_bc_conds, NameArray bc_cond_na);
void FreeAlquimiaDataPF(AlquimiaDataPF *alquimia_data, Grid *grid, ProblemData *problem_data);


typedef void (*BCConcentrationInvoke) (Problem *problem, Grid *grid, Vector ** concentrations, AlquimiaState *bc_chem_states, GrGeomSolid *gr_domain);

/* problem_bc_phase_saturation.c */
void BCConcentration(Problem *problem, Grid *grid, Vector ** concentrations, AlquimiaState *bc_chem_states, GrGeomSolid *gr_domain);
PFModule *BCConcentrationInitInstanceXtra(void);
void BCConcentrationFreeInstanceXtra(void);
PFModule *BCConcentrationNewPublicXtra(void);
void BCConcentrationFreePublicXtra(void);
int BCConcentrationSizeOfTempData(void);

void CopyConcenWithBoundary(Vector *x, Vector *y);

void PrintChemistryData(ChemPrintFlags *print_flags, AlquimiaSizes *chem_sizes, AlquimiaProblemMetaData *chem_metadata,
                        double t, int file_number, char* file_prefix, int *any_file_dumped, Vector **concentrations, 
                        Vector **total_immobilePF, Vector **mineral_specific_surfacePF, Vector **mineral_volume_fractionsPF, Vector **surface_site_densityPF, 
                        Vector **cation_exchange_capacityPF, Vector *pH, Vector **aqueous_kinetic_ratePF, Vector **mineral_saturation_indexPF, 
                        Vector **mineral_reaction_ratePF, Vector **primary_free_ion_concentrationPF, Vector **primary_activity_coeffPF, 
                        Vector **secondary_free_ion_concentrationPF, Vector **secondary_activity_coeffPF);


void SelectReactTransTimeStep(double max_velocity, double CFL, 
                double PF_dt, double *advect_react_dt, 
                double *num_rt_iterations);


#endif
