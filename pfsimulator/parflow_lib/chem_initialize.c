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
*
* Module for initializing the geochemical problem.
*
*-----------------------------------------------------------------------------
*
*****************************************************************************/




#include "parflow.h"
#include "pf_alquimia.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  int           time_index;
  int           num_geochem_conds;
  char*         engine_name;
  char*         chemistry_input_file;
  PFModule      *set_chem_data;
} PublicXtra;

typedef struct {
  Problem       *problem;
  Grid          *grid;
  PFModule      *set_chem_data;
} InstanceXtra;





/*--------------------------------------------------------------------------
 * InitializeChemistry
 *--------------------------------------------------------------------------*/

void InitializeChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  Problem       *problem = (instance_xtra->problem);
  Grid          *grid = (instance_xtra->grid);
  Subgrid          *subgrid;

  PFModule     *set_chem_data = (instance_xtra->set_chem_data);

  char *geochem_conds;
  int num_geochem_conds;
  int num_cells;

  BeginTiming(public_xtra->time_index);

  //num_cells = SubgridNX(subgrid) * SubgridNY(subgrid) * SubgridNZ(subgrid);
  
  // set chem data  
  PFModuleInvokeType(SetChemDataInvoke, set_chem_data, (problem_data));


  Vector        *geochemcond = ProblemDataGeochemCond(problem_data);

  AllocateAlquimiaEngineStatus(&alquimia_data->chem_status);
  CreateAlquimiaInterface(public_xtra->engine_name, &alquimia_data->chem, &alquimia_data->chem_status);

  if (alquimia_data->chem_status.error != 0) 
  {
    alquimia_error("Alquimia interface creation error: %s", alquimia_data->chem_status.message);
    exit(1);
  }
  else
  {
      if (!amps_Rank(amps_CommWorld))
      {
  	     amps_Printf("Successful creation of Alquimia interface\n");
      }
  }

  bool hands_off = true;
  AlquimiaEngineFunctionality chem_engine_functionality;
  alquimia_data->chem.Setup(public_xtra->chemistry_input_file,
                     hands_off,
                     &alquimia_data->chem_engine,
                     &alquimia_data->chem_sizes,
                     &chem_engine_functionality,
                     &alquimia_data->chem_status);

  if (alquimia_data->chem_status.error != 0) 
  {
    alquimia_error("Alquimia interface setup error: %s", alquimia_data->chem_status.message);
    exit(1);
  }
  else
  {
      if (!amps_Rank(amps_CommWorld))
      {
          amps_Printf("Successful setup() of Alquimia interface\n");
      }
  }

  num_cells = 100;

  AllocateAlquimiaProblemMetaData(&alquimia_data->chem_sizes, &alquimia_data->chem_metadata);
  alquimia_data->chem_properties = malloc(sizeof(AlquimiaProperties) * num_cells);
  alquimia_data->chem_state = malloc(sizeof(AlquimiaState) * num_cells);
  alquimia_data->chem_aux_data = malloc(sizeof(AlquimiaAuxiliaryData) * num_cells);
  alquimia_data->chem_aux_output = malloc(sizeof(AlquimiaAuxiliaryOutputData) * num_cells);


EndTiming(public_xtra->time_index);
}


/*--------------------------------------------------------------------------
 * InitializeChemistryInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *InitializeChemistryInitInstanceXtra(Problem *problem, Grid *grid)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra  *instance_xtra;


  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    instance_xtra = ctalloc(InstanceXtra, 1);

    (instance_xtra->set_chem_data) =
      PFModuleNewInstanceType(SetChemDataInitInstanceXtraInvoke,
                              (public_xtra->set_chem_data),
                              (problem, grid));
  }
  else
  {
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

    PFModuleReNewInstanceType(SetChemDataInitInstanceXtraInvoke,
                              (instance_xtra->set_chem_data),
                              (problem, grid));
  }

  /*-----------------------------------------------------------------------
   * Initialize data associated with argument `problem'
   *-----------------------------------------------------------------------*/

  if (problem != NULL)
    (instance_xtra->problem) = problem;

  /*-----------------------------------------------------------------------
   * Initialize data associated with argument `grid'
   *-----------------------------------------------------------------------*/

  if (grid != NULL)
  {
    /* free old data */
    if ((instance_xtra->grid) != NULL)
    {
    }

    /* set new data */
    (instance_xtra->grid) = grid;
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}



/*--------------------------------------------------------------------------
 * InitializeChemistryFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  InitializeChemistryFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);


  if (instance_xtra)
  {
    PFModuleFreeInstance((instance_xtra->set_chem_data));
    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * InitializeChemistryNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *InitializeChemistryNewPublicXtra()
{
  PFModule     *this_module = ThisPFModule;
  PublicXtra   *public_xtra;
  char key[IDB_MAX_KEY_LEN];

  public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->set_chem_data) = PFModuleNewModule(SetChemData, ());
  (public_xtra->time_index) = RegisterTiming("Chemistry Initialization");


  sprintf(key, "Chemistry.Engine");
  public_xtra->engine_name = GetStringDefault(key,"");

  if (strcmp(public_xtra->engine_name,kAlquimiaStringCrunchFlow) != 0 
  	 && strcmp(public_xtra->engine_name,kAlquimiaStringPFloTran) != 0)
  {
    InputError("Error: invalid value <%s> for key <%s>. Options are CrunchFlow or PFloTran\n",
               public_xtra->engine_name, key);
  }

  sprintf(key, "Chemistry.InputFile");
  public_xtra->chemistry_input_file = GetString(key);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * InitializeChemistryFreePublicXtra
 *--------------------------------------------------------------------------*/

void InitializeChemistryFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);


  if (public_xtra)
  {
    PFModuleFreeModule(public_xtra->set_chem_data);
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * InitializeChemistrySizeOfTempData
 *--------------------------------------------------------------------------*/

int InitializeChemistrySizeOfTempData()
{  
  return 0;
}




