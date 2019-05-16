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
  char*         engine_name;
  char*         chemistry_input_file;
  PFModule      *set_chem_data;
  NameArray     ic_cond_na;
  NameArray     bc_cond_na;
  int           num_ic_conds;
  int           num_bc_conds;
} PublicXtra;

typedef struct {
  Problem       *problem;
  Grid          *grid;
  PFModule      *set_chem_data;
  PFModule      *bc_concentration;
} InstanceXtra;





/*--------------------------------------------------------------------------
 * InitializeChemistry
 *--------------------------------------------------------------------------*/

void InitializeChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *phi)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  Problem       *problem = (instance_xtra->problem);
  Grid          *grid = (instance_xtra->grid);
  Subgrid       *subgrid;

  PFModule      *set_chem_data = (instance_xtra->set_chem_data);
  PFModule      *bc_concentration = (instance_xtra -> bc_concentration);

  GrGeomSolid   *gr_domain;

  char *geochem_conds;
  char* name;
  int num_cells;

  // containers for boundary condition storage
        AlquimiaState *chem_bc_state;
  AlquimiaAuxiliaryData *chem_bc_aux_data;
  AlquimiaProperties *chem_bc_properties;
    bool hands_off = true;

  BeginTiming(public_xtra->time_index);



  gr_domain = ProblemDataGrDomain(problem_data);
  
  // set chem data 
  // this ivokes the geochemcond function and makes available a pfvector of 
  // geochemcond indicators 
  PFModuleInvokeType(SetChemDataInvoke, set_chem_data, (problem_data));
  Vector *geochemcond = ProblemDataGeochemCond(problem_data);

  // find number of active cells for this subgrid
  num_cells = SubgridNumCells(grid, problem_data);
printf("num_cells: %d \n",num_cells);
  // fire alquimia up


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


  alquimia_data->chem.Setup(public_xtra->chemistry_input_file,
                     hands_off,
                     &alquimia_data->chem_engine,
                     &alquimia_data->chem_sizes,
                     &alquimia_data->chem_engine_functionality,
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





  AllocateAlquimiaProblemMetaData(&alquimia_data->chem_sizes, &alquimia_data->chem_metadata);
  alquimia_data->chem_properties = ctalloc(AlquimiaProperties, num_cells);
  alquimia_data->chem_state = ctalloc(AlquimiaState, num_cells);
  alquimia_data->chem_aux_data = ctalloc(AlquimiaAuxiliaryData, num_cells);
  alquimia_data->chem_aux_output = ctalloc(AlquimiaAuxiliaryOutputData, num_cells);


  AllocateChemCells(alquimia_data,grid,problem_data);

  alquimia_data->chem.GetProblemMetaData(&alquimia_data->chem_engine, 
                                         &alquimia_data->chem_metadata, 
                                         &alquimia_data->chem_status);

  printf("here\n");

  if (alquimia_data->chem_status.error != 0) 
  {
    alquimia_error("Alquimia GetProblemMetaData() error: %s", alquimia_data->chem_status.message);
    exit(1);
  }
  else
  {
      if (!amps_Rank(amps_CommWorld))
      {
          amps_Printf("Successful GetProblemMetaData() of Alquimia interface\n");
      }
  }


// assign interior geochemical conditions
  PrintAlquimiaEngineFunctionality(&alquimia_data->chem_engine_functionality, stdout);


  if (public_xtra->num_ic_conds > 0)
  {
    AllocateAlquimiaGeochemicalConditionVector(public_xtra->num_ic_conds, &alquimia_data->ic_condition_list);
    for (int i = 0; i < public_xtra->num_ic_conds; i++)
    {
      name = NA_IndexToName(public_xtra->ic_cond_na, i);
      printf("NAME: %s \n", name);
      AllocateAlquimiaGeochemicalCondition(strlen(name), 0, 0, &alquimia_data->ic_condition_list.data[i]);
      strcpy(alquimia_data->ic_condition_list.data[i].name, name);
            printf("top: %s \n",alquimia_data->ic_condition_list.data[i].name);

    }
  }

 static const double water_density = 999.9720;    // density of water in kg/m**3
 static const double aqueous_pressure = 201325.0; // pressure in Pa.
 double volume = 1.0;
/*
      alquimia_data->chem_properties[i].volume = volume;
      alquimia_data->chem_properties[i].saturation = driver->saturation;

      // Set the thermodynamic state.
      alquimia_data->chem_state[i].water_density = water_density;
      alquimia_data->chem_state[i].temperature = driver->temperature;
      alquimia_data->chem_state[i].porosity = driver->porosity;
      alquimia_data->chem_state[i].aqueous_pressure = aqueous_pressure;
*/


//AlquimiaSizes *chem_sizes = &alquimia_data->chem_sizes;
//printf("chem_sizes.num_primary: %d\n",chem_sizes->num_primary);

 // Grid          *grid = VectorGrid(geochemcond);


    SubgridArray  *subgrids = GridSubgrids(grid);
  
  //Subgrid       *subgrid;
  Subvector     *chem_ind_sub;
  Subvector     *por_sub;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int dx, dy, dz;
  int r;
  int chem_index, pf_index;
  double *chem_ind, *por;
 
    subgrid = SubgridArraySubgrid(subgrids, is);


 

  ForSubgridI(is, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, is);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    dx = SubgridDX(subgrid);
    dy = SubgridDY(subgrid);
    dz = SubgridDZ(subgrid);

    double vol = dx * dy * dz;

    // RDF: assume resolution is the same in all 3 directions 
    r = SubgridRX(subgrid);

    chem_ind_sub = VectorSubvector(geochemcond, is);
    chem_ind = SubvectorData(chem_ind_sub);

    por_sub = VectorSubvector(phi, is);
    por = SubvectorData(por_sub);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      pf_index = SubvectorEltIndex(chem_ind_sub, i, j, k);
      chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
      //chem_var[chem_index] = chem_ind[pf_index];
      printf("right before processcondition");


      alquimia_data->chem_properties[chem_index].volume = vol;
      alquimia_data->chem_properties[chem_index].saturation = 1.0;

      // Set the thermodynamic state.
      alquimia_data->chem_state[chem_index].water_density = 998.0;
      alquimia_data->chem_state[chem_index].temperature = 20.0;
      alquimia_data->chem_state[chem_index].porosity = por[pf_index];
      alquimia_data->chem_state[chem_index].aqueous_pressure = aqueous_pressure;

      

      // Invoke the chemical initial condition.
     alquimia_data->chem.ProcessCondition(&alquimia_data->chem_engine,
                                    &alquimia_data->ic_condition_list.data[(int)chem_ind[pf_index]], 
                                    &alquimia_data->chem_properties[chem_index],
                                    &alquimia_data->chem_state[chem_index],
                                    &alquimia_data->chem_aux_data[chem_index],
                                    &alquimia_data->chem_status);
      if (alquimia_data->chem_status.error != 0)
      {
        printf("TransportDriver: initialization error: %s\n", 
               alquimia_data->chem_status.message);
        break;
      }





    });
  }


CopyChemStateToPF(alquimia_data->chem_state,alquimia_data->chem_sizes,concentrations,problem_data);



// boundary conditions:

    if (public_xtra->num_bc_conds > 0)
  {
    AllocateAlquimiaGeochemicalConditionVector(public_xtra->num_bc_conds, &alquimia_data->bc_condition_list);
    chem_bc_properties = ctalloc(AlquimiaProperties, public_xtra->num_bc_conds);
  	chem_bc_state = ctalloc(AlquimiaState, public_xtra->num_bc_conds);
  	chem_bc_aux_data = ctalloc(AlquimiaAuxiliaryData, public_xtra->num_bc_conds);
    for (int i = 0; i < public_xtra->num_bc_conds; i++)
    {
      name = NA_IndexToName(public_xtra->bc_cond_na, i);
      AllocateAlquimiaGeochemicalCondition(strlen(name), 0, 0, &alquimia_data->bc_condition_list.data[i]);
      strcpy(alquimia_data->bc_condition_list.data[i].name, name);



      AllocateAlquimiaState(&alquimia_data->chem_sizes, &chem_bc_state[i]);
      AllocateAlquimiaProperties(&alquimia_data->chem_sizes, &chem_bc_properties[i]);
      AllocateAlquimiaAuxiliaryData(&alquimia_data->chem_sizes, &chem_bc_aux_data[i]);
    }

    for (int i = 0; i < public_xtra->num_bc_conds; i++)
    {
    	chem_bc_state[i].water_density = water_density;
		chem_bc_state[i].temperature = 25.0;
		chem_bc_state[i].porosity = 0.25;
		chem_bc_state[i].aqueous_pressure = aqueous_pressure;

		alquimia_data->chem.ProcessCondition(&alquimia_data->chem_engine,
                                    &alquimia_data->bc_condition_list.data[i], 
                                    &chem_bc_properties[i],
                                    &chem_bc_state[i],
                                    &chem_bc_aux_data[i],
                                    &alquimia_data->chem_status);

    }


     PFModuleInvokeType(BCConcentrationInvoke, bc_concentration, (problem, grid, concentrations, chem_bc_state, gr_domain));


  }





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


      (instance_xtra->bc_concentration) =
      PFModuleNewInstance(ProblemBCConcentration(problem), ());
  }
  else
  {
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

    PFModuleReNewInstanceType(SetChemDataInitInstanceXtraInvoke,
                              (instance_xtra->set_chem_data),
                              (problem, grid));
    PFModuleReNewInstance((instance_xtra -> bc_concentration), ());
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
    PFModuleFreeInstance(instance_xtra -> bc_concentration);

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
  char* ic_cond_names, *bc_cond_names;

  public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->set_chem_data) = PFModuleNewModule(SetChemData, ());
  (public_xtra->time_index) = RegisterTiming("Chemistry Initialization");


  sprintf(key, "Chemistry.Engine");
  public_xtra->engine_name = GetStringDefault(key,"");

  if (strcmp(public_xtra->engine_name,kAlquimiaStringCrunchFlow) != 0 
  	 && strcmp(public_xtra->engine_name,kAlquimiaStringPFloTran) != 0)
  {
    InputError("Error: Invalid value <%s> for key <%s>. Options are 'CrunchFlow' or 'PFloTran'.\n",
               public_xtra->engine_name, key);
  }

  sprintf(key, "Chemistry.InputFile");
  public_xtra->chemistry_input_file = GetString(key);

  
  bc_cond_names = GetStringDefault("BCConcentration.GeochemCondition.Names","");
  public_xtra->bc_cond_na = NA_NewNameArray(bc_cond_names);
  public_xtra->num_bc_conds = NA_Sizeof(public_xtra->bc_cond_na);

  ic_cond_names = GetStringDefault("GeochemCondition.Names","");
  public_xtra->ic_cond_na = NA_NewNameArray(ic_cond_names);
  public_xtra->num_ic_conds = NA_Sizeof(public_xtra->ic_cond_na);

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
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  int max_nx = 1000;
  int max_ny = 1000;

  int sz = 0;

  /* add local TempData size to `sz' */
  sz += (max_nx + 2 + 2) * (max_ny + 2 + 2);
  sz += (max_nx + 2 + 2) * (max_ny + 2 + 2);
  sz += (max_nx + 2 + 2) * (max_ny + 2 + 2) * 3;

  sz += (max_nx + 3 + 3) * (max_ny + 3 + 3);
  sz += (max_nx + 3 + 3) * (max_ny + 3 + 3);
  sz += (max_nx + 3 + 3) * (max_ny + 3 + 3);
  sz += (max_nx + 3 + 3) * (max_ny + 3 + 3);
  sz += (max_nx + 3 + 3) * (max_ny + 3 + 3);
  sz += (max_nx + 3 + 3);
  sz += (max_nx + 3 + 3);
  sz += (max_nx + 3 + 3);
  sz += (max_nx + 3 + 3) * 4;
  sz += (max_ny + 3 + 3) * 4;
  sz += (max_nx + 3 + 3) * 3;
  sz += (max_nx + 3 + 3) * 3; 
  return sz;
}




