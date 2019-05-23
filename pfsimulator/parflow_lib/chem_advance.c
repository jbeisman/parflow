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
* Module for initializing the geochemical problem. This code reads the PF input 
* file chemistry options, starts the alquimia interface, allocates the alquimia 
* and PF data storage, and processes and assigns geochemical initial and boundary
* conditions. The data is saved in a large struct, AlquimiaDataPF 
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
  int time_index;
  double time_conversion_factor;
} PublicXtra;

typedef struct {
  Problem       *problem;
  Grid          *grid;
} InstanceXtra;





/*--------------------------------------------------------------------------
 * AdvanceChemistry
 *--------------------------------------------------------------------------*/

void AdvanceChemistry(ProblemData *problem_data, AlquimiaDataPF *alquimia_data, Vector **concentrations, Vector *saturation, double dt)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  Problem       *problem = (instance_xtra->problem);
  Grid          *grid = (instance_xtra->grid);
  GrGeomSolid   *gr_domain;

  int num_cells;
  bool hands_off = true;
  double field_sum;
  double dt_seconds;

  BeginTiming(public_xtra->time_index);

  gr_domain = ProblemDataGrDomain(problem_data);


  double water_density = 900.0;    // density of water in kg/m**3
  double aqueous_pressure = 101325.0; // pressure in Pa.
  Subgrid       *subgrid;
  Subvector     *por_sub, *sat_sub;
  int is = 0;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int dx, dy, dz;
  int r;
  int chem_index, por_index, sat_index;
  double *por, *sat;
  char* name;
 
  SubgridArray  *subgrids = GridSubgrids(grid);
  subgrid = SubgridArraySubgrid(subgrids, is);

  dt_seconds = dt * public_xtra->time_conversion_factor;


  AdvectedPrimaryToChem(alquimia_data->chem_state, &alquimia_data->chem_sizes, concentrations, problem_data);


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

    por_sub = VectorSubvector(ProblemDataPorosity(problem_data), is);
    por = SubvectorData(por_sub);

    sat_sub = VectorSubvector(saturation, is);
    sat = SubvectorData(sat_sub);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      por_index = SubvectorEltIndex(por_sub, i, j, k);
      sat_index = SubvectorEltIndex(sat_sub, i, j, k);

      chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;

      alquimia_data->chem_properties[chem_index].volume = vol;
      alquimia_data->chem_properties[chem_index].saturation = sat[sat_index];

      // Set the thermodynamic state.
      alquimia_data->chem_state[chem_index].water_density = water_density;
      alquimia_data->chem_state[chem_index].temperature = 20.0;
      alquimia_data->chem_state[chem_index].porosity = por[por_index];
      alquimia_data->chem_state[chem_index].aqueous_pressure = aqueous_pressure;

      // Invoke the chemical initial condition.
      alquimia_data->chem.ReactionStepOperatorSplit(&alquimia_data->chem_engine,
                                             dt_seconds, &alquimia_data->chem_properties[chem_index],
                                             &alquimia_data->chem_state[chem_index],
                                             &alquimia_data->chem_aux_data[chem_index],
                                             &alquimia_data->chem_status);
      if (alquimia_data->chem_status.error != 0)
      {
        printf("ProcessGeochemICs: initialization error: %s\n", 
               alquimia_data->chem_status.message);
        break;
      }

    });
  }

  ChemDataToPFVectors(alquimia_data,concentrations,problem_data);


  

  EndTiming(public_xtra->time_index);
}


/*--------------------------------------------------------------------------
 * AdvanceChemistryInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *AdvanceChemistryInitInstanceXtra(Problem *problem, Grid *grid)
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra;

  SubgridArray *subgrids;

  Subgrid      *subgrid;

  int max_nx, max_ny, max_nz;
  int sg;


  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

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
 * AdvanceChemistryFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  AdvanceChemistryFreeInstanceXtra()
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * AdvanceChemistryNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *AdvanceChemistryNewPublicXtra()
{
  PFModule     *this_module = ThisPFModule;
  PublicXtra   *public_xtra;
  char key[IDB_MAX_KEY_LEN];
  NameArray      switch_na;
  char*          switch_name;
  int            switch_value;
  switch_na = NA_NewNameArray("s S SECONDS Seconds seconds m M MINUTES Minutes minutes h H HOURS Hours hours d D DAYS Days days y Y YEARS Years years");


  public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->time_index) = RegisterTiming("Chemistry Solver");

  sprintf(key, "Chemistry.ParFlowTimeUnits");
  switch_name = GetString(key);
  switch_value = NA_NameToIndex(switch_na, switch_name);
  if(switch_value < 0)
  {
 InputError("Error: invalid print switch value <%s> for key <%s>. Available options are one of <s S SECONDS Seconds seconds m M MINUTES Minutes minutes h H HOURS Hours hours d D DAYS Days days y Y YEARS Years years>\n",
       switch_name, key );
  }

  if (switch_value < 5)
    {
      //PF in seconds
      public_xtra->time_conversion_factor = 1.0;
    }
    else if (switch_value > 4 && switch_value < 10)
    {
      //PF in minutes
      public_xtra->time_conversion_factor = 60.0;
    }
    else if (switch_value > 9 && switch_value < 15)
    {
      //PF in hours
      public_xtra->time_conversion_factor = 3600.0; 
    }
    else if (switch_value > 14 && switch_value < 20)
    {
      //PF in days
      public_xtra->time_conversion_factor = 3600.0 * 24.0;
    }
    else if (switch_value > 19)
    {
      //PF in years
      public_xtra->time_conversion_factor = 3600.0 * 24.0 * 365.25;
    }

  NA_FreeNameArray(switch_na);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * AdvanceChemistryFreePublicXtra
 *--------------------------------------------------------------------------*/

void AdvanceChemistryFreePublicXtra()
{
  PFModule     *this_module = ThisPFModule;
  PublicXtra   *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * AdvanceChemistrySizeOfTempData
 *--------------------------------------------------------------------------*/

int AdvanceChemistrySizeOfTempData()
{

  int sz = 0;

  return sz;
}




