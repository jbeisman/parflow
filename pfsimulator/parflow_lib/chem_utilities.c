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
 **********************************************************************EHEADER

*****************************************************************************
* Collection of utility functions for reactive transport coupling
*
*-----------------------------------------------------------------------------
*
*****************************************************************************/
#include "parflow.h"
#include "pf_alquimia.h"

#ifdef HAVE_ALQUIMIA
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"
#endif



void CopyConcenWithBoundary(Vector *x, Vector *y)
{
  Grid       *grid = VectorGrid(x);
  Subgrid    *subgrid;

  Subvector  *y_sub;
  Subvector  *x_sub;

  double     *yp, *xp;

  int ix, iy, iz;
  int nx, ny, nz;
  int nx_x, ny_x, nz_x;
  int nx_y, ny_y, nz_y;

  int i_s, i, j, k, i_x, i_y;


  ForSubgridI(i_s, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, i_s);

    ix = SubgridIX(subgrid)-3;
    iy = SubgridIY(subgrid)-3;
    iz = SubgridIZ(subgrid)-3;

    nx = SubgridNX(subgrid)+6;
    ny = SubgridNY(subgrid)+6;
    nz = SubgridNZ(subgrid)+6;

    x_sub = VectorSubvector(x, i_s);
    y_sub = VectorSubvector(y, i_s);

    nx_x = SubvectorNX(x_sub);
    ny_x = SubvectorNY(x_sub);
    nz_x = SubvectorNZ(x_sub);

    nx_y = SubvectorNX(y_sub);
    ny_y = SubvectorNY(y_sub);
    nz_y = SubvectorNZ(y_sub);

    yp = SubvectorElt(y_sub, ix, iy, iz);
    xp = SubvectorElt(x_sub, ix, iy, iz);

    i_x = 0;
    i_y = 0;
    BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
              i_x, nx_x, ny_x, nz_x, 1, 1, 1,
              i_y, nx_y, ny_y, nz_y, 1, 1, 1,
    {
      yp[i_y] = xp[i_x];
    });
  }
}



int SubgridNumCells(Grid *grid, ProblemData *problem_data)
{
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int num_cells = 0;
  int ai = 0;

  gr_domain = ProblemDataGrDomain(problem_data);

  ForSubgridI(is, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, is);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    BoxLoopI1(i, j, k,
            ix, iy, iz, nx, ny, nz,
            ai, nx, ny, nz, 1, 1, 1,
    {
      num_cells++;
    });
  }
  return num_cells;
}




void SelectReactTransTimeStep(double max_velocity, double CFL, 
                double PF_dt, double *advect_react_dt, 
                int *num_rt_iterations)
{
  double cfl_dt;

  cfl_dt = CFL/max_velocity;

  if (PF_dt <= cfl_dt)
  {
    *advect_react_dt = PF_dt;
    *num_rt_iterations = 1;
  }
  else
  {
    *num_rt_iterations = (int)ceil(PF_dt / cfl_dt);
    *advect_react_dt = PF_dt / *num_rt_iterations;
  }

}



//void CutTimeStepandSolveRecursively(AlquimiaInterface chem, AlquimiaState *chem_state, AlquimiaProperties *chem_properties, void *chem_engine, AlquimiaAuxiliaryData *chem_aux_data, AlquimiaEngineStatus *chem_status, double original_dt, int level)

//void CutTimeStepandSolveRecursively(AlquimiaDataPF * alquimia_data, double original_dt, int level, int chem_index)
//{
//  double half_dt;
//  double dt_epsilon;
//  double t;
//  
//  half_dt = 0.5 * original_dt;
//  dt_epsilon = original_dt * 1e-9;
//  t = 0.0;
//
//  if (level < 10)
//  {
//
//    while (fabs(t - original_dt) > dt_epsilon)
//    {
//      // Solve the geochemical system with a dt of half the original_dt
//     // chem.ReactionStepOperatorSplit(&chem_engine,
//     //                                      half_dt, chem_properties,
//     //                                      chem_state,
//     //                                      chem_aux_data,
//     //                                      chem_status);
//
//      alquimia_data->chem.ReactionStepOperatorSplit(&alquimia_data->chem_engine,
//                                             half_dt, &alquimia_data->chem_properties[chem_index],
//                                             &alquimia_data->chem_state[chem_index],
//                                             &alquimia_data->chem_aux_data[chem_index],
//                                             &alquimia_data->chem_status);
//      //if (!chem_status->converged)
//      if (!(alquimia_data->chem_status.converged))
//      {
//        //CutTimeStepandSolveRecursively(chem, chem_state, chem_properties, chem_engine, chem_aux_data, chem_status, half_dt, level+1);
//
//        CopyAlquimiaState(&alquimia_data->chem_state_temp, &alquimia_data->chem_state[chem_index]);
//        CopyAlquimiaAuxiliaryData(&alquimia_data->chem_aux_data_temp ,&alquimia_data->chem_aux_data[chem_index]);
//        CopyAlquimiaProperties(&alquimia_data->chem_properties_temp, &alquimia_data->chem_properties[chem_index]);
//
//
//                CutTimeStepandSolveRecursively(alquimia_data, half_dt, level+1, chem_index);
//
//      }
//
//      CopyAlquimiaState(&alquimia_data->chem_state[chem_index], &alquimia_data->chem_state_temp);
//      CopyAlquimiaAuxiliaryData(&alquimia_data->chem_aux_data[chem_index], &alquimia_data->chem_aux_data_temp);
//      CopyAlquimiaProperties(&alquimia_data->chem_properties[chem_index], &alquimia_data->chem_properties_temp);
//
//
//      t += half_dt;
//    }
//  }
//  else
//  {
//    amps_Printf("Geochemical engine failed to converge. Timestep halved 10 times to %e seconds.\n",half_dt);
//    PARFLOW_ERROR("Geochemical engine error, exiting simulation.\n");
//  }
//}


#ifdef HAVE_ALQUIMIA
void CutTimeStepandSolveSingleCell(AlquimiaInterface chem, AlquimiaState *chem_state, AlquimiaProperties *chem_properties, void *chem_engine, AlquimiaAuxiliaryData *chem_aux_data, AlquimiaEngineStatus *chem_status, double original_dt)
{
  double tenth_dt = 0.1 * original_dt;
  amps_Printf("Original dt of %e seconds failed. Attempting to solve 10 smaller timesteps of %e seconds.\n",original_dt,tenth_dt);

  for (int i = 0; i < 10; i++)
  {
      // reduce dt by factor of 10, subcycle geochem solve 10 timesu
      chem.ReactionStepOperatorSplit(&chem_engine,
                                    tenth_dt, chem_properties,
                                    chem_state,
                                    chem_aux_data,
                                    chem_status);
   if (!chem_status->converged)
    {
      amps_Printf("Geochemical engine failed to converge at reduced timestep of %e seconds.\n",tenth_dt);
      PARFLOW_ERROR("Geochemical engine error, exiting simulation.\n");
    }
  }
}


void WriteChemChkpt(Grid *grid, 
                   ProblemData *problem_data,
                   AlquimiaSizes *chem_sizes,
                   AlquimiaState *chem_state,
                   AlquimiaAuxiliaryData *chem_aux_data,
                   AlquimiaProperties *chem_properties,
                   char *file_prefix,
                   char *file_suffix)
{
  SubgridArray   *subgrids = GridSubgrids(grid);
  Subgrid        *subgrid;
  GrGeomSolid    *gr_domain;

  int num_chem_vars;
  int is;
  int p;
  int num_cells;
  long size;
  double tmp;
  int ai = 0;

  char file_extn[7] = "pfb";
  char filename[255];
  amps_File file;

  BeginTiming(PFBTimingIndex);

  gr_domain = ProblemDataGrDomain(problem_data);
  p = amps_Rank(amps_CommWorld);

  if (p == 0)
    size = 6 * amps_SizeofDouble + 5 * amps_SizeofInt;
  else
    size = 0;


  // how big is the chemistry problem?
  num_chem_vars = chem_sizes->num_primary + // total mobile
                  chem_sizes->num_minerals * 3 + // volfx, ssa, rate constant
                  chem_sizes->num_surface_sites + // surface sites
                  chem_sizes->num_ion_exchange_sites + // ion exchange sites 
                  chem_sizes->num_isotherm_species * 3 + // isotherm, Freundlich, Langmuir
                  chem_sizes->num_aqueous_kinetics + // aqueous rate constant
                  chem_sizes->num_aux_doubles + // aux data doubles
                  chem_sizes->num_aux_integers; // aux num ints

  if (chem_sizes->num_sorbed > 0)
  {
    num_chem_vars += chem_sizes->num_primary; // sorbed immobile primary
  }

  ForSubgridI(is, subgrids)
  {
      size += 9 * amps_SizeofInt;
  }

  num_cells = SubgridNumCells(grid, problem_data);
  size += num_cells * num_chem_vars * amps_SizeofDouble;

  /* open file */
  sprintf(filename, "%s.%s.%s", file_prefix, file_suffix, file_extn);

  if ((file = amps_FFopen(amps_CommWorld, filename, "wb", size)) == NULL)
  {
    amps_Printf("Error: can't open output file %s\n", filename);
    exit(1);
  }

  /* Compute number of patches to write */
  int num_subgrids = GridNumSubgrids(grid);
  {
    amps_Invoice invoice = amps_NewInvoice("%i", &num_subgrids);

    amps_AllReduce(amps_CommWorld, invoice, amps_Add);

    amps_FreeInvoice(invoice);
  }


  if (p == 0)
  {
    amps_WriteDouble(file, &BackgroundX(GlobalsBackground), 1);
    amps_WriteDouble(file, &BackgroundY(GlobalsBackground), 1);
    amps_WriteDouble(file, &BackgroundZ(GlobalsBackground), 1);

#if 0
    amps_WriteInt(file, &BackgroundNX(GlobalsBackground), 1);
    amps_WriteInt(file, &BackgroundNY(GlobalsBackground), 1);
    amps_WriteInt(file, &BackgroundNZ(GlobalsBackground), 1);
#else
    amps_WriteInt(file, &SubgridNX(GridBackground(grid)), 1);
    amps_WriteInt(file, &SubgridNY(GridBackground(grid)), 1);
    amps_WriteInt(file, &SubgridNZ(GridBackground(grid)), 1);
#endif

    amps_WriteDouble(file, &BackgroundDX(GlobalsBackground), 1);
    amps_WriteDouble(file, &BackgroundDY(GlobalsBackground), 1);
    amps_WriteDouble(file, &BackgroundDZ(GlobalsBackground), 1);

    amps_WriteInt(file, &num_subgrids, 1);
    amps_WriteInt(file, &num_chem_vars, 1);
  }


  ForSubgridI(is, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, is);

    int i, j, k;
    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);

    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    amps_WriteInt(file, &ix, 1);
    amps_WriteInt(file, &iy, 1);
    amps_WriteInt(file, &iz, 1);

    amps_WriteInt(file, &nx, 1);
    amps_WriteInt(file, &ny, 1);
    amps_WriteInt(file, &nz, 1);

    amps_WriteInt(file, &SubgridRX(subgrid), 1);
    amps_WriteInt(file, &SubgridRY(subgrid), 1);
    amps_WriteInt(file, &SubgridRZ(subgrid), 1);

    BoxLoopI1(i, j, k,
            ix, iy, iz, nx, ny, nz,
            ai, nx, ny, nz, 1, 1, 1,
    {
      int chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;

      for (int idx = 0; idx < chem_sizes->num_primary; idx++)
      {
        amps_WriteDouble(file, &chem_state[chem_index].total_mobile.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_minerals; idx++)
      {
        amps_WriteDouble(file, &chem_state[chem_index].mineral_volume_fraction.data[idx], 1);
        amps_WriteDouble(file, &chem_state[chem_index].mineral_specific_surface_area.data[idx], 1);
        amps_WriteDouble(file, &chem_properties[chem_index].mineral_rate_cnst.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_surface_sites; idx++)
      {
        amps_WriteDouble(file, &chem_state[chem_index].surface_site_density.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_ion_exchange_sites; idx++)
      {
        amps_WriteDouble(file, &chem_state[chem_index].cation_exchange_capacity.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_isotherm_species; idx++)
      {
        amps_WriteDouble(file, &chem_properties[chem_index].isotherm_kd.data[idx], 1);
        amps_WriteDouble(file, &chem_properties[chem_index].freundlich_n.data[idx], 1);
        amps_WriteDouble(file, &chem_properties[chem_index].langmuir_b.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aqueous_kinetics; idx++)
      {
        amps_WriteDouble(file, &chem_properties[chem_index].aqueous_kinetic_rate_cnst.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aux_doubles; idx++)
      {
        amps_WriteDouble(file, &chem_aux_data[chem_index].aux_doubles.data[idx], 1);
      }


      if (chem_sizes->num_sorbed > 0)
      {
        for (int idx = 0; idx < chem_sizes->num_primary; idx++)
          amps_WriteDouble(file, &chem_state[chem_index].total_immobile.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aux_integers; idx++)
      {
        tmp = (double)chem_aux_data[chem_index].aux_ints.data[idx];
        amps_WriteDouble(file, &tmp, 1);
      }

    });

  }
  amps_FFclose(file);
  EndTiming(PFBTimingIndex);
}



void ReadChemChkpt(Grid *grid, 
                   ProblemData *problem_data,
                   AlquimiaSizes *chem_sizes,
                   AlquimiaState *chem_state,
                   AlquimiaAuxiliaryData *chem_aux_data,
                   AlquimiaProperties *chem_properties,
                   char *filename)
{
  SubgridArray   *subgrids = GridSubgrids(grid);
  Subgrid        *subgrid;
  GrGeomSolid    *gr_domain;

  int is;
  int p, P;
  int num_chars;
  
  amps_File file;
  double X, Y, Z;
  int NX, NY, NZ;
  double DX, DY, DZ;
  int num_vars;
  double tmp;
  int ai = 0;

  BeginTiming(PFBTimingIndex);
  
  gr_domain = ProblemDataGrDomain(problem_data);
  p = amps_Rank(amps_CommWorld);
  P = amps_Size(amps_CommWorld);

  if (((num_chars = strlen(filename)) < 4) ||
      (strcmp(".pfb", &filename[num_chars - 4])))
  {
    amps_Printf("Error: %s is not in pfb format\n", filename);
    exit(1);
  }

  if ((file = amps_FFopen(amps_CommWorld, filename, "rb", 0)) == NULL)
  {
    amps_Printf("Error: can't open input file %s\n", filename);
    exit(1);
  }

  if (p == 0)
  {
    amps_ReadDouble(file, &X, 1);
    amps_ReadDouble(file, &Y, 1);
    amps_ReadDouble(file, &Z, 1);

    amps_ReadInt(file, &NX, 1);
    amps_ReadInt(file, &NY, 1);
    amps_ReadInt(file, &NZ, 1);

    amps_ReadDouble(file, &DX, 1);
    amps_ReadDouble(file, &DY, 1);
    amps_ReadDouble(file, &DZ, 1);

    amps_ReadInt(file, &P, 1);
    amps_ReadInt(file, &num_vars, 1);
  }


  ForSubgridI(is, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, is);

    int ix, iy, iz;
    int nx, ny, nz;
    int rx, ry, rz;
    int i, j, k;

    (void)subgrid;

    amps_ReadInt(file, &ix, 1);
    amps_ReadInt(file, &iy, 1);
    amps_ReadInt(file, &iz, 1);

    amps_ReadInt(file, &nx, 1);
    amps_ReadInt(file, &ny, 1);
    amps_ReadInt(file, &nz, 1);

    amps_ReadInt(file, &rx, 1);
    amps_ReadInt(file, &ry, 1);
    amps_ReadInt(file, &rz, 1);


    BoxLoopI1(i, j, k,
            ix, iy, iz, nx, ny, nz,
            ai, nx, ny, nz, 1, 1, 1,
    {
      int chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;

      for (int idx = 0; idx < chem_sizes->num_primary; idx++)
      {
        amps_ReadDouble(file, &chem_state[chem_index].total_mobile.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_minerals; idx++)
      {
        amps_ReadDouble(file, &chem_state[chem_index].mineral_volume_fraction.data[idx], 1);
        amps_ReadDouble(file, &chem_state[chem_index].mineral_specific_surface_area.data[idx], 1);
        amps_ReadDouble(file, &chem_properties[chem_index].mineral_rate_cnst.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_surface_sites; idx++)
      {
        amps_ReadDouble(file, &chem_state[chem_index].surface_site_density.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_ion_exchange_sites; idx++)
      {
        amps_ReadDouble(file, &chem_state[chem_index].cation_exchange_capacity.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_isotherm_species; idx++)
      {
        amps_ReadDouble(file, &chem_properties[chem_index].isotherm_kd.data[idx], 1);
        amps_ReadDouble(file, &chem_properties[chem_index].freundlich_n.data[idx], 1);
        amps_ReadDouble(file, &chem_properties[chem_index].langmuir_b.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aqueous_kinetics; idx++)
      {
        amps_ReadDouble(file, &chem_properties[chem_index].aqueous_kinetic_rate_cnst.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aux_doubles; idx++)
      {
        amps_ReadDouble(file, &chem_aux_data[chem_index].aux_doubles.data[idx], 1);
      }


      if (chem_sizes->num_sorbed > 0)
      {
        for (int idx = 0; idx < chem_sizes->num_primary; idx++)
          amps_ReadDouble(file, &chem_state[chem_index].total_immobile.data[idx], 1);
      }


      for (int idx = 0; idx < chem_sizes->num_aux_integers; idx++)
      {
        amps_ReadDouble(file, &tmp, 1);
        chem_aux_data[chem_index].aux_ints.data[idx] = (int)tmp;
      }

    });

  }
  amps_FFclose(file);
  EndTiming(PFBTimingIndex);
}
#endif
