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
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"

void Chem2PF_Single(Vector *pf_vector, double *chem_var, ProblemData *problem_data)
{
  Grid          *grid = VectorGrid(pf_vector);
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  Subvector     *subvector;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index, pf_index;
  double *pf_data;

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

    /* RDF: assume resolution is the same in all 3 directions */
    r = SubgridRX(subgrid);

    subvector = VectorSubvector(pf_vector, is);
    pf_data = SubvectorData(subvector);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      pf_index = SubvectorEltIndex(subvector, i, j, k);
      chem_index = (i-ix ) + (j- iy) * (nx) + (k- iz) * (nx) * (ny);
      pf_data[pf_index] = chem_var[chem_index];
    });
  }
}


void PF2Chem_Single(Vector *pf_vector, double *chem_var, ProblemData *problem_data)
{
  Grid          *grid = VectorGrid(pf_vector);
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  Subvector     *subvector;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index, pf_index;
  double *pf_data;

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

    /* RDF: assume resolution is the same in all 3 directions */
    r = SubgridRX(subgrid);

    subvector = VectorSubvector(pf_vector, is);
    pf_data = SubvectorData(subvector);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      pf_index = SubvectorEltIndex(subvector, i, j, k);
      chem_index = (i-ix ) + (j- iy) * (nx) + (k- iz) * (nx) * (ny);
      chem_var[chem_index] = pf_data[pf_index];
    });
  }
}





void Chem2PF_Multi(Vector *pf_vector, double *chem_var, int num_var, ProblemData *problem_data)
{
  Grid          *grid = VectorGrid(pf_vector);
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  Subvector     *subvector;
  int is, var;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index, pf_index;
  double *pf_data;

  gr_domain = ProblemDataGrDomain(problem_data);

  for(var = 0; var < num_var; var++)
  {

    ForSubgridI(is, subgrids)
    {
      subgrid = SubgridArraySubgrid(subgrids, is);
  
      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);
  
      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);
      nz = SubgridNZ(subgrid);
  
      /* RDF: assume resolution is the same in all 3 directions */
      r = SubgridRX(subgrid);
  
      subvector = VectorSubvector(pf_vector, is);
      pf_data = SubvectorData(subvector);
  
      GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
      {
        pf_index = SubvectorEltIndex(subvector, i, j, k);
        chem_index = var + (i-ix ) * num_var 
                         + (j- iy) * (nx) * num_var 
                         + (k- iz) * (nx) * (ny) * num_var;
        pf_data[pf_index] = chem_var[chem_index];
      });
    }
  }
}





int SubgridNumCells(Grid *grid, ProblemData *problem_data)
{
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  int is, var;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int num_cells;

  gr_domain = ProblemDataGrDomain(problem_data);
  num_cells = 0;

  ForSubgridI(is, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, is);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    /* RDF: assume resolution is the same in all 3 directions */
    r = SubgridRX(subgrid);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      num_cells++;
    });
  }
  return num_cells;
}

void AllocateChemCells(AlquimiaDataPF *alquimia_data, Grid *grid, ProblemData *problem_data)
{
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index;

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

    /* RDF: assume resolution is the same in all 3 directions */
    r = SubgridRX(subgrid);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      
      chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
     
      AllocateAlquimiaState(&alquimia_data->chem_sizes, &alquimia_data->chem_state[chem_index]);
      AllocateAlquimiaProperties(&alquimia_data->chem_sizes, &alquimia_data->chem_properties[chem_index]);
      AllocateAlquimiaAuxiliaryData(&alquimia_data->chem_sizes, &alquimia_data->chem_aux_data[chem_index]);
      AllocateAlquimiaAuxiliaryOutputData(&alquimia_data->chem_sizes, &alquimia_data->chem_aux_output[chem_index]);
    });
  }
}



void AllocatePFChemData(AlquimiaDataPF *alquimia_data, Grid *grid)
{
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index;

  //pH
  alquimia_data->pH = NewVectorType( grid, 1, 0, vector_cell_centered );
  InitVectorAll(alquimia_data->pH, 0.0);
  
  // num_primary - primary activity and free ion concen
  if (alquimia_data->chem_sizes.num_primary > 0)
  {
    alquimia_data->primary_free_ion_concentrationPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_primary);
    alquimia_data->primary_activity_coeffPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_primary);

    for (i = 0; i < alquimia_data->chem_sizes.num_primary; i++)
    {
      alquimia_data->primary_free_ion_concentrationPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->primary_free_ion_concentrationPF[i], 0.0);

      alquimia_data->primary_activity_coeffPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->primary_activity_coeffPF[i], 0.0);
    } 
  }

  // num_sorbed - total_immobile - allocate for num_primary
  if (alquimia_data->chem_sizes.num_sorbed > 0)
  {
    alquimia_data->total_immobilePF = ctalloc(Vector *, alquimia_data->chem_sizes.num_primary);

    for (i = 0; i < alquimia_data->chem_sizes.num_primary; i++)
    {
      alquimia_data->total_immobilePF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->total_immobilePF[i], 0.0);
    } 
  }

  // num_minerals - volfx, SSA, SI, rate
  if (alquimia_data->chem_sizes.num_minerals > 0)
  {
    alquimia_data->mineral_volume_fractionsPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_minerals);
    alquimia_data->mineral_specific_surfacePF = ctalloc(Vector *, alquimia_data->chem_sizes.num_minerals);
    alquimia_data->mineral_saturation_indexPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_minerals);
    alquimia_data->mineral_reaction_ratePF = ctalloc(Vector *, alquimia_data->chem_sizes.num_minerals);

    for (i = 0; i < alquimia_data->chem_sizes.num_minerals; i++)
    {
      alquimia_data->mineral_volume_fractionsPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->mineral_volume_fractionsPF[i], 0.0);

      alquimia_data->mineral_specific_surfacePF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->mineral_specific_surfacePF[i], 0.0);

      alquimia_data->mineral_saturation_indexPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->mineral_saturation_indexPF[i], 0.0);

      alquimia_data->mineral_reaction_ratePF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->mineral_reaction_ratePF[i], 0.0);
    } 
  }

  // num_aqueous_complexes - secondary activity and free ion
  if (alquimia_data->chem_sizes.num_aqueous_complexes > 0)
  {
    alquimia_data->secondary_free_ion_concentrationPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_aqueous_complexes);
    alquimia_data->secondary_activity_coeffPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_aqueous_complexes);

    for (i = 0; i < alquimia_data->chem_sizes.num_aqueous_complexes; i++)
    {
      alquimia_data->secondary_free_ion_concentrationPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->secondary_free_ion_concentrationPF[i], 0.0);

      alquimia_data->secondary_activity_coeffPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->secondary_activity_coeffPF[i], 0.0);
    } 
  }

  // num_surface_sites - site density
  if (alquimia_data->chem_sizes.num_surface_sites > 0)
  {
    alquimia_data->surface_site_densityPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_surface_sites);

    for (i = 0; i < alquimia_data->chem_sizes.num_surface_sites; i++)
    {
      alquimia_data->surface_site_densityPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->surface_site_densityPF[i], 0.0);
    } 
  }

  // num_ion_exchange_sites - cation exchange capacity
  if (alquimia_data->chem_sizes.num_ion_exchange_sites > 0)
  {
    alquimia_data->cation_exchange_capacityPF = ctalloc(Vector *, alquimia_data->chem_sizes.num_ion_exchange_sites);

    for (i = 0; i < alquimia_data->chem_sizes.num_ion_exchange_sites; i++)
    {
      alquimia_data->cation_exchange_capacityPF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->cation_exchange_capacityPF[i], 0.0);
    } 
  }

  // num_aqueous_kinetics - aqueous kinetic rate
  if (alquimia_data->chem_sizes.num_aqueous_kinetics > 0)
  {
    alquimia_data->aqueous_kinetic_ratePF = ctalloc(Vector *, alquimia_data->chem_sizes.num_aqueous_kinetics);

    for (i = 0; i < alquimia_data->chem_sizes.num_aqueous_kinetics; i++)
    {
      alquimia_data->aqueous_kinetic_ratePF[i] = NewVectorType( grid, 1, 0, vector_cell_centered );
      InitVectorAll(alquimia_data->aqueous_kinetic_ratePF[i], 0.0);
    } 
  }
}




void FindIndexFromNameCaseInsensitive(const char* const name,
                               const AlquimiaVectorString* const names,
                               int* index) {
  int i;
  *index = -1;
  for (i = 0; i < names->size; ++i) {
    if (AlquimiaCaseInsensitiveStringCompare(name, names->data[i])) {
      *index = i;
      break;
    }
  }
} 


void CopyChemDataToPF(AlquimiaDataPF *alquimia_data, Vector **concentrations, ProblemData *problem_data)
{
  Grid          *grid = VectorGrid(concentrations[0]);
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index, pf_index1, pf_index2, pf_index3, pf_index4;
  double *pf_dat1, *pf_dat2, *pf_dat3, *pf_dat4;
  Subvector *pf_sub1, *pf_sub2, *pf_sub3, *pf_sub4;


  gr_domain = ProblemDataGrDomain(problem_data);
  int num_primary = alquimia_data->chem_sizes.num_primary;
  int num_minerals = alquimia_data->chem_sizes.num_minerals;
  int num_surf_sites = alquimia_data->chem_sizes.num_surface_sites;
  int num_ion_exchange_sites = alquimia_data->chem_sizes.num_ion_exchange_sites;
  int num_aqueous_kinetics = alquimia_data->chem_sizes.num_aqueous_kinetics;
  int num_aqueous_complexes = alquimia_data->chem_sizes.num_aqueous_complexes;
  int num_sorbed = alquimia_data->chem_sizes.num_sorbed;

    ForSubgridI(is, subgrids)
    {
      subgrid = SubgridArraySubgrid(subgrids, is);
  
      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);
  
      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);
      nz = SubgridNZ(subgrid);
  
      /* RDF: assume resolution is the same in all 3 directions */
      r = SubgridRX(subgrid);

      for(int concen = 0; concen < num_primary; concen++)
      {
        pf_sub1 = VectorSubvector(concentrations[concen],is);
        pf_dat1 = SubvectorData(pf_sub1);

        pf_sub2 = VectorSubvector(alquimia_data->primary_free_ion_concentrationPF[concen],is);
        pf_dat2 = SubvectorData(pf_sub2);

        pf_sub3 = VectorSubvector(alquimia_data->primary_activity_coeffPF[concen],is);
        pf_dat3 = SubvectorData(pf_sub3);
  
        GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
        {
          pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
          pf_index2 = SubvectorEltIndex(pf_sub2, i, j, k);
          pf_index3 = SubvectorEltIndex(pf_sub3, i, j, k);

          chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;

          pf_dat1[pf_index1] = alquimia_data->chem_state[chem_index].total_mobile.data[concen];
          pf_dat2[pf_index2] = alquimia_data->chem_aux_output[chem_index].primary_free_ion_concentration.data[concen];
          pf_dat3[pf_index3] = alquimia_data->chem_aux_output[chem_index].primary_activity_coeff.data[concen];
        });
      }



      for(int mineral = 0; mineral < num_minerals; mineral++)
        {
          pf_sub1 = VectorSubvector(alquimia_data->mineral_volume_fractionsPF[mineral],is);
          pf_dat1 = SubvectorData(pf_sub1);
  
          pf_sub2 = VectorSubvector(alquimia_data->mineral_specific_surfacePF[mineral],is);
          pf_dat2 = SubvectorData(pf_sub2);
  
          pf_sub3 = VectorSubvector(alquimia_data->mineral_saturation_indexPF[mineral],is);
          pf_dat3 = SubvectorData(pf_sub3);
  
          pf_sub4 = VectorSubvector(alquimia_data->mineral_reaction_ratePF[mineral],is);
          pf_dat4 = SubvectorData(pf_sub4);
    
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
          {
            pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
            pf_index2 = SubvectorEltIndex(pf_sub2, i, j, k);
            pf_index3 = SubvectorEltIndex(pf_sub3, i, j, k);
            pf_index4 = SubvectorEltIndex(pf_sub4, i, j, k);
  
            chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
  
            pf_dat1[pf_index1] = alquimia_data->chem_state[chem_index].mineral_volume_fraction.data[mineral];
            pf_dat2[pf_index2] = alquimia_data->chem_state[chem_index].mineral_specific_surface_area.data[mineral];
            pf_dat3[pf_index3] = alquimia_data->chem_aux_output[chem_index].mineral_saturation_index.data[mineral];
            pf_dat4[pf_index4] = alquimia_data->chem_aux_output[chem_index].mineral_reaction_rate.data[mineral];
        });
      }
  
  
  
  
  
  
      for(int surf = 0; surf < num_surf_sites; surf++)
        {
          pf_sub1 = VectorSubvector(alquimia_data->surface_site_densityPF[surf],is);
          pf_dat1 = SubvectorData(pf_sub1);
    
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
          {
            pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
            chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
            pf_dat1[pf_index1] = alquimia_data->chem_state[chem_index].surface_site_density.data[surf];
        });
      }
  
  
  
  
      for(int ion = 0; ion < num_ion_exchange_sites; ion++)
      {
          pf_sub1 = VectorSubvector(alquimia_data->cation_exchange_capacityPF[ion],is);
          pf_dat1 = SubvectorData(pf_sub1);
    
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
          {
            pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
            chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
            pf_dat1[pf_index1] = alquimia_data->chem_state[chem_index].cation_exchange_capacity.data[ion];
        });
      }
  
  
  
  
      for(int rate = 0; rate < num_aqueous_kinetics; rate++)
      {
          pf_sub1 = VectorSubvector(alquimia_data->aqueous_kinetic_ratePF[rate],is);
          pf_dat1 = SubvectorData(pf_sub1);
    
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
          {
            pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
            chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
            pf_dat1[pf_index1] = alquimia_data->chem_aux_output[chem_index].aqueous_kinetic_rate.data[rate];
        });
      }
  
  
  
  
      for(int complex = 0; complex < num_aqueous_complexes; complex++)
      {
          pf_sub1 = VectorSubvector(alquimia_data->secondary_free_ion_concentrationPF[complex],is);
          pf_dat1 = SubvectorData(pf_sub1);
  
          pf_sub2 = VectorSubvector(alquimia_data->secondary_activity_coeffPF[complex],is);
          pf_dat2 = SubvectorData(pf_sub2);
    
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
          {
            pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
            pf_index2 = SubvectorEltIndex(pf_sub2, i, j, k);
            chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
            pf_dat1[pf_index1] = alquimia_data->chem_aux_output[chem_index].secondary_free_ion_concentration.data[complex];
            pf_dat2[pf_index2] = alquimia_data->chem_aux_output[chem_index].secondary_activity_coeff.data[complex];
        });
      }
  
      if (num_sorbed > 0)
      {
        for(int sorbed = 0; sorbed < num_primary; sorbed++)
        {
            pf_sub1 = VectorSubvector(alquimia_data->total_immobilePF[sorbed],is);
            pf_dat1 = SubvectorData(pf_sub1);
      
            GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
            {
              pf_index1 = SubvectorEltIndex(pf_sub1, i, j, k);
              chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;
              pf_dat1[pf_index1] = alquimia_data->chem_state[chem_index].total_immobile.data[sorbed];
          });
        }
      }
    }
}



void CopyPFStateToChem(AlquimiaState* chem_state, AlquimiaSizes chem_sizes, Vector **concentrations, ProblemData *problem_data)
{
  Grid          *grid = VectorGrid(concentrations[0]);
  SubgridArray  *subgrids = GridSubgrids(grid);
  GrGeomSolid   *gr_domain;
  Subgrid       *subgrid;
  Subvector     *concen_sub;
  int is;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int r;
  int chem_index, pf_index;
  int concen;
  double *concen_dat;

  gr_domain = ProblemDataGrDomain(problem_data);
  int num_primary = chem_sizes.num_primary;



    ForSubgridI(is, subgrids)
    {
      subgrid = SubgridArraySubgrid(subgrids, is);
  
      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);
  
      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);
      nz = SubgridNZ(subgrid);
  
      /* RDF: assume resolution is the same in all 3 directions */
      r = SubgridRX(subgrid);

      for(concen = 0; concen < num_primary; concen++)
      {

        concen_sub = VectorSubvector(concentrations[concen],is);
        concen_dat = SubvectorData(concen_sub);
  
        GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
        {
          pf_index = SubvectorEltIndex(concen_sub, i, j, k);
          chem_index = (i-ix) + (j-iy) * nx + (k-iz) * nx * ny;

          chem_state[chem_index].total_mobile.data[concen] = concen_dat[pf_index];
      });
      }
    }
}


void     CopyConcenWithBoundary(
              Vector *x,
              Vector *y)
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




void ProcessGeochemICs(AlquimiaDataPF *alquimia_data, Grid *grid, ProblemData *problem_data, int num_ic_conds, NameArray ic_cond_na, Vector * saturation)
{
  double water_density = 998.0;    // density of water in kg/m**3
  double aqueous_pressure = 101325.0; // pressure in Pa.
  Subgrid       *subgrid;
  Subvector     *chem_ind_sub;
  Subvector     *por_sub, *sat_sub;
  int is = 0;
  int i, j, k;
  int ix, iy, iz;
  int nx, ny, nz;
  int dx, dy, dz;
  int r;
  int chem_index, cond_index, por_index, sat_index;
  double *chem_ind, *por, *sat;
  GrGeomSolid   *gr_domain;
  char* name;
 
  SubgridArray  *subgrids = GridSubgrids(grid);
  gr_domain = ProblemDataGrDomain(problem_data);
  subgrid = SubgridArraySubgrid(subgrids, is);



    // assign interior geochemical conditions
  if (num_ic_conds > 0)
  {
    AllocateAlquimiaGeochemicalConditionVector(num_ic_conds, &alquimia_data->ic_condition_list);
    for (int i = 0; i < num_ic_conds; i++)
    {
      name = NA_IndexToName(ic_cond_na, i);
      AllocateAlquimiaGeochemicalCondition(strlen(name), 0, 0, &alquimia_data->ic_condition_list.data[i]);
      strcpy(alquimia_data->ic_condition_list.data[i].name, name);
    }
  }


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

    chem_ind_sub = VectorSubvector(ProblemDataGeochemCond(problem_data), is);
    chem_ind = SubvectorData(chem_ind_sub);

    por_sub = VectorSubvector(ProblemDataPorosity(problem_data), is);
    por = SubvectorData(por_sub);

    sat_sub = VectorSubvector(saturation, is);
    sat = SubvectorData(sat_sub);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      cond_index = SubvectorEltIndex(chem_ind_sub, i, j, k);
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
      alquimia_data->chem.ProcessCondition(&alquimia_data->chem_engine,
                                    &alquimia_data->ic_condition_list.data[(int)chem_ind[cond_index]], 
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
}




void ProcessGeochemBCs(AlquimiaDataPF *alquimia_data, int num_bc_conds, NameArray bc_cond_na)
{
  char* name;
  double water_density = 998.0;    // density of water in kg/m**3
  double aqueous_pressure = 101325.0; // pressure in Pa.
  // boundary conditions:

    if (num_bc_conds > 0)
  {
    AllocateAlquimiaGeochemicalConditionVector(num_bc_conds, &alquimia_data->bc_condition_list);
    alquimia_data->chem_bc_properties = ctalloc(AlquimiaProperties, num_bc_conds);
    alquimia_data->chem_bc_state = ctalloc(AlquimiaState, num_bc_conds);
    alquimia_data->chem_bc_aux_data = ctalloc(AlquimiaAuxiliaryData, num_bc_conds);
    
    for (int i = 0; i < num_bc_conds; i++)
    {
      name = NA_IndexToName(bc_cond_na, i);
      AllocateAlquimiaGeochemicalCondition(strlen(name), 0, 0, &alquimia_data->bc_condition_list.data[i]);
      strcpy(alquimia_data->bc_condition_list.data[i].name, name);

      AllocateAlquimiaState(&alquimia_data->chem_sizes, &alquimia_data->chem_bc_state[i]);
      AllocateAlquimiaProperties(&alquimia_data->chem_sizes, &alquimia_data->chem_bc_properties[i]);
      AllocateAlquimiaAuxiliaryData(&alquimia_data->chem_sizes, &alquimia_data->chem_bc_aux_data[i]);
    }

    for (int i = 0; i < num_bc_conds; i++)
    {
      alquimia_data->chem_bc_state[i].water_density = water_density;
      alquimia_data->chem_bc_state[i].temperature = 25.0;
      alquimia_data->chem_bc_state[i].porosity = 0.25;
      alquimia_data->chem_bc_state[i].aqueous_pressure = aqueous_pressure;

      alquimia_data->chem.ProcessCondition(&alquimia_data->chem_engine,
                                    &alquimia_data->bc_condition_list.data[i], 
                                    &alquimia_data->chem_bc_properties[i],
                                    &alquimia_data->chem_bc_state[i],
                                    &alquimia_data->chem_bc_aux_data[i],
                                    &alquimia_data->chem_status);
    }
  }
}

