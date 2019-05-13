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
      
      chem_index = (i-ix ) + (j- iy) * (nx) + (k- iz) * (nx) * (ny);
     
      AllocateAlquimiaState(&alquimia_data->chem_sizes, &alquimia_data->chem_state[chem_index]);
      AllocateAlquimiaProperties(&alquimia_data->chem_sizes, &alquimia_data->chem_properties[chem_index]);
      AllocateAlquimiaAuxiliaryData(&alquimia_data->chem_sizes, &alquimia_data->chem_aux_data[chem_index]);
      AllocateAlquimiaAuxiliaryOutputData(&alquimia_data->chem_sizes, &alquimia_data->chem_aux_output[chem_index]);
    });
  }
}



void FindIndexFromNameCaseInsensitive(const char* const name,
                               const AlquimiaVectorString* const names,
                               int* index) {
  int i;
  *index = -1;
  for (i = 0; i < names->size; ++i) {
    printf ("!!!!!!! i: %d",i);
    printf("alquimianame: %s \n", names->data[i]);
    if (AlquimiaCaseInsensitiveStringCompare(name, names->data[i])) {
      *index = i;
      break;
    }
  }
} 


void CopyChemStateToPF(AlquimiaState* chem_state, AlquimiaSizes chem_sizes, Vector **concentrations, ProblemData *problem_data)
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
          chem_index = concen + (i-ix ) * num_primary 
                         + (j- iy) * (nx) * num_primary 
                         + (k- iz) * (nx) * (ny) * num_primary;
          concen_dat[pf_index] = chem_state[chem_index].total_mobile.data[concen];





      });

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
          chem_index = concen + (i-ix ) * num_primary 
                         + (j- iy) * (nx) * num_primary 
                         + (k- iz) * (nx) * (ny) * num_primary;
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

