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


/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  int           time_index;
  int           num_geochem_conds;
} PublicXtra;

typedef struct {
  Problem       *problem;
  PFModule      *geochemcond;
  Grid          *grid;
  int           site_data_not_formed;
} InstanceXtra;





/*--------------------------------------------------------------------------
 * InitializeChemistry
 *--------------------------------------------------------------------------*/

void          InitializeChemistry(
                             ProblemData *problem_data)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  Problem       *problem = (instance_xtra->problem);
  Grid          *grid = (instance_xtra->grid);

  PFModule      *geochemcond = (instance_xtra->geochemcond);

  char *geochem_conds;
  int num_geochem_conds;




if ((instance_xtra->site_data_not_formed))
  {
        PFModuleInvokeType(GeochemCondInvoke, geochemcond,
                       (problem_data,
                        ProblemDataGeochemCond(problem_data),
                        public_xtra->num_geochem_conds));
        (instance_xtra->site_data_not_formed) = 0;
  }

/*
          ForSubgridI(is, GridSubgrids(grid))
          {
              int         ix,   iy,   iz, j_x, p_x;
              int         i,j,k, CF_index;
              double      *por, *jinit;
              
                           
              subgrid = GridSubgrid(grid, is);
              Subvector  *por_sub, *jinit_sub;
              
              nx = SubgridNX(subgrid);
              ny = SubgridNY(subgrid);
              nz = SubgridNZ(subgrid);
                            
              ix = SubgridIX(subgrid);
              iy = SubgridIY(subgrid);
              iz = SubgridIZ(subgrid);
              
          por_sub   = VectorSubvector(ProblemDataPorosity(problem_data),is);
          jinit_sub = VectorSubvector(ProblemDataCrunchContam(problem_data),is);
              por       = SubvectorData(por_sub);
          jinit     = SubvectorData(jinit_sub);
              p_x = 0;
          j_x = 0;
          
          instance_xtra -> porCRUNCH = ctalloc(double, nx*ny*nz);
          instance_xtra -> jinitCRUNCH = ctalloc(int, nx*ny*nz);
             
              
              for (i = ix; i < ix + nx; i++)
              {
                  for (j = iy; j < iy + ny; j++)
                  {
                      for (k = iz; k < iz + nz; k++)
                      {
                          p_x = SubvectorEltIndex(por_sub, i,j,k);
               j_x = SubvectorEltIndex(jinit_sub, i,j,k);
                         
                          CF_index = (i-ix ) +
                          (j- iy) * (nx)  +
                          (k- iz) * (nx) * (ny);
              
                          instance_xtra -> porCRUNCH[CF_index]   = por[p_x];
               instance_xtra -> jinitCRUNCH[CF_index] = (int)jinit[j_x];
              // printf("i: %d j: %d k: %d  jinitCRUNCH: %d \n",i,j,k, instance_xtra -> jinitCRUNCH[CF_index]);
                          
                      }
                  }
              }
          }

*/

}


/*--------------------------------------------------------------------------
 * InitializeChemistryInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *InitializeChemistryInitInstanceXtra(
                                          Problem *problem,
                                          Grid *   grid)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  /*-----------------------------------------------------------------------
   * Initialize data associated with `problem'
   *-----------------------------------------------------------------------*/

  if (problem != NULL)
    (instance_xtra->problem) = problem;

  /*-----------------------------------------------------------------------
   * Initialize data associated with `grid'
   *-----------------------------------------------------------------------*/

  if (grid != NULL)
    (instance_xtra->grid) = grid;

  /*-----------------------------------------------------------------------
   * Initialize module instances
   *-----------------------------------------------------------------------*/

  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    (instance_xtra -> geochemcond) =                                 
       PFModuleNewInstance(ProblemGeochemCond(problem), ());
  }
  else
  {
    PFModuleReNewInstance((instance_xtra -> geochemcond), ());
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
    PFModuleFreeInstance(instance_xtra -> geochemcond);
    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * InitializeChemistryNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *InitializeChemistryNewPublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;

  public_xtra = ctalloc(PublicXtra, 1);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * InitializeChemistryFreePublicXtra
 *--------------------------------------------------------------------------*/

void  InitializeChemistryFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);


  if (public_xtra)
  {
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * InitializeChemistrySizeOfTempData
 *--------------------------------------------------------------------------*/

int       InitializeChemistrySizeOfTempData()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  
  return 0;
}




