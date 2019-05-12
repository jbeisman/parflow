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

/***************************************************************************
*
*
*---------------------------------------------------------------------------
*
***************************************************************************/

#include "parflow.h"
#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"

#include <stdio.h>

typedef struct {
  int time_index;
} PublicXtra;

typedef struct {
  /* InitInstanceXtra arguments */
  Problem *problem;
  Grid    *grid;
  double  *temp_data;
} InstanceXtra;






PFModule  *ChemAdvanceInstanceXtra(
                                   Problem *problem,
                                   Grid *   grid,
                                   double * temp_data)
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra;

  SubgridArray *subgrids;

  Subgrid      *subgrid;

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
      // Does something need to be freed here?
    }

    /* set new data */
    (instance_xtra->grid) = grid;
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}



PFModule  *ChemAdvanceNewPublicXtra()
{
  PFModule     *this_module = ThisPFModule;
  PublicXtra   *public_xtra;


  public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->time_index) = RegisterTiming("Geochemical Engine");

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


void  ChemAdvanceFreeInstanceXtra()
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

void ChemAdvanceFreePublicXtra()
{
  PFModule     *this_module = ThisPFModule;
  PublicXtra   *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}

int  ChemAdvanceSizeOfTempData()
{
  return 0;
}
