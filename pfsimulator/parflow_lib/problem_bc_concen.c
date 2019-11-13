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
*****************************************************************************/

#include "parflow.h"
#include "pf_alquimia.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  int     num_geochemconds;
  int     num_patches;
  int    *geochemcond_indexes;
  int    *input_types;    /* num_patches input types */
  int     num_domain_patches;
  
  void  **data;           /* num_patches pointers to Type structures */

  NameArray geochem_cond_na;
  NameArray bc_patches_na;
} PublicXtra;

typedef void InstanceXtra;

typedef struct {
  int condition;
} Type0;               /* Dirichlet, constant */

typedef struct {
  char  *filename;
} Type1;                      /* .pfb file */

/*--------------------------------------------------------------------------
 * BCConcentration
 *   This routine implements the concentration boundary conditions
 *   (Dirichlet only) by setting the concentration of 3 ghost layers
 *   outside of the boundary.
 *--------------------------------------------------------------------------*/
#ifdef HAVE_ALQUIMIA
void BCConcentration(Problem *problem,
                     Grid *grid,
                     Vector **concentrations,
                     AlquimiaState *chem_bc_state)
{
  PFModule       *this_module = ThisPFModule;
  PublicXtra     *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Type0          *dummy0;
  Type1          *dummy1;

  int            num_domain_patches = (public_xtra->num_domain_patches);
  int            *input_types = (public_xtra->input_types);

  SubgridArray   *subgrids = GridSubgrids(grid);
  Subgrid        *subgrid;
  Subvector      *concen_sub, *subvector;
  double         *concen_dat, *tmpp;
  Vector         *pfb_indicator;

  int condition;
  int nx_v, ny_v, nz_v;
  int ipatch, is, i, j, k;
  int num_concen, itmp;
  int ix,iy,iz,nx,ny,nz;
  int dir[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
  int ci;
  int iv1, iv2, iv3;

  /*-----------------------------------------------------------------------
   * Implement BC's
   *-----------------------------------------------------------------------*/

  num_concen =  ProblemNumContaminants(problem);

  for (ipatch = 0; ipatch < num_domain_patches; ipatch++)
  {
    switch (input_types[ipatch])
    {
      case -1:
      {

        BCConcenCopyPatch(problem, grid, concentrations, ipatch);

        break;
      }
      case 0:
      {
        dummy0 = (Type0*)(public_xtra->data[ipatch]);
        condition = (dummy0->condition);
        
        for (int concen = 0; concen < num_concen; concen++)
        {
          ForSubgridI(is, subgrids)
          {
            subgrid = SubgridArraySubgrid(subgrids, is);

            concen_sub = VectorSubvector(concentrations[concen],is);
            concen_dat = SubvectorData(concen_sub);

            nx_v = SubvectorNX(concen_sub);
            ny_v = SubvectorNY(concen_sub);
            nz_v = SubvectorNZ(concen_sub);

            BCConcenPatchExtent(subgrid,&ix,&iy,&iz,&nx,&ny,&nz,ipatch);

            ci = 0;
            BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                      ci, nx_v, ny_v, nz_v, 1, 1, 1,
            {
              if (BoundaryCell(ipatch,i,j,k))
              {
                iv1 = SubvectorEltIndex(concen_sub, i + dir[ipatch][0], j + dir[ipatch][1], k + dir[ipatch][2]);
                iv2 = SubvectorEltIndex(concen_sub, i + 2*dir[ipatch][0], j + 2*dir[ipatch][1], k + 2*dir[ipatch][2]);
                iv3 = SubvectorEltIndex(concen_sub, i + 3*dir[ipatch][0], j + 3*dir[ipatch][1], k + 3*dir[ipatch][2]);
                concen_dat[iv1] = chem_bc_state[condition].total_mobile.data[concen];
                concen_dat[iv2] = chem_bc_state[condition].total_mobile.data[concen];
                concen_dat[iv3] = chem_bc_state[condition].total_mobile.data[concen];
              }
            });
          }
        }
        break;
      }
      case 1:
      {
        dummy1 = (Type1*)(public_xtra->data[ipatch]);
        pfb_indicator = NewVectorType(grid, 1, 0, vector_cell_centered);
        ReadPFBinary((dummy1->filename), pfb_indicator);

        subvector = VectorSubvector(pfb_indicator, is);
        tmpp = SubvectorData(subvector);
        
        for (int concen = 0; concen < num_concen; concen++)
        {
          ForSubgridI(is, subgrids)
          {
            subgrid = SubgridArraySubgrid(subgrids, is);

            concen_sub = VectorSubvector(concentrations[concen],is);
            concen_dat = SubvectorData(concen_sub);

            nx_v = SubvectorNX(concen_sub);
            ny_v = SubvectorNY(concen_sub);
            nz_v = SubvectorNZ(concen_sub);

            BCConcenPatchExtent(subgrid,&ix,&iy,&iz,&nx,&ny,&nz,ipatch);

            ci = 0;
            BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                      ci, nx_v, ny_v, nz_v, 1, 1, 1,
            {
              if (BoundaryCell(ipatch,i,j,k))
              {
                itmp = SubvectorEltIndex(subvector, i, j, k);
                iv1 = SubvectorEltIndex(concen_sub, i + dir[ipatch][0], j + dir[ipatch][1], k + dir[ipatch][2]);
                iv2 = SubvectorEltIndex(concen_sub, i + 2*dir[ipatch][0], j + 2*dir[ipatch][1], k + 2*dir[ipatch][2]);
                iv3 = SubvectorEltIndex(concen_sub, i + 3*dir[ipatch][0], j + 3*dir[ipatch][1], k + 3*dir[ipatch][2]);

                concen_dat[iv1] = chem_bc_state[(int)tmpp[itmp]].total_mobile.data[concen];
                concen_dat[iv2] = chem_bc_state[(int)tmpp[itmp]].total_mobile.data[concen];
                concen_dat[iv3] = chem_bc_state[(int)tmpp[itmp]].total_mobile.data[concen];
              }
            });
          }
        }
        FreeVector(pfb_indicator);
        break;
      }
    }
  }  
}

#endif


/*--------------------------------------------------------------------------
 * BCConcenCopyPatch
 *   
 *   Copies concentration values from interior cell on boundary
 *   into adjacent boundary cells
 *   3 layers deep 
 *--------------------------------------------------------------------------*/
void BCConcenCopyPatch(Problem *problem, Grid *grid, 
                       Vector **concentrations, 
                       int ipatch)
{
  Subvector      *concen_sub;
  double         *concen_dat;

  int nx_v, ny_v, nz_v;
  int ci;
  int is, i, j, k;
  int iv, iv1, iv2, iv3;
  int num_concen;
  int ix,iy,iz,nx,ny,nz;
  int dir[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };

  Subgrid *subgrid;
  SubgridArray *subgrids = GridSubgrids(grid);
  num_concen = ProblemNumContaminants(problem);

  /*-----------------------------------------------------------------------
   * Implement BC's
   *-----------------------------------------------------------------------*/
  for (int concen = 0; concen < num_concen; concen++)
  {
    ForSubgridI(is, subgrids)
    {

      subgrid = SubgridArraySubgrid(subgrids, is);

      concen_sub = VectorSubvector(concentrations[concen],is);
      concen_dat = SubvectorData(concen_sub);

      nx_v = SubvectorNX(concen_sub);
      ny_v = SubvectorNY(concen_sub);
      nz_v = SubvectorNZ(concen_sub);

      BCConcenPatchExtent(subgrid,&ix,&iy,&iz,&nx,&ny,&nz,ipatch);

      ci = 0;
      BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                ci, nx_v, ny_v, nz_v, 1, 1, 1,
      {
        if (BoundaryCell(ipatch,i,j,k))
        {
          iv = SubvectorEltIndex(concen_sub, i, j, k);
          iv1 = SubvectorEltIndex(concen_sub, i + dir[ipatch][0], j + dir[ipatch][1], k + dir[ipatch][2]);
          iv2 = SubvectorEltIndex(concen_sub, i + 2*dir[ipatch][0], j + 2*dir[ipatch][1], k + 2*dir[ipatch][2]);
          iv3 = SubvectorEltIndex(concen_sub, i + 3*dir[ipatch][0], j + 3*dir[ipatch][1], k + 3*dir[ipatch][2]);

          concen_dat[iv1] = concen_dat[iv];
          concen_dat[iv2] = concen_dat[iv];
          concen_dat[iv3] = concen_dat[iv];
        }
      });
    }
  }
}

/*-----------------------------------------------------------------------
 * BCConcenPatchExtent
 * Determine extent of subgrid patch 
 * - one cell thick in the direction normal to the face
 * - three boundary cell thick in the other two directions
 * - like other ParFlow methods, relies on "left right front back bottom top" ordering of patches in TCL script
 *-----------------------------------------------------------------------*/
void BCConcenPatchExtent(Subgrid *subgrid, int *ix, int *iy, int *iz, int *nx, int *ny, int *nz, int ipatch)
{
      *ix = (ipatch > 1) ? SubgridIX(subgrid)-3 : (ipatch == 0) ? SubgridIX(subgrid) : SubgridIX(subgrid) + SubgridNX(subgrid) - 1;
      *iy = (ipatch < 2 || ipatch > 3) ? SubgridIY(subgrid)-3 : (ipatch == 2) ? SubgridIY(subgrid) : SubgridIY(subgrid) + SubgridNY(subgrid) - 1;
      *iz = (ipatch < 4) ? SubgridIZ(subgrid)-3 : (ipatch == 4) ? SubgridIZ(subgrid) : SubgridIZ(subgrid) + SubgridNZ(subgrid) - 1;

      *nx = (ipatch < 2) ? 1 : SubgridNX(subgrid)+6;
      *ny = (ipatch > 1 && ipatch < 4) ? 1 : SubgridNY(subgrid)+6;
      *nz = (ipatch > 3) ? 1 : SubgridNZ(subgrid)+6;
}

/*--------------------------------------------------------------------------
 * BCConcenCopyAdjacent
 * copy concen of adjacent interior cell into 3 ghost boundary layers
 *
 * Alternative access to BCConcen routines when not built with Alquimia
 *--------------------------------------------------------------------------*/
void BCConcenCopyAdjacent(Problem *problem, Grid *grid, 
                          Vector **concentrations)
{
  int num_domain_patches;
  int domain_index;
  char *switch_name;

  switch_name = GetString("Domain.GeomName");
  domain_index = NA_NameToIndex(GlobalsGeomNames, switch_name);
  
  if (domain_index < 0)
    InputError("Error: invalid geometry name <%s> for key <%s>\n",
                switch_name, "Domain.GeomName");

  num_domain_patches = NA_Sizeof(GlobalsGeometries[domain_index]->patches);

  for (int ipatch = 0; ipatch < num_domain_patches; ipatch++)
  {
    BCConcenCopyPatch(problem, grid, concentrations, 
                      ipatch);
  }
}

/*--------------------------------------------------------------------------
 * BoundaryCell
 * return 1 if subgrid touches current patch, 0 otherwise
 *--------------------------------------------------------------------------*/
int BoundaryCell (int ipatch, int i, int j, int k)
{
  if ((ipatch < 2 && (i == BackgroundX(GlobalsBackground) || i ==
    BackgroundX(GlobalsBackground) + BackgroundNX(GlobalsBackground) - 1))
    ||
    ((ipatch == 2 || ipatch == 3) && (j == BackgroundY(GlobalsBackground) || j ==
      BackgroundY(GlobalsBackground) + BackgroundNY(GlobalsBackground) - 1))
    ||
    ((ipatch == 4 || ipatch == 5) && (k == BackgroundZ(GlobalsBackground) || k ==
      BackgroundZ(GlobalsBackground) + BackgroundNZ(GlobalsBackground) - 1)))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


/*--------------------------------------------------------------------------
 * BCConcentrationInitInstanceXtra
 *--------------------------------------------------------------------------*/
#ifdef HAVE_ALQUIMIA
PFModule *BCConcentrationInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;


#if 0
  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
#endif
  instance_xtra = NULL;

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * BCConcentrationFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void BCConcentrationFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);


  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * BCConcentrationNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *BCConcentrationNewPublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;
  int num_patches = 0;

  Type0          *dummy0;
  Type1          *dummy1;

  int  j;
  char key[IDB_MAX_KEY_LEN];

  char *switch_name;
  int  domain_index;
  char *geochem_cond_names;
  char *patch_name;
  char *bc_patches;

  NameArray type_na;

  type_na = NA_NewNameArray("Constant PFBFile");
  public_xtra = ctalloc(PublicXtra, 1);
  
  geochem_cond_names = GetStringDefault("BCConcentration.GeochemCondition.Names","");
  public_xtra->geochem_cond_na = NA_NewNameArray(geochem_cond_names);
  public_xtra->num_geochemconds = NA_Sizeof(public_xtra->geochem_cond_na); 

  bc_patches = GetStringDefault("BCConcentration.PatchNames","");
  public_xtra->bc_patches_na = NA_NewNameArray(bc_patches);
  public_xtra->num_patches = num_patches = NA_Sizeof(public_xtra->bc_patches_na);
  

    /* Determine the domain geom index from domain name */
    switch_name = GetString("Domain.GeomName");
    domain_index = NA_NameToIndex(GlobalsGeomNames, switch_name);
  
    if (domain_index < 0)
      InputError("Error: invalid geometry name <%s> for key <%s>\n",
                 switch_name, "Domain.GeomName");

    public_xtra->num_domain_patches = NA_Sizeof(GlobalsGeometries[domain_index]->patches);
    
    (public_xtra->input_types) = ctalloc(int, public_xtra->num_domain_patches);
    (public_xtra->data) = ctalloc(void *, public_xtra->num_domain_patches);
  
  
    for (j = 0; j < public_xtra->num_domain_patches; j++)
    {
      patch_name = NA_IndexToName(GlobalsGeometries[domain_index]->patches, j);
  
      sprintf(key, "Patch.%s.BCConcentration.Type", patch_name);
      switch_name = GetStringDefault(key,"");
      public_xtra->input_types[j] = NA_NameToIndex(type_na, switch_name);

      switch ((public_xtra->input_types[j]))
      {
        case 0: //Constant
        {
  
          dummy0 = ctalloc(Type0, 1);
  
          sprintf(key, "Patch.%s.BCConcentration.Value",
                  patch_name);
          dummy0->condition = NA_NameToIndex(public_xtra->geochem_cond_na, GetString(key));
          
          if (dummy0->condition < 0)
          {
            InputError("Error: invalid value <%s> for key <%s>. The geochemical condition must be listed in <BCConcentration.GeochemCondition.Names>\n",
                     GetString(key), key);
          }
  
          (public_xtra->data[j]) = (void*)dummy0;
          break;
        }
  
        case 1: //PFB input of integer corresponding to geochemcondname list
        {  
          dummy1 = ctalloc(Type1, 1);
  
          sprintf(key, "Patch.%s.BCConcentration.FileName",
                  patch_name);
          dummy1->filename = GetString(key);
          (public_xtra->data[j]) = (void*)dummy1;
          break;
        }
  
        case -1: // nothing, defaults to copying adjacent interior cell
        {  
          break;
        }
  
        default: // can't get here
        {
          InputError("Error: invalid type <%s> for key <%s>\n",
                     switch_name, key);
        }
      }
    }

  NA_FreeNameArray(type_na);
  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * BCConcentrationFreePublicXtra
 *--------------------------------------------------------------------------*/

void  BCConcentrationFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Type0  *dummy0;
  Type1  *dummy1;

  int j;

  if (public_xtra)
  {

    NA_FreeNameArray(public_xtra->bc_patches_na);
    NA_FreeNameArray(public_xtra->geochem_cond_na);

      for (j = 0; j < (public_xtra->num_domain_patches); j++)
      {
        switch ((public_xtra->input_types[j]))
        {
          case 0:
            dummy0 = (Type0*)(public_xtra->data[j]);
            tfree(dummy0);
            break;

          case 1:
            dummy1 = (Type1*)(public_xtra->data[j]);
            tfree(dummy1);
            break;

          case -1:
          	break;
        }
      }

    tfree(public_xtra->data);
    tfree(public_xtra->input_types);
    


    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * BCConcentrationSizeOfTempData
 *--------------------------------------------------------------------------*/

int  BCConcentrationSizeOfTempData()
{
  return 0;
}

#endif
