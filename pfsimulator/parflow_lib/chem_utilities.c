

void Chem2PF_Single(Vector *pf_vector, double *chem_var)
{
  Grid       *grid = VectorGrid(pf_vector);
  int sg;

  ForSubgridI(sg, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, sg);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);

    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    int chem_index, pf_index;

    Subvector *subvector = VectorSubvector(pf_vector, sg);
    double *subvector_data = SubvectorData(subvector);

    int i, j, k;

    for (i = ix; i < ix + nx; i++)
    {
      for (j = iy; j < iy + ny; j++)
      {
        for (k = iz; k < iz + wrf_depth; k++)
        {
          pf_index   = SubvectorEltIndex(subvector, i, j, k);
          chem_index = (i-ix ) + (j- iy) * (nx) + (k- iz) * (nx) * (ny);

          subvector_data[pf_index] = chem_var[chem_index];
        }
      }
    }
  }
}






void PF2Chem_Single(Vector *pf_vector, double *chem_var)
{
  Grid       *grid = VectorGrid(pf_vector);
  int sg;

  ForSubgridI(sg, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, sg);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);

    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    int chem_index, pf_index;

    Subvector *subvector = VectorSubvector(pf_vector, sg);
    double *subvector_data = SubvectorData(subvector);

    int i, j, k;

    for (i = ix; i < ix + nx; i++)
    {
      for (j = iy; j < iy + ny; j++)
      {
        for (k = iz; k < iz + wrf_depth; k++)
        {
          pf_index   = SubvectorEltIndex(subvector, i, j, k);
          chem_index = (i-ix ) 
                       + (j- iy) * (nx) 
                       + (k- iz) * (nx) * (ny);

          chem_var[chem_index] = subvector_data[pf_index];
        }
      }
    }
  }
}






void Chem2PF_Multi(Vector *pf_vector, double *chem_var, int num_var)
{
  Grid       *grid = VectorGrid(pf_vector);
  int sg, var;

  for(var = 0; var < num_var; var++)
  {
    ForSubgridI(sg, GridSubgrids(grid))
    {
      Subgrid *subgrid = GridSubgrid(grid, sg);

      int ix = SubgridIX(subgrid);
      int iy = SubgridIY(subgrid);
      int iz = SubgridIZ(subgrid);

      int nx = SubgridNX(subgrid);
      int ny = SubgridNY(subgrid);
      int nz = SubgridNZ(subgrid);

      int chem_index, pf_index;

      Subvector *subvector = VectorSubvector(pf_vector, sg);
      double *subvector_data = SubvectorData(subvector);

      int i, j, k;

      for (i = ix; i < ix + nx; i++)
      {
        for (j = iy; j < iy + ny; j++)
        {
          for (k = iz; k < iz + wrf_depth; k++)
          {
            pf_index   = SubvectorEltIndex(subvector, i, j, k);
            chem_index = var + (i-ix ) * num_var 
                         + (j- iy) * (nx) * num_var 
                         + (k- iz) * (nx) * (ny) * num_var;

            subvector_data[pf_index] = chem_var[chem_index];
          }
        }
      }
    }
  }
}




void PF2Chem_Multi(Vector *pf_vector, double *chem_var, int num_var)
{
  Grid       *grid = VectorGrid(pf_vector);
  int sg, var;

  for(var = 0; var < num_var; var++)
  {
    ForSubgridI(sg, GridSubgrids(grid))
    {
      Subgrid *subgrid = GridSubgrid(grid, sg);

      int ix = SubgridIX(subgrid);
      int iy = SubgridIY(subgrid);
      int iz = SubgridIZ(subgrid);

      int nx = SubgridNX(subgrid);
      int ny = SubgridNY(subgrid);
      int nz = SubgridNZ(subgrid);

      int chem_index, pf_index;

      Subvector *subvector = VectorSubvector(pf_vector, sg);
      double *subvector_data = SubvectorData(subvector);

      int i, j, k;

      for (i = ix; i < ix + nx; i++)
      {
        for (j = iy; j < iy + ny; j++)
        {
          for (k = iz; k < iz + wrf_depth; k++)
          {
            pf_index   = SubvectorEltIndex(subvector, i, j, k);
            chem_index = var + (i-ix ) * num_var 
                         + (j- iy) * (nx) * num_var 
                         + (k- iz) * (nx) * (ny) * num_var;

            chem_var[chem_index] = subvector_data[pf_index];
          }
        }
      }
    }
  }
}

