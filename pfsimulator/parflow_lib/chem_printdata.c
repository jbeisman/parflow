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
* Print all of the chemistry data the user asked for
*
*-----------------------------------------------------------------------------
*
*****************************************************************************/
#include "parflow.h"
#include "pf_alquimia.h"


#ifdef HAVE_ALQUIMIA
void PrintChemistryData(ChemPrintFlags *print_flags, AlquimiaSizes *chem_sizes, AlquimiaProblemMetaData *chem_metadata, double t, int file_number, 
  char* file_prefix, int *any_file_dumped, Vector **concentrations, Vector **total_immobilePF, Vector **mineral_specific_surfacePF, Vector **mineral_volume_fractionsPF,
  Vector **surface_site_densityPF, Vector **cation_exchange_capacityPF, Vector *pH, Vector **aqueous_kinetic_ratePF, Vector **mineral_saturation_indexPF, 
  Vector **mineral_reaction_ratePF, Vector **primary_free_ion_concentrationPF, Vector **primary_activity_coeffPF, Vector **secondary_free_ion_concentrationPF,
  Vector **secondary_activity_coeffPF)
{
	char file_type[2048], file_postfix[2048];

	//primary mobile
	if (print_flags->print_primary_mobile)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "PrimaryMobile.%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 concentrations[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_primary_mobile)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  sprintf(file_type, "PrimaryMobile");
                  WriteSilo(file_prefix, file_type, file_postfix, concentrations[i],
                            t, file_number, "Concen");
		}
		*any_file_dumped = 1;
	}



	//primary sorbed
	if (print_flags->print_sorbed)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "PrimarySorbed.%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 total_immobilePF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_sorbed)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  sprintf(file_type, "PrimarySorbed");
                  WriteSilo(file_prefix, file_type, file_postfix, total_immobilePF[i],
                            t, file_number, "Sorbed");
		}
		*any_file_dumped = 1;
	}



	//mineral volfx
	if (print_flags->print_mineral_volfx)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "MineralVolfx.%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 mineral_volume_fractionsPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_mineral_volfx)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  sprintf(file_type, "MineralVolfx");
                  WriteSilo(file_prefix, file_type, file_postfix, mineral_volume_fractionsPF[i],
                            t, file_number, "Volfx");
		}
		*any_file_dumped = 1;
	}



	//mineral surface area
	if (print_flags->print_mineral_surfarea)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "MineralSurfArea.%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 mineral_specific_surfacePF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_mineral_surfarea)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  sprintf(file_type, "MineralSurfArea");
                  WriteSilo(file_prefix, file_type, file_postfix, mineral_specific_surfacePF[i],
                            t, file_number, "SurfArea");
		}
		*any_file_dumped = 1;
	}



	//mineral saturation index
	if (print_flags->print_mineral_SI)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "MineralSI.%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 mineral_saturation_indexPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_mineral_SI)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  sprintf(file_type, "MineralSI");
                  WriteSilo(file_prefix, file_type, file_postfix, mineral_saturation_indexPF[i],
                            t, file_number, "SatIndex");
		}
		*any_file_dumped = 1;
	}



	//mineral kinetic rate
	if (print_flags->print_mineral_rate)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "MineralRate.%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 mineral_reaction_ratePF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_mineral_rate)
	{
		for (int i = 0; i < chem_sizes->num_minerals; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->mineral_names.data[i], file_number);
                  sprintf(file_type, "MineralRate");
                  WriteSilo(file_prefix, file_type, file_postfix, mineral_reaction_ratePF[i],
                            t, file_number, "MineralRate");
		}
		*any_file_dumped = 1;
	}



	//surface site density
	if (print_flags->print_surf_dens)
	{
		for (int i = 0; i < chem_sizes->num_surface_sites; i++)
		{
			sprintf(file_postfix, "SurfSiteDens.%02d.%s.%05d", i, chem_metadata->surface_site_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 surface_site_densityPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_surf_dens)
	{
		for (int i = 0; i < chem_sizes->num_surface_sites; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->surface_site_names.data[i], file_number);
                  sprintf(file_type, "SurfSiteDens");
                  WriteSilo(file_prefix, file_type, file_postfix, surface_site_densityPF[i],
                            t, file_number, "SurfSiteDens");
		}
		*any_file_dumped = 1;
	}



	//cation exchange capacity
	if (print_flags->print_CEC)
	{
		for (int i = 0; i < chem_sizes->num_ion_exchange_sites; i++)
		{
			sprintf(file_postfix, "CEC.%02d.%s.%05d", i, chem_metadata->ion_exchange_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 cation_exchange_capacityPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_CEC)
	{
		for (int i = 0; i < chem_sizes->num_ion_exchange_sites; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->ion_exchange_names.data[i], file_number);
                  sprintf(file_type, "CEC");
                  WriteSilo(file_prefix, file_type, file_postfix, cation_exchange_capacityPF[i],
                            t, file_number, "CEC");
		}
		*any_file_dumped = 1;
	}



	//pH
	if (print_flags->print_pH)
	{
			sprintf(file_postfix, "pH.%05d", file_number);
        		  WritePFBinary(file_prefix, file_postfix,
                                pH);
		*any_file_dumped = 1;
	}

	if (print_flags->silo_pH)
	{
			sprintf(file_postfix, "%05d", file_number);
                  sprintf(file_type, "pH");
                  WriteSilo(file_prefix, file_type, file_postfix, pH,
                            t, file_number, "pH");
		*any_file_dumped = 1;
	}



	//aqueous kinetic rate
	if (print_flags->print_aqueous_rate)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_kinetics; i++)
		{
			sprintf(file_postfix, "AqueousRate.%02d.%s.%05d", i, chem_metadata->aqueous_kinetic_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 aqueous_kinetic_ratePF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_aqueous_rate)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_kinetics; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->aqueous_kinetic_names.data[i], file_number);
                  sprintf(file_type, "AqueousRate");
                  WriteSilo(file_prefix, file_type, file_postfix, aqueous_kinetic_ratePF[i],
                            t, file_number, "AqueousRate");
		}
		*any_file_dumped = 1;
	}



	//primary free ion concentration
	if (print_flags->print_primary_freeion)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "PrimaryFreeIon.%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 primary_free_ion_concentrationPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_primary_freeion)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  sprintf(file_type, "PrimaryFreeIon");
                  WriteSilo(file_prefix, file_type, file_postfix, primary_free_ion_concentrationPF[i],
                            t, file_number, "PrimaryFreeIon");
		}
		*any_file_dumped = 1;
	}



	//primary activity
	if (print_flags->print_primary_activity)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "PrimaryActivity.%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 primary_activity_coeffPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_primary_activity)
	{
		for (int i = 0; i < chem_sizes->num_primary; i++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", i, chem_metadata->primary_names.data[i], file_number);
                  sprintf(file_type, "PrimaryActivity");
                  WriteSilo(file_prefix, file_type, file_postfix, primary_activity_coeffPF[i],
                            t, file_number, "PrimaryActivity");
		}
		*any_file_dumped = 1;
	}



	//secondary free ion concentration
	if (print_flags->print_secondary_freeion)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_complexes; i++)
		{
			sprintf(file_postfix, "SecondaryFreeIon.%02d.%05d", i, file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 secondary_free_ion_concentrationPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_secondary_freeion)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_complexes; i++)
		{
			sprintf(file_postfix, "%02d.%05d", i, file_number);
                  sprintf(file_type, "SecondaryFreeIon");
                  WriteSilo(file_prefix, file_type, file_postfix, secondary_free_ion_concentrationPF[i],
                            t, file_number, "SecondaryFreeIon");
		}
		*any_file_dumped = 1;
	}



	//secondary activity
	if (print_flags->print_secondary_activity)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_complexes; i++)
		{
			sprintf(file_postfix, "SecondaryActivity.%02d.%05d", i, file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 secondary_activity_coeffPF[i]);
		}
		*any_file_dumped = 1;
	}

	if (print_flags->silo_secondary_activity)
	{
		for (int i = 0; i < chem_sizes->num_aqueous_complexes; i++)
		{
			sprintf(file_postfix, "%02d.%05d", i, file_number);
                  sprintf(file_type, "SecondaryActivity");
                  WriteSilo(file_prefix, file_type, file_postfix, secondary_activity_coeffPF[i],
                            t, file_number, "SecondaryActivity");
		}
		*any_file_dumped = 1;
	}
}
#endif

