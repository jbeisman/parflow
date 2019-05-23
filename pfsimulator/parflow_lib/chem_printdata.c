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
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"




void PrintChemistryData(ChemPrintFlags *print_flags, AlquimiaSizes *chem_sizes, AlquimiaProblemMetaData *chem_metadata, double t, int file_number, 
  char* file_prefix, Vector **concentrations, Vector **total_immobilePF, Vector **mineral_specific_surfacePF, 
  Vector **surface_site_densityPF, Vector **cation_exchange_capacityPF, Vector *pH, Vector **aqueous_kinetic_ratePF, Vector **mineral_saturation_indexPF, 
  Vector **mineral_reaction_ratePF, Vector **primary_free_ion_concentrationPF, Vector **primary_activity_coeffPF, Vector **secondary_free_ion_concentrationPF,
  Vector **secondary_activity_coeffPF)
{
	char file_type[2048], file_postfix[2048];

	printf("WE ARE INSIDE PRINT FUNCTION\n");

	//primary mobile
	if (print_flags->print_primary_mobile)
	{
		for (int concen = 0; concen < chem_sizes->num_primary; concen++)
		{
			sprintf(file_postfix, "TotalMobile.%02d.%s.%05d", concen, chem_metadata->primary_names.data[concen], file_number);
                  WritePFBinary(file_prefix, file_postfix,
                                 concentrations[concen]);
                  	printf("WE ARE INSIDE PRINT PFB\n");

		}
	}


	if (print_flags->silo_primary_mobile)
	{
		for (int concen = 0; concen < chem_sizes->num_primary; concen++)
		{
			sprintf(file_postfix, "%02d.%s.%05d", concen, chem_metadata->primary_names.data[concen], file_number);
                  sprintf(file_type, "TotalMobile");
                  WriteSilo(file_prefix, file_type, file_postfix, concentrations[concen],
                            t, file_number, "Concentration");
                  	printf("WE ARE INSIDE PRINT SILO\n");

		}
	}



























}




