#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
pfset FileVersion 4

#-----------------------------------------------------------------------------
# Process Topology
#-----------------------------------------------------------------------------

pfset Process.Topology.P        1
pfset Process.Topology.Q        1
pfset Process.Topology.R        1

#-----------------------------------------------------------------------------
# Computational Grid
#-----------------------------------------------------------------------------
pfset ComputationalGrid.Lower.X                0.0
pfset ComputationalGrid.Lower.Y                 0.0
pfset ComputationalGrid.Lower.Z                  0.0

pfset ComputationalGrid.DX	                 2.0
pfset ComputationalGrid.DY                   2.0
pfset ComputationalGrid.DZ	                 3.0

pfset ComputationalGrid.NX                      50
pfset ComputationalGrid.NY                      50
pfset ComputationalGrid.NZ                      1

#-----------------------------------------------------------------------------
# The Names of the GeomInputs
#-----------------------------------------------------------------------------
pfset GeomInput.Names "domain_input background_input source_region_input concen_region_input"

#-----------------------------------------------------------------------------
# Domain Geometry Input
#-----------------------------------------------------------------------------
pfset GeomInput.domain_input.InputType            Box
pfset GeomInput.domain_input.GeomName             domain

#-----------------------------------------------------------------------------
# Domain Geometry
#-----------------------------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0 
pfset Geom.domain.Lower.Y                         0.0
pfset Geom.domain.Lower.Z                          0.0

pfset Geom.domain.Upper.X                        100.0
pfset Geom.domain.Upper.Y                        100.0
pfset Geom.domain.Upper.Z                          3.0

pfset Geom.domain.Patches "left right front back bottom top"

#-----------------------------------------------------------------------------
# Background Geometry Input
#-----------------------------------------------------------------------------
pfset GeomInput.background_input.InputType         Box
pfset GeomInput.background_input.GeomName          background

#-----------------------------------------------------------------------------
# Background Geometry
#-----------------------------------------------------------------------------
pfset Geom.background.Lower.X -99999999.0
pfset Geom.background.Lower.Y -99999999.0
pfset Geom.background.Lower.Z -99999999.0

pfset Geom.background.Upper.X  99999999.0
pfset Geom.background.Upper.Y  99999999.0
pfset Geom.background.Upper.Z  99999999.0

#-----------------------------------------------------------------------------
# Source_Region Geometry Input
#-----------------------------------------------------------------------------
pfset GeomInput.source_region_input.InputType      Box
pfset GeomInput.source_region_input.GeomName       source_region

#-----------------------------------------------------------------------------
# Source_Region Geometry
#-----------------------------------------------------------------------------
pfset Geom.source_region.Lower.X    0.0
pfset Geom.source_region.Lower.Y    0.0
pfset Geom.source_region.Lower.Z    0.0

pfset Geom.source_region.Upper.X    100.0
pfset Geom.source_region.Upper.Y    100.0
pfset Geom.source_region.Upper.Z     3.0



#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names "background"

pfset Geom.background.Perm.Type     Constant
pfset Geom.background.Perm.Value    4.0

pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "background"

pfset Geom.background.Perm.TensorValX  1.0
pfset Geom.background.Perm.TensorValY  1.0
pfset Geom.background.Perm.TensorValZ  1.0


#-----------------------------------------------------------------------------
# Concen_Region Geometry Input
#-----------------------------------------------------------------------------
pfset GeomInput.concen_region_input.InputType       Box
pfset GeomInput.concen_region_input.GeomName        concen_region

#-----------------------------------------------------------------------------
# Concen_Region Geometry
#-----------------------------------------------------------------------------
pfset Geom.concen_region.Lower.X  10.0
pfset Geom.concen_region.Lower.Y  10.0
pfset Geom.concen_region.Lower.Z  0.0

pfset Geom.concen_region.Upper.X  30.0
pfset Geom.concen_region.Upper.Y  30.0
pfset Geom.concen_region.Upper.Z   3.0

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	Constant
pfset Phase.water.Density.Value	1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names			"tce dummy dummy2"
pfset Contaminants.tce.Degradation.Value	 0.0

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

pfset TimingInfo.BaseUnit		0.5
pfset TimingInfo.StartCount		0
pfset TimingInfo.StartTime		0.0
pfset TimingInfo.StopTime            60.0
pfset TimingInfo.DumpInterval	     1

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          background

pfset Geom.background.Porosity.Type    Constant
pfset Geom.background.Porosity.Value   0.25

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Mobility
#-----------------------------------------------------------------------------
pfset Phase.water.Mobility.Type        Constant
pfset Phase.water.Mobility.Value       1.0

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames           background
pfset Geom.background.tce.Retardation.Type     Linear
pfset Geom.background.tce.Retardation.Rate     0.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names constant
pfset Cycle.constant.Names		"alltime"
pfset Cycle.constant.alltime.Length	 20.0
pfset Cycle.constant.Repeat		-1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames "left right front back bottom top"

#pfset Patch.top.BCPressure.Type			DirEquilRefPatch
#pfset Patch.top.BCPressure.Cycle			"constant"
#pfset Patch.top.BCPressure.RefGeom			domain
#pfset Patch.top.BCPressure.RefPatch			bottom
#pfset Patch.top.BCPressure.alltime.Value		3.0
#
#pfset Patch.bottom.BCPressure.Type			DirEquilRefPatch
#pfset Patch.bottom.BCPressure.Cycle			"constant"
#pfset Patch.bottom.BCPressure.RefGeom			domain
#pfset Patch.bottom.BCPressure.RefPatch			bottom
#pfset Patch.bottom.BCPressure.alltime.Value		0.0
pfset Patch.front.BCPressure.Type			DirEquilRefPatch
pfset Patch.front.BCPressure.Cycle			"constant"
pfset Patch.front.BCPressure.RefGeom			domain
pfset Patch.front.BCPressure.RefPatch			bottom
pfset Patch.front.BCPressure.alltime.Value		15.0

pfset Patch.back.BCPressure.Type			DirEquilRefPatch
pfset Patch.back.BCPressure.Cycle			"constant"
pfset Patch.back.BCPressure.RefGeom			domain
pfset Patch.back.BCPressure.RefPatch			bottom
pfset Patch.back.BCPressure.alltime.Value		10.0

pfset Patch.left.BCPressure.Type			DirEquilRefPatch
pfset Patch.left.BCPressure.Cycle			"constant"
pfset Patch.left.BCPressure.RefGeom			domain
pfset Patch.left.BCPressure.RefPatch			bottom
pfset Patch.left.BCPressure.alltime.Value		15.0

pfset Patch.right.BCPressure.Type			DirEquilRefPatch
pfset Patch.right.BCPressure.Cycle			"constant"
pfset Patch.right.BCPressure.RefGeom			domain
pfset Patch.right.BCPressure.RefPatch			bottom
pfset Patch.right.BCPressure.alltime.Value		10.0

pfset Patch.top.BCPressure.Type				FluxConst
pfset Patch.top.BCPressure.Cycle			"constant"
pfset Patch.top.BCPressure.alltime.Value		0.0

pfset Patch.bottom.BCPressure.Type			FluxConst
pfset Patch.bottom.BCPressure.Cycle			"constant"
pfset Patch.bottom.BCPressure.alltime.Value		0.0

#pfset Patch.left.BCPressure.Type			FluxConst
#pfset Patch.left.BCPressure.Cycle			"constant"
#pfset Patch.left.BCPressure.alltime.Value		-0.25
#
#pfset Patch.right.BCPressure.Type			 FluxConst
#pfset Patch.right.BCPressure.Cycle			"constant"
#pfset Patch.right.BCPressure.alltime.Value		-0.25
#
#pfset Patch.back.BCPressure.Type			FluxConst
#pfset Patch.back.BCPressure.Cycle			"constant"
#pfset Patch.back.BCPressure.alltime.Value		0.0
#
#pfset Patch.front.BCPressure.Type			FluxConst
#pfset Patch.front.BCPressure.Cycle			"constant"
#pfset Patch.front.BCPressure.alltime.Value		-0.25






#-----------------------------------------------------------------------------
# Boundary Conditions: Temperature 
#-----------------------------------------------------------------------------
pfset BCTemperature.PatchNames        "left right front back bottom top"
 
pfset Patch.left.BCTemperature.Type                      DirConst
pfset Patch.left.BCTemperature.Cycle                     "constant"
pfset Patch.left.BCTemperature.alltime.Value             288.15
 
pfset Patch.right.BCTemperature.Type                     DirConst
pfset Patch.right.BCTemperature.Cycle                    "constant"
pfset Patch.right.BCTemperature.alltime.Value            295.
pfset Patch.right.BCTemperature.alltime.Value            293.15
 
pfset Patch.front.BCTemperature.Type                     DirConst 
pfset Patch.front.BCTemperature.Type                     FluxConst 
pfset Patch.front.BCTemperature.Cycle                    "constant"
pfset Patch.front.BCTemperature.alltime.Value            305.
pfset Patch.front.BCTemperature.alltime.Value            0.0 
 
pfset Patch.back.BCTemperature.Type                      DirConst 
pfset Patch.back.BCTemperature.Type                      FluxConst 
pfset Patch.back.BCTemperature.Cycle                     "constant"
pfset Patch.back.BCTemperature.alltime.Value             295.
pfset Patch.back.BCTemperature.alltime.Value             0.0
 
pfset Patch.bottom.BCTemperature.Type                    FluxConst 
pfset Patch.bottom.BCTemperature.Cycle                   "constant"
pfset Patch.bottom.BCTemperature.alltime.Value           0.0
 
pfset Patch.top.BCTemperature.Type                       FluxConst 
pfset Patch.top.BCTemperature.Cycle                      "constant"
pfset Patch.top.BCTemperature.alltime.Value              0.0

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    background
pfset PhaseSources.water.Geom.background.Value        0.0

pfset PhaseSources.Type                         Constant
pfset PhaseSources.GeomNames                    background
pfset PhaseSources.Geom.background.Value               0.0

pfset PhaseConcen.water.tce.Type                      Constant
pfset PhaseConcen.water.tce.GeomNames                 "concen_region"
pfset PhaseConcen.water.tce.Geom.concen_region.Value  0.1


#-----------------------------------------------------------------------------
# Temperature sources:
#-----------------------------------------------------------------------------
pfset TempSources.Type                         Constant
pfset TempSources.GeomNames                   "background"
pfset TempSources.Geom.background.Value        0.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "background"
pfset Geom.background.SpecificStorage.Value          1.0e-5

#-----------------------------------------------------------------------------
# Heat Capacity 
#-----------------------------------------------------------------------------

pfset Phase.water.HeatCapacity.Type                      Constant
pfset Phase.water.HeatCapacity.GeomNames                 "background"
pfset Phase.water.Geom.background.HeatCapacity.Value        4000. 

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.0

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.0

#---------------------------------------------------------
# Mannings coefficient 
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 2.3e-7



#---------------------------------------------------------
# ALQUIMIA INPUT VARS
#---------------------------------------------------------
pfset Solver.Chemistry True
pfset Chemistry.Engine CrunchFlow
pfset Chemistry.InputFile calcite-1d-crunch.in


# order of geomnames matters
# just like everything else w/ PF geometries,
# geominputs listed later will overwrite those
# listed earlier if they overlap
pfset GeochemCondition.Type "Constant"
pfset GeochemCondition.GeomNames "source_region concen_region"
pfset GeochemCondition.Names "initial west"
pfset GeochemCondition.Geom.source_region.Value "initial"
pfset GeochemCondition.Geom.concen_region.Value "west"

pfset BCConcentration.GeochemCondition.Names "west"
pfset BCConcentration.PatchNames "left front"
pfset Patch.left.BCConcentration.Type Constant
pfset Patch.left.BCConcentration.Value west

pfset Patch.front.BCConcentration.Type Constant
pfset Patch.front.BCConcentration.Value west

pfset Chemistry.PrintMineralRate True
pfset Chemistry.ParFlowTimeUnits minutes


pfset Solver.WriteSiloConcentration True
pfset Solver.WriteSiloPressure True

#-----------------------------------------------------------------------------
# The Solver Impes MaxIter default value changed so to get previous
# results we need to set it back to what it was
#-----------------------------------------------------------------------------
pfset Solver.MaxIter 5000
pfset Solver.CFL 0.9
pfset Solver.AdvectOrder 2
pfset Solver.RelTol 1.0e-35
pfset Solver.AbsTol 1.0e-50
pfset Solver.Nonlinear.ResidualTol 1.0e-15

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfrun alquimia_driver
pfundist alquimia_driver

#-----------------------------------------------------------------------------
# If running as test; check output.
# You do not need this for normal PF input files; this is done so the examples
# are run and checked as part of our testing process.
#-----------------------------------------------------------------------------

