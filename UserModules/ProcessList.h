/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessList.h
 *
 * Created on 2012-07-04 by Norihiro Watanabe
 */

#pragma once

//------------------------------------------------------------------------------
// List of active modules
// * This file should be called in body files (e.g. cpp)
//------------------------------------------------------------------------------
#include "FemGroundwaterFlow/Head.h"
#include "FemGroundwaterFlow/HeadToElementVelocity.h"
#include "FemGroundwaterFlow/PressureToHead.h"
#include "FemGroundwaterFlow/PressureToElementVelocity.h"
#include "FemGroundwaterFlow/LiquidPressure.h"
#include "FemRichardsFlow/FunctionRichards.h"
#include "FemMassTransport/Concentration.h"
#include "FemHeatTransport/Temperature.h"
#include "FemKinReactGIA/Concentrations.h"
#include "FemReactGIAReduct/ReductConc.h"
#include "FemReactOPS/OPSConc.h"
#include "FemCompMultiPhase/FunctionCMP_2P2C.h"
#include "FemCompMultiPhasePressureForm/FunctionCMP_PressureForm.h"
#include "FemCompMultiPhaseToTalMassDensityForm/FunctionCMP_TotalDensityForm.h"
#include "FemCompMultiPhaseTotalMassEXPForm/FunctionCMP_TotalMassEXPForm.h"
#include "FemCompNonIsoTOTALMassDensityForm/FunctionCMP_NonIso_TotalDensityForm.h"
#include "FemCompNonIsoCapPressureForm/FunctionCMP_NonIso_CapPressureForm.h"
#include "FemCompNonIsoGlobalNCPForm/FunctionCMP_NonIsoGlobalNCPForm.h"
#include "FemCompNonIsoLocalNCPForm/FunctionCMP_NonIso_LocalNCPForm.h"
#include "Fem3Comp2PhaseForm/FunctionCMP_2P3CGlobalNCPForm.h"
#include "FemCompNonIsoTotalMassEXPForm/FunctionCMP_NonIsoTotalMassEXPForm.h"
#include "FemDeformationTotalForm/Displacement.h"
#include "FemDeformationTotalForm/ElementStressStrain.h"
#include "FemDeformationTotalForm/NodalStressStrain.h"
#include "FemPoroelastic/DisplacementPressure.h"
#include "Xfem/XFEM_EXAMPLE_CRACK1.h"
#include "THMmf/Pmf.h"
#include "THMmf/Tmf.h"
#include "THMmf/Umf.h"
#include "FemDeformationIncrementalForm/IncrementalDisplacement.h"
//#include "ProcessBuilder.h"

//OGS_LINK_PROCESS(GROUNDWATER_FLOW, FunctionHead);


