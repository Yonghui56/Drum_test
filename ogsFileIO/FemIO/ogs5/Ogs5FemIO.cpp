/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FemIO.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Ogs5FemIO.h"

#include "logog.hpp"

#include "GeoIO/Ogs4GeoIO.h"
#include "MeshIO/MeshIoOgs5.h"


#include "BaseLib/CodingTools.h"

namespace ogs5
{

bool Ogs5FemIO::read(const std::string &proj_path, Ogs5FemData &ogs5data)
{
    PCSRead(proj_path, ogs5data.pcs_vector);
    MFPRead(proj_path, ogs5data.mfp_vector);
    MSPRead(proj_path, ogs5data.msp_vector);
    MMPRead(proj_path, ogs5data.mmp_vector);
    CPRead( proj_path, ogs5data.cp_vector);
    BCRead( proj_path, ogs5data.bc_vector);
    STRead( proj_path, ogs5data.st_vector);
    ICRead( proj_path, ogs5data.ic_vector);
    OUTRead(proj_path, ogs5data.out_vector);
    TIMRead(proj_path, ogs5data.time_vector);
    NUMRead(proj_path, ogs5data.num_vector);
	KRRead( proj_path, ogs5data.KinReact_vector, 
		               ogs5data.KinReactData_vector, 
					   ogs5data.KinBlob_vector); 
	CURRead(proj_path, ogs5data.kurven_vector);

    // set primary variable name
    size_t mass_transport_count = 0;
    for (size_t i=0; i<ogs5data.pcs_vector.size(); i++) {
        CRFProcess* pcs = ogs5data.pcs_vector[i];
        switch (pcs->getProcessType()) {
        case FiniteElement::GROUNDWATER_FLOW:
			pcs->primary_variable_name.push_back("HEAD");
            break;
        case FiniteElement::LIQUID_FLOW:
            pcs->primary_variable_name.push_back("PRESSURE1");
            break;
		case FiniteElement::RICHARDS_FLOW:
            pcs->primary_variable_name.push_back("PRESSURE1"); 
            break; 
        case FiniteElement::HEAT_TRANSPORT:
            pcs->primary_variable_name.push_back("TEMPERATURE1");
            break;
        case FiniteElement::MASS_TRANSPORT:
			pcs->primary_variable_name.push_back(ogs5data.cp_vector[mass_transport_count]->compname);
            mass_transport_count++;
            break;
        case FiniteElement::DEFORMATION:
            pcs->primary_variable_name.push_back("DISPLACEMENT");
            break;
        case FiniteElement::DEFORMATION_FLOW:
            pcs->primary_variable_name.push_back("DISPLACEMENT");
            pcs->primary_variable_name.push_back("PRESSURE1");
            break;
		case FiniteElement::KIN_REACT_GIA:
			for (size_t i=0; i<ogs5data.cp_vector.size() ; i++)
				pcs->primary_variable_name.push_back( ogs5data.cp_vector[i]->compname); 
			break;
        case FiniteElement::REACT_GIA:
			for (size_t i=0; i<ogs5data.cp_vector.size() ; i++)
				pcs->primary_variable_name.push_back( ogs5data.cp_vector[i]->compname);
			break;
        case FiniteElement::REACT_TRANS_OPS:
			for (size_t i=0; i<ogs5data.cp_vector.size() ; i++)
				pcs->primary_variable_name.push_back( ogs5data.cp_vector[i]->compname);
			break; 
		case FiniteElement::CMP_2P2C:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("MOLAR_FRACTION");
			break;
		case FiniteElement::CMP_PressureForm:
			pcs->primary_variable_name.push_back("LIQUID_PRESSURE");
			pcs->primary_variable_name.push_back("MOLAR_FRACTION");
			break;
		case FiniteElement::CMP_TotalDensityForm:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("TOTAL_MASS_DENSITY");
			break;
		case FiniteElement::CMP_NonIso_TotalDensityForm:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("TOTAL_MASS_DENSITY");
			pcs->primary_variable_name.push_back("TEMPERATURE");
			break;
		case FiniteElement::CMP_NonIso_CapPressureForm:
			pcs->primary_variable_name.push_back("GAS_PRESSURE");
			pcs->primary_variable_name.push_back("CAPILLARY_PRESSURE");
			pcs->primary_variable_name.push_back("TEMPERATURE");
			break;
		case FiniteElement::CMP_NonIso_GlobalNCPForm:
			pcs->primary_variable_name.push_back("GAS_PRESSURE");
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_AIR");
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_VAPOR");
			pcs->primary_variable_name.push_back("SATURATION");
			pcs->primary_variable_name.push_back("TEMPERATURE");
			break;
		case FiniteElement::CMP_NonIso_LocalNCPForm:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("TOTAL_MASS_DENSITY");
			pcs->primary_variable_name.push_back("TEMPERATURE");
			break;
		case FiniteElement::CMP_2P3CGlobalNCPForm:
			pcs->primary_variable_name.push_back("GAS_PRESSURE");
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_GAS_H");//molar fraction of hydrogen in gas phase
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_GAS_C");//molar fraction of co2 in gas phase
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_GAS_W");// molar fraction of water vapor in gas phase
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_GAS_AIR");// molar fraction of water vapor in gas phase
			pcs->primary_variable_name.push_back("MOLAR_FRACTION_GAS_CO2");// molar fraction of water vapor in gas phase
			pcs->primary_variable_name.push_back("SATURATION");
			break;
		case FiniteElement::CMP_TotalMassEXPForm:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("TOTAL_MASS_DENSITY");//molar fraction of hydrogen in gas phase
			break;
		case FiniteElement::CMP_NonIsoTotalMassEXPForm:
			pcs->primary_variable_name.push_back("MEAN_PRESSURE");
			pcs->primary_variable_name.push_back("TOTAL_MASS_DENSITY");//molar fraction of hydrogen in gas phase
			pcs->primary_variable_name.push_back("TEMPERATURE");//molar fraction of hydrogen in gas phase
			break;
		default:
		    break;
        }
    }

    // gli
    std::vector<std::string> geo_errors;
    ogs5data.geo_obj = new GeoLib::GEOObjects();
    Ogs4GeoIO::readGLIFileV4(proj_path+".gli", ogs5data.geo_obj, ogs5data.geo_unique_name, geo_errors);
    if (geo_errors.size()>0) {
        ERR(" Error when reading GLI");
        for (size_t i=0; i<geo_errors.size(); i++)
            ERR(" -> %s", geo_errors[i].c_str());
//        return false;
    }

    // mesh
    MeshIoOgs5::readMesh(proj_path+".msh", ogs5data.list_mesh);
    if (ogs5data.list_mesh.size()==0) {
        ERR(" Error: Cannot find mesh - %s.msh", proj_path.c_str());
//        return false;
    }

    return true;


//    RCRead(proj_path);
//    KRRead(proj_path, geo_obj, unique_name);
//    KRWrite(proj_path);
//
//
//    PCTRead(proj_path);                   // PCH
//    FMRead(proj_path);                    // PCH
//    FCTRead(proj_path);                   //OK
//    CURRead(proj_path);                   //OK
}

}
