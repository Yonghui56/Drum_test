/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2015-03-19 by Yonghui HUANG
*/

#ifndef NON_LINEAR_CMP_NONISO_LOCALNCPFORM_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_NONISO_LOCALNCPFORM_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "Ogs6FemData.h"
//#include "EOS_NonIso_TotalDensityForm.h"
#include "EOS_NonIso_HeatPipeProb.h"

#include "LocalProblem_EOS_NonIso_LocalNCP.h"


/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_NonIso_LocalNCPForm_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_NonIso_LocalNCPForm_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
	{
		_EOS = new EOS_NonIso_HeatPipeProb();
		_LP_EOS = new LocalProblem_EOS_NonIso_LocalNCP();
	};
	
	virtual ~NonLinearCMP_NonIso_LocalNCPForm_TimeODELocalAssembler() {
		BaseLib::releaseObject(_EOS);
		BaseLib::releaseObject(_LP_EOS);
		M.resize(0, 0);
		D.resize(0, 0);
		
	};

protected:
	
	virtual void assembleODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
	{
		// -------------------------------------------------
		// current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		// integration method 
		FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
		// initialize the local EOS problem
		// now get the sizes
		const std::size_t n_dim = e.getDimension();
		const std::size_t mat_id = e.getGroupID();
		const std::size_t ele_id = e.getID();
		const std::size_t n_dof = u1.size();
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
		const std::size_t n_gsp = q->getNumberOfSamplingPoints(); // number of Gauss points
		std::size_t node_id(0);  // index of the node

		// get the instance of porous media class
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		// get the fluid property for gas phase
		MaterialLib::Fluid* fluid1 = Ogs6FemData::getInstance()->list_fluid[0];
		// get the fluid property for liquid phase
		MaterialLib::Fluid* fluid2 = Ogs6FemData::getInstance()->list_fluid[1];
		// get the first (light) component - hydrogen
		MaterialLib::Compound* component1 = Ogs6FemData::getInstance()->list_compound[0];
		// get the second (heavy) component - water
		MaterialLib::Compound* component2 = Ogs6FemData::getInstance()->list_compound[1];
		MaterialLib::Solid* solid = Ogs6FemData::getInstance()->list_solid[0];
		const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
		double geo_area = 1.0;
		pm->geo_area->eval(e_pos, geo_area);
		const bool hasGravityEffect =  _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		//initialize

		P_gp = LocalMatrixType::Zero(1,1);
		X_gp = LocalMatrixType::Zero(1, 1);
		T_gp = LocalMatrixType::Zero(1, 1);

		Input = LocalVectorType::Zero(3);
		Output = LocalVectorType::Zero(3);

		S_gp = LocalVectorType::Zero(1);
		rho_G_h_gp = LocalVectorType::Zero(1);
		rho_G_w_gp =  LocalVectorType::Zero(1);
		rho_L_h_gp = LocalVectorType::Zero(1);
		P_Gw_gp = LocalVectorType::Zero(1);
		P_Gh_gp = LocalVectorType::Zero(1);
		PC_gp = LocalVectorType::Zero(1);

		massfractionXair_gp = LocalVectorType::Zero(1);
		Xvap_gp = LocalVectorType::Zero(1);

		drho_G_hdP_gp = LocalVectorType::Zero(1);
		drho_G_hdX_gp = LocalVectorType::Zero(1);
		drho_G_h_dT_gp =  LocalVectorType::Zero(1);
		drho_G_w_dP_gp =  LocalVectorType::Zero(1);
		drho_G_w_dX_gp = LocalVectorType::Zero(1);
		drho_G_w_dT_gp = LocalVectorType::Zero(1);
		drho_L_h_dP_gp = LocalVectorType::Zero(1);
		drho_L_h_dX_gp = LocalVectorType::Zero(1);
		drho_L_h_dT_gp = LocalVectorType::Zero(1);
		dSgdP_gp = LocalVectorType::Zero(1);
		dSgdX_gp = LocalVectorType::Zero(1);
		dSgdT_gp = LocalVectorType::Zero(1);
		dPGdSg_gp = LocalVectorType::Zero(1);
		dPGw_dP_gp = LocalVectorType::Zero(1);
		dPGw_dT_gp = LocalVectorType::Zero(1);
		dPGh_dP_gp = LocalVectorType::Zero(1);
		dPGh_dT_gp = LocalVectorType::Zero(1);
		dPC_dSg_gp = LocalVectorType::Zero(1);
		EnthalpyH_G_gp = LocalVectorType::Zero(1);
		EnthalpyH_L_gp = LocalVectorType::Zero(1);
		dmassXair_dP_gp = LocalVectorType::Zero(1);
		dmassXair_dX_gp = LocalVectorType::Zero(1);
		dmassXair_dT_gp = LocalVectorType::Zero(1);
		//Charact_gp = LocalVectorType::Zero(1);
		Lam_pm_gp = LocalVectorType::Zero(1);

		M = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//MASS MATRIX
		D = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//
		A1 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		A2 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		H = MathLib::LocalVector::Zero(n_dof / n_nodes);//for gravity term
		tmp = MathLib::LocalMatrix::Zero(1, 1);
		tmp_adv = MathLib::LocalMatrix::Zero(1, n_dim);
		localMass_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localAdvection_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localGravity_tmp = MathLib::LocalVector::Zero(n_nodes);//tmp matrix for gravity 
		for (j = 0; j < n_gsp; j++)
		{
			// calc on each gauss point
			// get the sampling point
			q->getSamplingPoint(j, gp_x);
			// compute basis functions
			fe->computeBasisFunctions(gp_x);
			// get the coordinates of current location 
			fe->getRealCoordinates(real_x);
			// transfer to Gauss point position
			NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
			double fac = geo_area*fe->getDetJ()*q->getWeight(j);
			// get the porosity
			pm->porosity->eval(gp_pos, poro);
			// get the intrinsic permeability
			pm->permeability->eval(gp_pos, K_perm);
			// evaluate viscosity of each fluid phase
			fluid1->dynamic_viscosity->eval(gp_pos, mu_G);
			fluid2->dynamic_viscosity->eval(gp_pos, mu_L);

			//fluid1->molar_mass->eval(gp_pos, M_G);
			//fluid2->molar_mass->eval(gp_pos, M_L);

			fluid1->density->eval(gp_pos, rho_G_std);
			fluid2->density->eval(gp_pos, rho_L_std);
			fluid1->specific_heat->eval(gp_pos, C_pg);
			fluid2->specific_heat->eval(gp_pos, C_pl);

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			solid->density->eval(gp_pos, rho_S_std);
			solid->specific_heat->eval(gp_pos, C_solid);

			// evaluation of the shape function
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();
			
			P_gp = Np*u1.head(n_nodes);// Gas phase pressure
			X_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);// rho_G^h S_G+rho_L^h S_L
			T_gp = Np*u1.tail(n_nodes);

			C_v = M_G / R / T_gp(0, 0);
			C_w = M_L / R / T_gp(0, 0);
			//+++++++++++++++++++++++++Calculate the secondary variable+++++++++++++++++++++++++++++
			Input(0) = P_gp(0, 0);
			Input(1) = X_gp(0, 0);
			Input(2) = T_gp(0, 0);
			//deliver the initial guess
			
			//_function_data->getS()->setIntegrationPointValue(ele_id, j, S_gp);
			_function_data->getS()->eval(gp_pos, test);
			Output(0) = test(0,0);
			_function_data->get_rhoLh()->eval(gp_pos, test);
			Output(1) = test(0, 0);
			_function_data->get_rhoGh()->eval(gp_pos, test);
			Output(2) = test(0, 0);

			_EOS->set_env_condition(Input);
 			_LP_EOS->solve(Input, Output);

			S_gp(0) = Output(0);
			rho_L_h_gp(0) = Output(1);
			rho_G_h_gp(0) = Output(2);

			_function_data->getS()->setIntegrationPointValue(ele_id, j, S_gp);			
			_function_data->get_rhoLh()->setIntegrationPointValue(ele_id, j, rho_L_h_gp);
			_function_data->get_rhoGh()->setIntegrationPointValue(ele_id, j, rho_G_h_gp);
			//+++++++++++++++++++++++++End Calculate++++++++++++++++++++++++++++++++++++++++++++++++
			PC_gp(0) = pm->getPc_bySat(S_gp(0));
			dPC_dSg_gp(0) = pm->Deriv_dPCdS(S_gp(0));
			P_Gh_gp(0) = _EOS->getPGH(P_gp(0, 0), T_gp(0, 0));
			//P_Gw_gp(0) = P_gp(0, 0) - P_Gh_gp(0);
			P_Gw_gp(0) = _EOS->getPGW(P_gp(0, 0), P_Gh_gp(0), T_gp(0, 0));
			rho_G_w_gp(0) = P_Gw_gp(0)*M_L / R / T_gp(0, 0);
			
			Xvap_gp(0) = P_Gw_gp(0) / P_gp(0, 0);//MOLAR FRACTION
			//massfractionXair_gp(0) = _EOS->get_massfraction(1 - Xvap_gp(0));
			massfractionXair_gp(0) = rho_G_h_gp(0) / (rho_G_w_gp(0) + rho_G_h_gp(0));
			EnthalpyH_G_gp(0) = C_pg*(T_gp(0, 0) - T0)*massfractionXair_gp(0) + (C_pl*(T_gp(0, 0) - T0) + delta_h_vap)*(1 - massfractionXair_gp(0));
			EnthalpyH_L_gp(0) = C_pl*(T_gp(0, 0) - T0);
			//+++++++++++++++++++++++++Calculate the derivatives +++++++++++++++++++++++++++++++++++
			dPGh_dP_gp(0) = _EOS->Deriv_PGH_dPG(P_gp(0, 0), P_Gh_gp(0), T_gp(0, 0));
			dPGh_dT_gp(0) = _EOS->Deriv_PGH_dT(P_gp(0, 0), P_Gh_gp(0), T_gp(0, 0));

			dPGw_dP_gp(0) = 1 - dPGh_dP_gp(0);
			dPGw_dT_gp(0) = -dPGh_dT_gp(0);
			
			drho_G_hdP_gp(0) = dPGh_dP_gp(0)*C_v;
			drho_G_h_dT_gp(0) = dPGh_dT_gp(0)*C_v - P_Gh_gp(0)*C_v / T_gp(0, 0);

			drho_L_h_dP_gp(0) = Hen*M_G*dPGh_dP_gp(0, 0);
			drho_L_h_dT_gp(0) = Hen*M_G*dPGh_dT_gp(0, 0);

			drho_G_w_dP_gp(0) = dPGw_dP_gp(0)*C_w;
			drho_G_w_dT_gp(0) = dPGw_dT_gp(0)*C_w - P_Gw_gp(0)*C_w/T_gp(0,0);

			dSgdX_gp(0) = _EOS->Deriv_dSgdX(S_gp(0), rho_L_h_gp(0), rho_G_h_gp(0), drho_L_h_dX_gp(0), drho_G_hdX_gp(0));
			dSgdP_gp(0) = _EOS->Deriv_dSgdP(S_gp(0), rho_L_h_gp(0), rho_G_h_gp(0), drho_L_h_dP_gp(0), drho_G_hdP_gp(0));
			dSgdT_gp(0) = _EOS->Deriv_dSgdT(S_gp(0), rho_L_h_gp(0), rho_G_h_gp(0), drho_L_h_dT_gp(0), drho_G_h_dT_gp(0));

			dmassXair_dP_gp(0) = M_G*M_L*(P_Gw_gp(0) - P_gp(0, 0)*dPGw_dP_gp(0)) / pow((M_G*P_gp(0, 0) + M_L*P_Gw_gp(0)), 2);
			dmassXair_dX_gp(0) = 0.0;
			dmassXair_dT_gp(0) = M_G*M_L*P_gp(0, 0)*dPGw_dT_gp(0) / pow((M_G*P_gp(0, 0) + M_L*P_Gw_gp(0)), 2);
			//+++++++++++++++++++++++++End Calculation++++++++++++++++++++++++++++++++++++++++++++++
			//+++++++++++++++++++++++++Calculate secondary variable for energy balance equations++++
			RHO_G = rho_G_h_gp(0) + rho_G_w_gp(0);// mass density of gas phase
			RHO_L = rho_l_std + rho_L_h_gp(0);
			Lam_pm_gp(0) = _EOS->get_overall_Heat_Capacity(S_gp(0));

			//+++++++++++++++++++++++++End calculate +++++++++++++++++++++++++++++++++++++++++++++++
			//Calc each entry of the mass matrix
			M(0, 0) = 0.0;
			M(0, 1) = poro;
			M(0, 2) = 0.0;
			M(1, 0) = poro*(S_gp(0)*drho_G_w_dP_gp(0) - (rho_l_std - rho_G_w_gp(0))*dSgdP_gp(0));
			M(1, 1) = poro*(1 + S_gp(0)*drho_G_w_dX_gp(0) - (rho_l_std - rho_G_w_gp(0))*dSgdX_gp(0));
			M(1, 2) = poro*(S_gp(0)*drho_G_w_dT_gp(0) - (rho_l_std - rho_G_w_gp(0))*dSgdT_gp(0));
			

			M(2, 0) = poro*S_gp(0)*(drho_G_hdP_gp(0) + drho_G_w_dP_gp(0))*EnthalpyH_G_gp(0) 
				+ poro*(1 - S_gp(0))*drho_L_h_dP_gp(0)*EnthalpyH_L_gp(0)
				+ poro*RHO_G*EnthalpyH_G_gp(0)*dSgdP_gp(0)
				- poro*RHO_L*EnthalpyH_L_gp(0)*dSgdP_gp(0)
				- poro*S_gp(0) - poro*P_gp(0, 0)*dSgdP_gp(0);
				

			M(2, 1) = poro*RHO_G*EnthalpyH_G_gp(0)*dSgdX_gp(0)
				- poro*RHO_L*EnthalpyH_L_gp(0)*dSgdX_gp(0)
				- poro*P_gp(0, 0)*dSgdX_gp(0);

			M(2, 2) = (1 - poro)*rho_S_std*C_solid
				+ poro*RHO_G*S_gp(0, 0)*(C_pg*massfractionXair_gp(0, 0) + C_pl*(1 - massfractionXair_gp(0, 0)))+ poro*RHO_L*C_pl*(1 - S_gp(0, 0))
				- poro*P_gp(0, 0)*dSgdT_gp(0)				
				+ poro*RHO_G*EnthalpyH_G_gp(0)*dSgdT_gp(0)
				- poro*RHO_L*EnthalpyH_L_gp(0)*dSgdT_gp(0);
				

			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp(0, 0) = M(ii, jj);
					localMass_tmp.setZero();
					fe->integrateWxN(j, tmp, localMass_tmp);
					localM.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localMass_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localM=" << std::endl;
			//std::cout << localM << std::endl;
			//-------------End debugging------------------------
			//calculate the relative permeability 
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0));
			lambda_L = K_perm*Kr_L_gp / mu_L;
			lambda_G = K_perm*Kr_G_gp / mu_G;
			//+++++++++++++++++++++++++Calculate the vecocity+++++++++++++++++++++++++++++++++++++++
			vel_L_gp = -lambda_L*((dNp)*u1.head(n_nodes)
				- dPC_dSg_gp(0)*(dSgdP_gp(0)*(dNp)*u1.head(n_nodes) + dSgdX_gp(0)*(dNp)*u1.block(n_nodes, 0, n_nodes, 1) + dSgdT_gp(0)*(dNp)*u1.tail(n_nodes)));
			vel_G_gp = -lambda_G*(dNp)*u1.head(n_nodes);
			//+++++++++++++++++++++++++ End Calculate+++++++++++++++++++++++++++++++++++++++
			isinf = _finite(dPC_dSg_gp(0));
			if (isinf==0)
			{
				dPC_dSg_gp(0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0) = dPC_dSg_gp(0);
			}

			A1(2, 2) = C_pl*RHO_L*vel_L_gp(0) + (C_pg*massfractionXair_gp(0) + C_pl*(1 - massfractionXair_gp(0))) *RHO_G*vel_G_gp(0);
			if (n_dim > 1){
				A2(2, 2) = C_pl*RHO_L*vel_L_gp(1) + (C_pg*massfractionXair_gp(0) + C_pl*(1 - massfractionXair_gp(0))) *RHO_G*vel_G_gp(1);
			}
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp_adv(0,0) = A1(ii, jj);
					if (n_dim > 1){
						tmp_adv(0, 1) = A2(ii, jj);
					}
					localAdvection_tmp.setZero();
					fe->integrateWxDN(j, tmp_adv , localAdvection_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) -= localAdvection_tmp;
				}
			}

			//Calc each entry of the Laplace Matrix
			D(0, 0) = (rho_G_h_gp(0)*lambda_G + rho_L_h_gp(0)*lambda_L) - lambda_L*rho_L_h_gp(0)*dPC_dSg_gp(0)*dSgdP_gp(0)
				+ poro*(1 - S_gp(0))*rho_l_std*(1/RHO_L)*D_L*drho_L_h_dP_gp(0) + poro*S_gp(0)*D_G*(1/RHO_G)*(rho_G_w_gp(0)*drho_G_hdP_gp(0)-rho_G_h_gp(0)*drho_G_w_dP_gp(0));
			D(0, 1) = -lambda_L*rho_L_h_gp(0)*dPC_dSg_gp(0)*dSgdX_gp(0);
			D(0, 2) = -lambda_L*rho_L_h_gp(0)*dPC_dSg_gp(0)*dSgdT_gp(0)
				+ poro*S_gp(0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_h_dT_gp(0) - rho_G_h_gp(0)*drho_G_w_dT_gp(0))
				+ poro*(1 - S_gp(0))*rho_l_std*(1 / RHO_L)*D_L*drho_L_h_dT_gp(0);
				
			D(1, 0) = (RHO_G*lambda_G + RHO_L*lambda_L) - lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdP_gp(0);
			D(1, 1) = -lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdX_gp(0);
			D(1, 2) = - lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdT_gp(0);

			D(2, 0) = (lambda_G*RHO_G*EnthalpyH_G_gp(0) + lambda_L*RHO_L*EnthalpyH_L_gp(0)) - lambda_L*RHO_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdP_gp(0);
			D(2, 1) = -lambda_L* RHO_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdX_gp(0);
			D(2, 2) = Lam_pm_gp(0)
				- lambda_L* RHO_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdT_gp(0);
			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp(0, 0) = D(ii, jj);
					localDispersion_tmp.setZero();
					fe->integrateDWxDN(j, tmp, localDispersion_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localDispersion_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localK=" << std::endl;
			//std::cout << localK << std::endl;
			//--------------end debugging-------------------
			//----------------assembly the gravity term--------------------------------
			

		}
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;

	}

private:
	FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data;
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];

	LocalVectorType Input, Output;

	LocalVectorType vec_g; // gravity term 
	LocalVectorType H;//this is for local gravity vector 
	LocalVectorType localGravity_tmp;

	LocalVectorType PC_gp;
	LocalVectorType P_Gw_gp;
	LocalVectorType P_Gh_gp;
	LocalVectorType rho_L_h_gp;
	LocalVectorType rho_G_h_gp;
	LocalVectorType rho_G_w_gp;
	LocalVectorType dPGh_dPG_gp;
	LocalVectorType dPC_dSg_gp;
	LocalVectorType PL_gp;
	LocalVectorType S_gp;
	LocalVectorType vel_L_gp;
	LocalVectorType vel_G_gp;
	LocalVectorType Lam_pm_gp;
	
	LocalMatrixType P_gp;
	LocalMatrixType X_gp;
	LocalMatrixType T_gp;

	LocalVectorType dSgdX_gp;
	LocalVectorType dSgdP_gp;
	LocalVectorType dSgdT_gp;
	LocalVectorType dPGdSg_gp;
	LocalVectorType drho_G_hdP_gp;
	LocalVectorType drho_G_hdX_gp;
	LocalVectorType drho_G_h_dT_gp;
	LocalVectorType drho_G_w_dT_gp;
	LocalVectorType drho_G_w_dP_gp;
	LocalVectorType drho_G_w_dX_gp;
	LocalVectorType drho_L_h_dP_gp;
	LocalVectorType drho_L_h_dX_gp;
	LocalVectorType drho_L_h_dT_gp;
	LocalVectorType dPGw_dP_gp;
	LocalVectorType dPGw_dX_gp;
	LocalVectorType dPGw_dT_gp;
	LocalVectorType dPGh_dP_gp;
	LocalVectorType dPGh_dX_gp;
	LocalVectorType dPGh_dT_gp;

	
	LocalVectorType Xvap_gp;
	LocalVectorType massfractionXair_gp;
	LocalVectorType dmassXair_dP_gp;
	LocalVectorType dmassXair_dX_gp;
	LocalVectorType dmassXair_dT_gp;
	
	LocalVectorType EnthalpyH_G_gp;
	LocalVectorType EnthalpyH_L_gp;

	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType A1;
	LocalMatrixType A2;
	LocalMatrixType tmp;//method 2
	LocalMatrixType tmp_adv;

	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localAdvection_tmp; //for method2
	LocalMatrixType localDispersion;
	LocalMatrixType matN;

	LocalMatrixType test;
	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;
	double C_pl;
	double C_pg;
	double C_solid;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w

	const double R = 8.314;
	const double C_h =0.0; 
	double C_v ;
	double C_w ;//M_L/RT
	const double rho_l_std = 1000.0;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	const double Hen = 1e-4;
	const double delta_h_vap = 2258000.0;//J/kg
	const double T0=273.15;
	double rho_S_std;
	double gamma_gp;
	std::size_t isinf;
	std::size_t Charact_func;
	EOS_NonIso_HeatPipeProb* _EOS;
	LocalProblem_EOS_NonIso_LocalNCP* _LP_EOS;//LOCAL PROBLEM
};



#endif  // end of ifndef