/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearCMP_2P2C_JacobianLocalAssembler.h
 *
 * 2015-03-21 by Yonghui Huang
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef NON_LINEAR_CMP_NONISOTOTALMASSEXPFORM_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_NONISOTOTALMASSEXPFORM_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "Ogs6FemData.h"
#include "NonLinearCMP_NonIsoTotalMassEXPForm_TimeODELocalAssembler.h"
#include "EOS_NonIsoTotalMassEXPForm.h"

//#include "LocalProblem_EOS_TotalMassEXPForm.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>//,
class NonLinearCMP_NonIsoTotalMassEXPForm_JacobianLocalAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_NonIsoTotalMassEXPForm_JacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
    {
		_EOS = new EOS_NonIsoTotalMassEXPForm();
		_LP_EOS = new LocalProblem_EOS_NonIso_TotalMassEXP();
	};

	virtual ~NonLinearCMP_NonIsoTotalMassEXPForm_JacobianLocalAssembler()
    {
		_function_data = NULL; 
		BaseLib::releaseObject(_EOS);
		BaseLib::releaseObject(_LP_EOS);
    };

	
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}
	
	void setTheta(double v)
	{
		assert(v >= .0 && v <= 1.0);
		_Theta = v;
	}
	/*
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		//clear the local Jacobian Matrix

		MathLib::LocalMatrix TMP_M(2*n_nodes,2* n_nodes);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		double dt = time.getTimeStepSize();
		LocalMatrixType _M = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
		LocalMatrixType _K = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
		_M = _function_data->get_elem_M_matrix().at(elem_id);
		_K = _function_data->get_elem_K_matrix().at(elem_id);
		double eps = 1e-7;// can also be defined as numerical way
		TMP_M = (1.0 / dt)*_M + _Theta*_K;
		//LocalMatrixType Local_LHS=(1/dt)*M+
		//LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
		//_local_assembler=assembleODE(time, e, u1, u0, M, K);
		for (size_t u_idx = 0; u_idx < 2 * n_nodes; u_idx++)
		{
			//clear the epsilon vectoe 
			LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));

			localJ.block(0, u_idx, u1.size(), 1) = (TMP_M*(u1 - e_vec) - TMP_M*(u1 + e_vec)) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}

		
    }  // end of function assembly
	*/
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
	{
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		const size_t n_dof = u1.size();
		double dt = time.getTimeStepSize();
		TMP_M_1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		TMP_M_2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		double eps = 1e-8;

		for (size_t u_idx = 0; u_idx < n_dof; u_idx++)
		{
			//clear M1,K1,F1,M2,K2,F2
			M1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F1 = MathLib::LocalVector::Zero(n_dof);
			M2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F2 = MathLib::LocalVector::Zero(n_dof);
			//clear local_r
			local_r_1 = MathLib::LocalVector::Zero(n_dof);
			local_r_2 = MathLib::LocalVector::Zero(n_dof);
			//clear the epsilon vectoe 
			MathLib::LocalVector e_vec = MathLib::LocalVector::Zero(n_dof);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
			assemble_jac_ODE(time, e, u1 - e_vec, u0, M1, K1, F1);
			TMP_M_1 = (1.0 / dt)*M1 + _Theta*K1;
			local_r_1 = TMP_M_1*(u1 - e_vec);
			TMP_M_1 = (1.0 / dt) * M1 - (1. - _Theta) * K1;
			local_r_1.noalias() -= TMP_M_1 * u0;
			local_r_1.noalias() -= F1;
			assemble_jac_ODE(time, e, u1 + e_vec, u0, M2, K2, F2);
			TMP_M_2 = (1.0 / dt)*M2 + _Theta*K2;
			local_r_2 = TMP_M_2*(u1 + e_vec);
			TMP_M_2 = (1.0 / dt) * M2 - (1. - _Theta) * K2;
			local_r_2.noalias() -= TMP_M_2 * u0;
			local_r_2.noalias() -= F2;
			localJ.block(0, u_idx, u1.size(), 1) = (local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}
	}
protected:
	virtual void assemble_jac_ODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix  & localM, MathLib::LocalMatrix & localK, MathLib::LocalMatrix & localF)
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
		const bool hasGravityEffect = _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES
		*/
		/*
		*/
		Input = LocalVectorType::Zero(3); //two primary variables
		Output = LocalVectorType::Zero(2); //seven primary variables
		/*
		*Derivatives of the secondary variables
		*/
		// define the secondary variable value at specific one gauss point
		//Scalar value

		dSgdX_gp = LocalMatrixType::Zero(1, 1);
		dSgdP_gp = LocalMatrixType::Zero(1, 1);
		dSgdT_gp = LocalMatrixType::Zero(1, 1);
		dPC_dSg_gp = LocalMatrixType::Zero(1, 1);
		/*
		* vector of tmp variables
		*Omega, M and Charact
		*/
		S_gp = LocalVectorType::Zero(1);
		PL_gp = LocalVectorType::Zero(1);
		PC_gp = LocalVectorType::Zero(1);
		PGH_gp = LocalVectorType::Zero(1);
		PGW_gp = LocalVectorType::Zero(1);
		rho_G_h_gp = LocalVectorType::Zero(1);
		rho_G_w_gp = LocalVectorType::Zero(1);
		dPGw_dPG_gp = LocalVectorType::Zero(1);
		dPGw_dT_gp = LocalVectorType::Zero(1);
		drho_G_hdP_gp = LocalVectorType::Zero(1);
		drho_G_hdX_gp = LocalVectorType::Zero(1);
		drho_G_hdT_gp = LocalVectorType::Zero(1);
		drho_G_wdP_gp = LocalVectorType::Zero(1);
		drho_G_wdX_gp = LocalVectorType::Zero(1);
		drho_G_wdT_gp = LocalVectorType::Zero(1);
		dPGh_dPG_gp = LocalVectorType::Zero(1);
		alpha_gp = LocalVectorType::Zero(1);
		beta_gp = LocalVectorType::Zero(1);
		Charact_gp = LocalVectorType::Zero(1);
		massfractionXair_gp = LocalVectorType::Zero(1);
		EnthalpyH_G_gp = LocalVectorType::Zero(1);
		EnthalpyH_L_gp = LocalVectorType::Zero(1);
		Lam_pm_gp = LocalVectorType::Zero(1);

		vel_L_gp = LocalVectorType::Zero(1);
		vel_G_gp = LocalVectorType::Zero(1);
		vel_X_air_gp = LocalVectorType::Zero(1);
		vel_X_vap_gp = LocalVectorType::Zero(1);
		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively

		M = LocalMatrixType::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		D = LocalMatrixType::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		A1 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		A2 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		H = LocalVectorType::Zero(n_dof / n_nodes);//for gravity term
		tmp = LocalMatrixType::Zero(1, 1);
		tmp_adv = MathLib::LocalMatrix::Zero(1, n_dim);
		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(n_dof, n_dof);
		localAdvection_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localGravity_tmp = LocalVectorType::Zero(n_nodes);//tmp matrix for gravity 
		test = LocalMatrixType::Zero(1, 1);

		
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		Var_a = 0.0;
		gamma_gp = 0.0;
		isinf = 0;
		Charact_func = 1;
		beta_func = 1;
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
			fluid1->thermal_conductivity->eval(gp_pos, lambda_n);
			fluid2->thermal_conductivity->eval(gp_pos, lambda_w);
			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			solid->density->eval(gp_pos, rho_S_std);
			solid->specific_heat->eval(gp_pos, C_solid);

			// evaluation of the shape function
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();


			RHO_L = 0.0;
			RHO_G = 0.0;


			P_gp = Np*u1.head(n_nodes);// assume to be gas phase pressure for simplicity
			X_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);// assume to be the total mass density of heavy component
			T_gp = Np*u1.tail(n_nodes);// assume to be the total mass density of light component
			//+++++++++++++++++++++++++Define constant values+++++++++++++++++++++++++++++
			C_v = M_G / R / T_gp(0, 0);
			C_w = M_L / R / T_gp(0, 0);
			//+++++++++++++++++++++++++Calculate the secondary variable+++++++++++++++++++++++++++++
			Input(0) = P_gp(0, 0);
			Input(1) = X_gp(0, 0);
			Input(2) = T_gp(0, 0);
			//deliver the initial guess

			//_function_data->getS()->setIntegrationPointValue(ele_id, j, S_gp);
			_function_data->getS()->eval(gp_pos, test);
			Output(0) = test(0, 0);
			_function_data->get_rhoGw()->eval(gp_pos, test);
			Output(1) = test(0, 0);

			_EOS->set_env_condition(Input);
			_LP_EOS->solve(Input, Output);

			S_gp(0) =  Output(0);
			rho_G_w_gp(0) = Output(1);

			/*if (rho_G_w_gp(0)<0)
				Charact_func = 0;
			else
				Charact_func = 1;
				*/
			rho_G_h_gp(0) = _EOS->get_RHO_G_H(P_gp(0, 0), T_gp(0, 0), rho_G_w_gp(0));
			PGW_gp(0) = _EOS->get_P_sat(P_gp(0, 0), T_gp(0, 0));
			if (_EOS->get_Vapor_Pressure(T_gp(0, 0)) > P_gp(0, 0)) {
				beta_func = 1;
			}
			else
				beta_func = 1;
			massfractionXair_gp(0) = rho_G_h_gp(0) / (rho_G_w_gp(0) + rho_G_h_gp(0));
			EnthalpyH_G_gp(0) = C_pg*(T_gp(0, 0) - T0)*massfractionXair_gp(0) + (C_pl*(T_gp(0, 0) - T0) + delta_h_vap)*(1 - massfractionXair_gp(0));
			EnthalpyH_L_gp(0) = C_pl*(T_gp(0, 0) - T0);
			//+++++++++++++++++++++++++End Calculate++++++++++++++++++++++++++++++++++++++++++++++++
			//+++++++++++++++++++++++++Calculate the characteristic function++++++++++++++++++++++++
			PC_gp(0) = pm->getPc_bySat(S_gp(0));
			dPC_dSg_gp(0) = pm->Deriv_dPCdS(S_gp(0));
			PL_gp(0) = P_gp(0, 0) - PC_gp(0);
			
			//+++++++++++++++++++++++++Calculate the derivatives +++++++++++++++++++++++++++++++++++
			dPGw_dPG_gp(0) = 0.0;
			dPGw_dT_gp(0) = _EOS->Deriv_dPsat_dT(P_gp(0, 0), T_gp(0, 0))*beta_func;

			drho_G_wdP_gp(0) = _EOS->Deriv_drhoGw_dP(dPGh_dPG_gp(0), Charact_func);
			drho_G_wdX_gp(0) = _EOS->Deriv_drhoGw_dX(Charact_func);
			drho_G_wdT_gp(0) = _EOS->Deriv_drhoGw_dT(PGW_gp(0), dPGw_dT_gp(0), Charact_func);

			dSgdX_gp(0) = _EOS->Deriv_dSgdX(rho_G_w_gp(0), drho_G_wdX_gp(0), S_gp(0));
			dSgdP_gp(0) = _EOS->Deriv_dSgdP(rho_G_w_gp(0), drho_G_wdP_gp(0), S_gp(0));
			dSgdT_gp(0) = _EOS->Deriv_dSgdT(rho_G_w_gp(0), drho_G_wdT_gp(0), S_gp(0));

			drho_G_hdP_gp(0) = _EOS->Deriv_drhoGh_dP(drho_G_wdP_gp(0), T_gp(0, 0));
			drho_G_hdX_gp(0) = _EOS->Deriv_drhoGh_dX(drho_G_wdX_gp(0));
			drho_G_hdT_gp(0) = _EOS->Deriv_drhoGh_dT(drho_G_wdT_gp(0), P_gp(0, 0), T_gp(0, 0));
			//+++++++++++++++++++++++++End Calculation++++++++++++++++++++++++++++++++++++++++++++++
			//Calc each entry of the mass matrix
			RHO_L = rho_L_std;
			RHO_G = rho_G_h_gp(0) + rho_G_w_gp(0);
			Lam_pm_gp(0) = _EOS->get_overall_Heat_Capacity(S_gp(0), lambda_n, lambda_w);

			M(0, 0) = poro*(drho_G_hdP_gp(0)*S_gp(0) + rho_G_h_gp(0)*dSgdP_gp(0));
			M(0, 1) = poro*(1 + drho_G_hdX_gp(0)*S_gp(0) + rho_G_h_gp(0)*dSgdX_gp(0));
			M(0, 2) = poro*(drho_G_hdT_gp(0)*S_gp(0) + rho_G_h_gp(0)*dSgdT_gp(0));
			M(1, 0) = 0.0;
			M(1, 1) = poro;
			M(1, 2) = 0.0;

			M(2, 0) = poro*S_gp(0)*(drho_G_hdP_gp(0) + drho_G_wdP_gp(0))*EnthalpyH_G_gp(0)
				+ poro*RHO_G*EnthalpyH_G_gp(0)*dSgdP_gp(0)
				- poro*RHO_L*EnthalpyH_L_gp(0)*dSgdP_gp(0)
				+ poro*RHO_G*S_gp(0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdP_gp(0) - rho_G_h_gp(0)*drho_G_wdP_gp(0)) / pow(RHO_G, 2)
				- poro*S_gp(0) - poro*P_gp(0, 0)*dSgdP_gp(0);
			M(2, 1) = poro*RHO_G*EnthalpyH_G_gp(0)*dSgdX_gp(0)
				- poro*RHO_L*EnthalpyH_L_gp(0)*dSgdX_gp(0)
				+ poro*RHO_G*S_gp(0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdX_gp(0) - rho_G_h_gp(0)*drho_G_wdX_gp(0)) / pow(RHO_G, 2)
				- poro*P_gp(0, 0)*dSgdX_gp(0);
			M(2, 2) = (1 - poro)*rho_S_std*C_solid + poro*S_gp(0)*(drho_G_hdT_gp(0) + drho_G_wdT_gp(0))*EnthalpyH_G_gp(0)
				+ poro*RHO_G*S_gp(0, 0)*(C_pg*massfractionXair_gp(0, 0) + C_pl*(1 - massfractionXair_gp(0, 0))) + poro*RHO_L*C_pl*(1 - S_gp(0, 0))
				+ poro*RHO_G*EnthalpyH_G_gp(0)*dSgdT_gp(0) - poro*RHO_L*EnthalpyH_L_gp(0)*dSgdT_gp(0)
				+ poro*RHO_G*S_gp(0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdT_gp(0) - rho_G_h_gp(0)*drho_G_wdT_gp(0)) / pow(RHO_G, 2)
				- poro*P_gp(0, 0)*dSgdT_gp(0);

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
			vel_X_air_gp = -(dNp)*((rho_G_w_gp(0)*drho_G_hdP_gp(0) - rho_G_h_gp(0)*drho_G_wdP_gp(0))*u1.head(n_nodes)
				+ (rho_G_w_gp(0)*drho_G_hdX_gp(0) - rho_G_h_gp(0)*drho_G_wdX_gp(0))*u1.block(n_nodes, 0, n_nodes, 1)
				+ (rho_G_w_gp(0)*drho_G_hdT_gp(0) - rho_G_h_gp(0)*drho_G_wdT_gp(0))*u1.tail(n_nodes));
			vel_X_vap_gp = -vel_X_air_gp;
			//+++++++++++++++++++++++++ End Calculate+++++++++++++++++++++++++++++++++++++++
			isinf = _finite(dPC_dSg_gp(0));
			if (isinf == 0)
			{
				dPC_dSg_gp(0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0) = dPC_dSg_gp(0);
			}

			A1(2, 2) = C_pl*RHO_L*vel_L_gp(0) + (C_pg*massfractionXair_gp(0) + C_pl*(1 - massfractionXair_gp(0))) *RHO_G*vel_G_gp(0)
				+ C_pg*poro*S_gp(0)*(1 / RHO_G)*D_L*vel_X_air_gp(0)
				+ C_pl*poro*S_gp(0)*(1 / RHO_G)*D_L*vel_X_vap_gp(0);
			if (n_dim > 1){
				A2(2, 2) = C_pl*RHO_L*vel_L_gp(1) + (C_pg*massfractionXair_gp(0) + C_pl*(1 - massfractionXair_gp(0))) *RHO_G*vel_G_gp(1)
					+ C_pg*poro*S_gp(0)*(1 / RHO_G)*D_L*vel_X_air_gp(1)
					+ C_pl*poro*S_gp(0)*(1 / RHO_G)*D_L*vel_X_vap_gp(1);
			}
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp_adv(0, 0) = A1(ii, jj);
					if (n_dim > 1){
						tmp_adv(0, 1) = A2(ii, jj);
					}
					localAdvection_tmp.setZero();
					fe->integrateWxDN(j, tmp_adv, localAdvection_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localAdvection_tmp;
				}
			}

			//Calc each entry of the Laplace Matrix
			D(0, 0) = lambda_G*RHO_G + lambda_L*RHO_L - lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdP_gp(0);
			D(0, 1) = -lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdX_gp(0);
			D(0, 2) = -lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdT_gp(0);

			D(1, 0) = RHO_L*lambda_L + rho_G_w_gp(0)*lambda_G - lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdP_gp(0)
				- poro*S_gp(0)*D_L*(rho_G_w_gp(0)*drho_G_hdP_gp(0) - rho_G_h_gp(0)*drho_G_wdP_gp(0)) / RHO_G;

			D(1, 1) = -lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdX_gp(0)
				- poro*S_gp(0)*D_L*(rho_G_w_gp(0)*drho_G_hdX_gp(0) - rho_G_h_gp(0)*drho_G_wdX_gp(0)) / RHO_G;
			D(1, 2) = -lambda_L*RHO_L*dPC_dSg_gp(0)*dSgdT_gp(0)
				- poro*S_gp(0)*D_L*(rho_G_w_gp(0)*drho_G_hdT_gp(0) - rho_G_h_gp(0)*drho_G_wdT_gp(0)) / RHO_G;;

			D(2, 0) = RHO_L*lambda_L*EnthalpyH_L_gp(0) + RHO_G*lambda_G*EnthalpyH_G_gp(0)
				- RHO_L*lambda_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdP_gp(0)
				+ poro*S_gp(0)*D_L*(1 / RHO_G)*C_pg*(T_gp(0, 0) - T0)*(rho_G_w_gp(0)*drho_G_hdP_gp(0) - rho_G_h_gp(0)*drho_G_wdP_gp(0))
				- poro*S_gp(0)*D_L*(1 / RHO_G)*(C_pl*(T_gp(0, 0) - T0) + delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdP_gp(0) - rho_G_h_gp(0)*drho_G_wdP_gp(0));
			D(2, 1) = -RHO_L*lambda_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdX_gp(0)
				+ poro*S_gp(0)*D_L*(1 / RHO_G)*C_pg*(T_gp(0, 0) - T0)*(rho_G_w_gp(0)*drho_G_hdX_gp(0) - rho_G_h_gp(0)*drho_G_wdX_gp(0))
				- poro*S_gp(0)*D_L*(1 / RHO_G)*(C_pl*(T_gp(0, 0) - T0) + delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdX_gp(0) - rho_G_h_gp(0)*drho_G_wdX_gp(0));
			D(2, 2) = Lam_pm_gp(0)
				- lambda_L* RHO_L*EnthalpyH_L_gp(0)*dPC_dSg_gp(0)*dSgdT_gp(0)
				+ poro*S_gp(0)*D_L*(1 / RHO_G)*C_pg*(T_gp(0, 0) - T0)*(rho_G_w_gp(0)*drho_G_hdT_gp(0) - rho_G_h_gp(0)*drho_G_wdT_gp(0))
				- poro*S_gp(0)*D_L*(1 / RHO_G)*(C_pl*(T_gp(0, 0) - T0) + delta_h_vap)*(rho_G_w_gp(0)*drho_G_hdT_gp(0) - rho_G_h_gp(0)*drho_G_wdT_gp(0));


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
			H(0) = -RHO_L*RHO_L*lambda_L - RHO_G*RHO_G*lambda_G; //-pow(RHO_G, 2)*lambda_G - pow(RHO_L, 2)*lambda_L;
			H(1) = -rho_G_w_gp(0)*RHO_G*lambda_G - rho_L_std*RHO_L*lambda_L;
			H(2) = -RHO_L*RHO_L*lambda_L*EnthalpyH_L_gp(0) - RHO_G*RHO_G*lambda_G*EnthalpyH_G_gp(0);
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			if (hasGravityEffect) {
				// F += dNp^T * H* gz

				for (int idx = 0; idx < 3; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0, 0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;
				}
			}

		}

	}
private:
    /**
      * FEM object
      */ 
    FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data; 
	//LocalMatrixType M = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	//LocalMatrixType K = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	double _Theta;
	//LocalMatrixType* _M;
	//LocalMatrixType* _K;
	std::size_t i, j, ii, jj;
	double real_x[3], gp_x[3];
	
	MathLib::LocalVector local_r_1;
	MathLib::LocalVector local_r_2;
	MathLib::LocalMatrix M1, M2, K1, K2, F1, F2, TMP_M_1, TMP_M_2;

	MathLib::LocalVector Input, Output;

	LocalVectorType vec_g; // gravity term 
	LocalVectorType H;//this is for local gravity vector 
	LocalVectorType localGravity_tmp;

	LocalVectorType PC_gp;
	//LocalVectorType PL_gp;
	LocalVectorType PGH_gp;
	LocalVectorType PGW_gp;
	LocalVectorType rho_L_h_gp;
	LocalVectorType rho_G_h_gp;
	LocalVectorType rho_G_w_gp;
	LocalVectorType dPGh_dPG_gp;
	LocalVectorType dPGw_dPG_gp;
	LocalVectorType dPGw_dT_gp;
	LocalVectorType dPC_dSg_gp;
	LocalVectorType PL_gp;
	LocalVectorType S_gp;


	LocalVectorType vel_L_gp;
	LocalVectorType vel_G_gp;
	LocalVectorType vel_X_air_gp;
	LocalVectorType vel_X_vap_gp;

	LocalVectorType Lam_pm_gp;

	LocalMatrixType P_gp;
	LocalMatrixType X_gp;
	LocalMatrixType T_gp;

	LocalVectorType dSgdX_gp;
	LocalVectorType dSgdP_gp;
	LocalVectorType dSgdT_gp;

	LocalVectorType drho_G_hdP_gp;
	LocalVectorType drho_G_hdX_gp;
	LocalVectorType drho_G_hdT_gp;

	LocalVectorType drho_G_wdP_gp;
	LocalVectorType drho_G_wdX_gp;
	LocalVectorType drho_G_wdT_gp;

	LocalVectorType alpha_gp;
	LocalVectorType beta_gp;
	LocalVectorType Charact_gp;

	LocalVectorType massfractionXair_gp;
	LocalVectorType EnthalpyH_G_gp;
	LocalVectorType EnthalpyH_L_gp;

	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType tmp;//method 2
	LocalMatrixType A1;
	LocalMatrixType A2;
	LocalMatrixType tmp_adv;
	//LocalMatrixType tmp;//method 2

	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localDispersion;
	LocalMatrixType localAdvection_tmp; //for method2
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
	double C_v;
	double C_w;
	const double R = 8.314;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	const double rho_l_std = 1000.0;
	const double delta_h_vap = 2258000.0;//J/kg
	const double T0 = 273.15;
	double gamma_gp;
	double rho_S_std;
	double lambda_w;
	double lambda_n;
	std::size_t isinf;
	std::size_t Charact_func;
	std::size_t beta_func;
	EOS_NonIsoTotalMassEXPForm* _EOS;
	LocalProblem_EOS_NonIso_TotalMassEXP* _LP_EOS;//LOCAL PROBLEM
};

#endif  // end of ifndef
