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

#ifndef NON_LINEAR_CMP_NONISO_CapPressureForm_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_NONISO_CapPressureForm_JACOBIAN_LOCAL_ASSEMBLER_H

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
#include "NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler.h"
#include "EOS_NonIso_HeatPipe_CapPressureForm.h"


template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>//,
class NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
    {
		_EOS = new EOS_NonIso_HeatPipe_CapPressureForm();
		
	};

	virtual ~NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler()
    {
		_function_data = NULL; 
		BaseLib::releaseObject(_EOS);
		
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
			vec_g = MathLib::LocalVector::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES
		*/
		/*
		*/
		Input = MathLib::LocalVector::Zero(3); //two primary variables
		Output = MathLib::LocalVector::Zero(2); //seven primary variables
		/*
		*Derivatives of the secondary variables
		*/
		// define the secondary variable value at specific one gauss point
		//Scalar value
		PG_gp = MathLib::LocalMatrix::Zero(1, 1);
		PC_gp = MathLib::LocalMatrix::Zero(1, 1);
		T_gp = MathLib::LocalMatrix::Zero(1, 1);

		dSgdPC_gp = MathLib::LocalMatrix::Zero(1, 1);
		dSgdT_gp = MathLib::LocalMatrix::Zero(1, 1);
		/*
		* vector of tmp variables
		*Omega, M and Charact
		*/
		S_gp = MathLib::LocalVector::Zero(1);
		PL_gp = MathLib::LocalVector::Zero(1);
		P_G_a_gp = MathLib::LocalVector::Zero(1);
		P_G_w_gp = MathLib::LocalVector::Zero(1);

		rho_G_a_gp = MathLib::LocalVector::Zero(1);
		rho_G_w_gp = MathLib::LocalVector::Zero(1);
		rho_G_gp = MathLib::LocalVector::Zero(1);
		vel_L_gp = MathLib::LocalVector::Zero(2);
		vel_G_gp = MathLib::LocalVector::Zero(2);
		vel_X_gp = MathLib::LocalVector::Zero(2);
		Lam_pm_gp = MathLib::LocalVector::Zero(1);

		drho_G_adPG_gp = MathLib::LocalVector::Zero(1);
		drho_G_adPC_gp = MathLib::LocalVector::Zero(1);
		drho_G_adT_gp = MathLib::LocalVector::Zero(1);

		drho_G_w_dPG_gp = MathLib::LocalVector::Zero(1);
		drho_G_w_dPC_gp = MathLib::LocalVector::Zero(1);
		drho_G_w_dT_gp = MathLib::LocalVector::Zero(1);

		dP_gw_dPC_gp = MathLib::LocalVector::Zero(1);
		dP_gw_dT_gp = MathLib::LocalVector::Zero(1);

		drho_GdT_gp = MathLib::LocalVector::Zero(1);
		drho_GdP_gp = MathLib::LocalVector::Zero(1);
		drho_GdPC_gp = MathLib::LocalVector::Zero(1);

		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively

		M = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		D = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2

		A1 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		A2 = MathLib::LocalMatrix::Zero(n_dof / n_nodes, n_dof / n_nodes);//method 2
		H = MathLib::LocalVector::Zero(n_dof / n_nodes);//for gravity term
		tmp = MathLib::LocalMatrix::Zero(1, 1);
		tmp_adv = MathLib::LocalMatrix::Zero(1, n_dim);
		localMass_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localAdvection_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localDispersion = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		localGravity_tmp = MathLib::LocalVector::Zero(n_nodes);//tmp matrix for gravity 

		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		Var_a = 0.0;
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
			MathLib::LocalMatrix &Np = *fe->getBasisFunction();
			MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();

			PG_gp = Np*u1.head(n_nodes);
			PC_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);
			T_gp = Np*u1.tail(n_nodes);
			C_v = M_G / R / T_gp(0, 0);
			C_w = M_L / R / T_gp(0, 0);
			//+++++++++++++++++++++++++Calculate the secondary variable+++++++++++++++++++++++++++++

			S_gp(0) = pm->getSat_byPC(PC_gp(0, 0));
			dSgdPC_gp(0) = pm->get_diriv_PC(PC_gp(0, 0));
			PL_gp(0) = PG_gp(0, 0) - PC_gp(0, 0);//calculate the liquid phase pressure

			P_G_w_gp(0) = _EOS->get_P_G_w(PC_gp(0, 0), T_gp(0, 0));//calculate the water vapor pressure
			P_G_a_gp(0) = PG_gp(0, 0) - P_G_w_gp(0);// partial pressure of air in gas phase
			if (P_G_a_gp(0) < 0)
				P_G_a_gp(0) = 0.0;

			rho_G_a_gp(0) = _EOS->get_RHO_G_a(P_G_a_gp(0),T_gp(0,0));// mass density of air in gas phase
			rho_G_w_gp(0) = _EOS->get_RHO_G_W(P_G_w_gp(0), T_gp(0, 0));// mas density of water vapor in gas phase
			C_pg = 1180 * P_G_a_gp(0, 0) / PG_gp(0) + 1996 * P_G_w_gp(0) / PG_gp(0);

			//+++++++++++++++++++++++++End Calculate++++++++++++++++++++++++++++++++++++++++++++++++

			//+++++++++++++++++++++++++Calculate the derivatives +++++++++++++++++++++++++++++++++++
			dSgdT_gp(0) = _EOS->Deriv_dSgdT(S_gp(0));

			drho_G_adPG_gp(0) = C_v;
			drho_G_adPC_gp(0) = _EOS->Deriv_drhoGa_dPC(PC_gp(0, 0), T_gp(0));
			drho_G_adT_gp(0) = _EOS->Deriv_drhoGa_dT(P_G_a_gp(0), PC_gp(0, 0), T_gp(0));
			
			if (P_G_a_gp(0) < 0.0){
				drho_G_adPG_gp(0) = 0.0;
				drho_G_adPC_gp(0) = 0.0;
				drho_G_adT_gp(0) = 0.0;
			}
			

			drho_G_w_dPG_gp(0) = 0.0;
			drho_G_w_dPC_gp(0) = _EOS->Deriv_drhoGW_dPC(PC_gp(0, 0), T_gp(0));
			drho_G_w_dT_gp(0) = _EOS->Deriv_drhoGW_dT(P_G_w_gp(0), PC_gp(0, 0), T_gp(0));

			dP_gw_dPC_gp(0) = _EOS->Deriv_dPgw_dPC(PC_gp(0, 0), T_gp(0));
			dP_gw_dT_gp(0) = _EOS->Deriv_dPgw_dT(PC_gp(0, 0), T_gp(0));

			drho_GdP_gp(0) = drho_G_adPG_gp(0) + drho_G_w_dPG_gp(0);
			drho_GdPC_gp(0) = drho_G_adPC_gp(0) + drho_G_w_dPC_gp(0);
			drho_GdT_gp(0) = drho_G_adT_gp(0) + drho_G_w_dT_gp(0);
			//+++++++++++++++++++++++++End Calculation++++++++++++++++++++++++++++++++++++++++++++++
			//+++++++++++++++++++++++++Calculate secondary variable for energy balance equations++++
			RHO_L =  rho_L_std;
			RHO_G = rho_G_a_gp(0) + rho_G_w_gp(0);
			Lam_pm_gp(0) = _EOS->get_overall_Heat_Capacity(S_gp(0));
			rho_G_gp(0) = RHO_G;// _EOS->getrho_G(PG_gp(0, 0), P_G_w_gp(0, 0), T_gp(0, 0));
			//+++++++++++++++++++++++++End calculate +++++++++++++++++++++++++++++++++++++++++++++++
			//Calc each entry of the mass matrix

			M(0, 0) = poro*(S_gp(0)*drho_G_adPG_gp(0));
			M(0, 1) = poro*(S_gp(0)*drho_G_adPC_gp(0) + rho_G_a_gp(0)*dSgdPC_gp(0));//0?
			M(0, 2) = poro*S_gp(0)*drho_G_adT_gp(0);

			M(1, 0) = poro*S_gp(0)*drho_GdP_gp(0);
			M(1, 1) = poro*((rho_G_gp(0) - rho_L_std)*dSgdPC_gp(0)) + poro*S_gp(0)*drho_G_adPC_gp(0);
			M(1, 2) = poro*(S_gp(0)*drho_GdT_gp(0));

			M(2, 0) = poro*C_pg*(T_gp(0, 0) - 273.15)*S_gp(0)*drho_GdP_gp(0);
			M(2, 1) = poro*(rho_G_gp(0)*C_pg - rho_L_std*C_pl)*(T_gp(0, 0) - 273.15)*dSgdPC_gp(0) + poro*(T_gp(0, 0) - 273.15)*S_gp(0)*drho_GdPC_gp(0)*C_pg;
			M(2, 2) = (1 - poro)*rho_S_std*C_solid + poro*(RHO_G*S_gp(0)*C_pg + RHO_L*(1 - S_gp(0))*C_pl)
				+ poro*S_gp(0)*C_pg*(T_gp(0, 0) - 273.15)*drho_GdT_gp(0);

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
			vel_G_gp = -lambda_G*(dNp)*u1.head(n_nodes);

			vel_L_gp = -(lambda_L*(dNp)*u1.head(n_nodes) - lambda_L*(dNp)*u1.block(n_nodes, 0, n_nodes, 1));

			vel_X_gp = -(rho_G_w_gp(0)*drho_G_adPG_gp(0)*(dNp)*u1.head(n_nodes)
				+ (rho_G_w_gp(0)*drho_G_adPC_gp(0) - rho_G_a_gp(0)*drho_G_w_dPC_gp(0))*(dNp)*u1.block(n_nodes, 0, n_nodes, 1)
				+ (rho_G_w_gp(0)*drho_G_adT_gp(0) - rho_G_a_gp(0)*drho_G_w_dT_gp(0))*(dNp)*u1.tail(n_nodes)
				);

			
			//+++++++++++++++++++++++++ End Calculate+++++++++++++++++++++++++++++++++++++++
			A1(0, 0) = 0.0;
			A1(0, 1) = 0.0;
			A1(0, 2) = 0.0;
			A1(1, 0) = 0.0;
			A1(1, 1) = 0.0;
			A1(1, 2) = 0.0;
			A1(2, 0) = 0.0;
			A1(2, 1) = 0.0;
			A1(2, 2) = C_pl*RHO_L*vel_L_gp(0) + C_pg*rho_G_gp(0, 0)*vel_G_gp(0) + 816 * poro*S_gp(0, 0)*D_G*(1 / RHO_G)*vel_X_gp(0);
			if (n_dim > 1){
				A2(2, 2) = C_pl*RHO_L*vel_L_gp(1) + C_pg*rho_G_gp(0, 0)*vel_G_gp(1) + 816 * poro*S_gp(0, 0)*D_G*(1 / RHO_G)*vel_X_gp(1);
			}
			


			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp_adv(0, 0) = A1(ii, jj);
					if (n_dim > 1){
						tmp_adv(0, 1) = A2(ii, jj);
					}
					localAdvection_tmp.setZero();
					fe->integrateWxDN(j, tmp_adv, localAdvection_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) -= localAdvection_tmp;
				}
			}
			//std::cout << "localK=" << std::endl;
			//std::cout << localK << std::endl;
			//Calc each entry of the advection matrix

			D(0, 0) = rho_G_a_gp(0)*lambda_G + poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adPG_gp(0));
			D(0, 1) = poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adPC_gp(0) - rho_G_a_gp(0)*drho_G_w_dPC_gp(0));
			D(0, 2) = poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adT_gp(0) - rho_G_a_gp(0)*drho_G_w_dT_gp(0));

			D(1, 0) = rho_L_std*lambda_L + rho_G_gp(0)*lambda_G;
			D(1, 1) = -rho_L_std*lambda_L;
			D(1, 2) = 0.0;

			D(2, 0) = rho_G_gp(0)*C_pg*(T_gp(0, 0) - 273.15)*lambda_G + rho_L_std*C_pl*(T_gp(0, 0) - 273.15)*lambda_L
				+ 816 * (T_gp(0, 0) - 273.15)*poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adPG_gp(0));
			D(2, 1) = -rho_L_std*C_pl*(T_gp(0, 0) - 273.15)*lambda_L
				+ 816 * (T_gp(0, 0) - 273.15)*poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adPC_gp(0) - rho_G_a_gp(0)*drho_G_w_dPC_gp(0));
			D(2, 2) = 2*Lam_pm_gp(0)
				+ 816 * (T_gp(0, 0) - 273.15)*poro*S_gp(0, 0)*D_G*(1 / RHO_G)*(rho_G_w_gp(0)*drho_G_adT_gp(0) - rho_G_a_gp(0)*drho_G_w_dT_gp(0));


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
			//H(0) = -rho_L_h_gp(0, 0)*RHO_L*lambda_L - rho_G_h_gp(0, 0)*RHO_G*lambda_G; //-pow(RHO_G, 2)*lambda_G - pow(RHO_L, 2)*lambda_L;
			//H(1) = -pow(RHO_G, 2)*lambda_G - pow(RHO_L, 2)*lambda_L;
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			/*
			if (hasGravityEffect) {
				// F += dNp^T * H* gz
				for (int idx = 0; idx < 2; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0, 0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;
				}
			}
			*/

		}
		//_function_data->get_elem_M_matrix()[ele_id] = localM;
		//_function_data->get_elem_K_matrix()[ele_id] = localK;

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

	double _Theta;

	std::size_t i, j, ii, jj;
	double real_x[3], gp_x[3];
	
	MathLib::LocalVector local_r_1;
	MathLib::LocalVector local_r_2;
	MathLib::LocalMatrix M1, M2, K1, K2, F1, F2, TMP_M_1, TMP_M_2;

	MathLib::LocalVector Input, Output;

	MathLib::LocalVector vec_g; // gravity term 
	MathLib::LocalVector H;//this is for local gravity vector 
	MathLib::LocalVector localGravity_tmp;

	MathLib::LocalVector rho_G_a_gp;
	MathLib::LocalVector rho_G_w_gp;
	MathLib::LocalVector rho_G_gp;
	MathLib::LocalVector PL_gp;
	MathLib::LocalVector S_gp;
	MathLib::LocalVector vel_L_gp;
	MathLib::LocalVector vel_G_gp;
	MathLib::LocalVector vel_X_gp;
	MathLib::LocalVector Lam_pm_gp;

	MathLib::LocalMatrix PG_gp;
	MathLib::LocalMatrix PC_gp;
	MathLib::LocalMatrix T_gp;

	MathLib::LocalVector P_G_w_gp;//partial pressure of water in gas phase
	MathLib::LocalVector P_G_a_gp;//partial pressure of air in gas phase

	MathLib::LocalVector drho_G_adPG_gp;
	MathLib::LocalVector drho_G_adPC_gp;
	MathLib::LocalVector drho_G_adT_gp;

	MathLib::LocalVector dSgdPC_gp;
	MathLib::LocalVector dSgdT_gp;

	MathLib::LocalVector drho_G_w_dT_gp;
	MathLib::LocalVector drho_G_w_dPG_gp;
	MathLib::LocalVector drho_G_w_dPC_gp;
	MathLib::LocalVector drho_GdT_gp;
	MathLib::LocalVector drho_GdPC_gp;
	MathLib::LocalVector drho_GdP_gp;

	MathLib::LocalVector dP_gw_dT_gp;
	MathLib::LocalVector dP_gw_dPC_gp;

	MathLib::LocalMatrix M;//Local Mass Matrix
	MathLib::LocalMatrix D;//local dispersion matrix
	MathLib::LocalMatrix A1;//local advection matrix
	MathLib::LocalMatrix A2;
	MathLib::LocalMatrix tmp;//method 2
	MathLib::LocalMatrix tmp_adv;

	MathLib::LocalMatrix localMass_tmp; //for method2
	MathLib::LocalMatrix localDispersion_tmp; //for method2
	MathLib::LocalMatrix localAdvection_tmp; //for method2
	MathLib::LocalMatrix localDispersion;
	MathLib::LocalMatrix matN;


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

	double RHO_L;// 
	double RHO_G;// 

	const double R = 8.314;
	double C_v;
	double C_w;//M_L/RT
	const double rho_l_std = 1000.0;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	double M_K;
	const double delta_h_vap = 2258000;//J/kg
	//double delta_h_vap = 0.0;
	double rho_S_std;
	EOS_NonIso_HeatPipe_CapPressureForm* _EOS;

};

#endif  // end of ifndef
