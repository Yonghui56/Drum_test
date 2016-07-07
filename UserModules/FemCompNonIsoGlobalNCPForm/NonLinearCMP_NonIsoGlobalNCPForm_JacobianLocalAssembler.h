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

#ifndef NON_LINEAR_CMP_GlobalComplementaryForm_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_GlobalComplementaryForm_JACOBIAN_LOCAL_ASSEMBLER_H

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
#include "NonLinearCMP_NonIsoGlobalNCPForm_TimeODELocalAssembler.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA >//,
class NonLinearCMP_NonIsoGlobalNCPForm_JacobianLocalAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_NonIsoGlobalNCPForm_JacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
    {
		_EOS = new EOS_NonIsoGlobalNCPForm();
	};

	virtual ~NonLinearCMP_NonIsoGlobalNCPForm_JacobianLocalAssembler()
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
	/*
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		//clear the local Jacobian Matrix

		MathLib::LocalMatrix TMP_M(3*n_nodes,3* n_nodes);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		double dt = time.getTimeStepSize();
		LocalMatrixType _M = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		LocalMatrixType _K = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		LocalMatrixType _J = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		_M = _function_data->get_elem_M_matrix().at(elem_id);
		_K = _function_data->get_elem_K_matrix().at(elem_id);
		//_J = _function_data->get_elem_J_matrix().at(elem_id);
		
		double eps = 1e-9;// can also be defined as numerical way
		TMP_M = (1.0 / dt)*_M + _Theta*_K;
		//std::cout << TMP_M;
		//LocalMatrixType Local_LHS=(1/dt)*M+
		//LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
		//_local_assembler=assembleODE(time, e, u1, u0, M, K);
		for (size_t u_idx = 0; u_idx < u1.size(); u_idx++)
		{
			//clear the epsilon vectoe 
			LocalVectorType e_vec = LocalVectorType::Zero(3 * n_nodes);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
			localJ.block(0, u_idx, u1.size(), 1) = (TMP_M*(u1 - e_vec) - TMP_M*(u1 + e_vec)) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}

		//localJ.noalias() -= _J;
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
		//FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		//MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		double eps = 1e-7;

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
			//e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
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
			localJ.block(0, u_idx, u1.size(), 1) += (local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(u1(u_idx)));
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

		P_gp = MathLib::LocalMatrix::Zero(1, 1);
		X_G_a_gp = MathLib::LocalMatrix::Zero(1, 1);
		X_G_w_gp = MathLib::LocalMatrix::Zero(1, 1);
		S_gp = MathLib::LocalMatrix::Zero(1, 1);
		T_gp = MathLib::LocalMatrix::Zero(1, 1);

		/*
		*SECONDARY VARIABLES initialize
		*/
		PC = MathLib::LocalVector::Zero(n_nodes);
		rho_L_h = MathLib::LocalVector::Zero(n_nodes);
		rho_G_h = MathLib::LocalVector::Zero(n_nodes);

		dPC_dSg = MathLib::LocalVector::Zero(n_nodes);
		dMassfraction_air = MathLib::LocalVector::Zero(n_nodes);
		dMassfraction_vap = MathLib::LocalVector::Zero(n_nodes);
		rho_mol_G_gp = LocalMatrixType::Zero(1, 1);
		rho_mol_L_gp = LocalMatrixType::Zero(1, 1);
		rho_mass_G_gp = LocalMatrixType::Zero(1, 1);
		rho_mass_L_gp = LocalMatrixType::Zero(1, 1);

		PC_gp = LocalMatrixType::Zero(1, 1);
		dPC_dSg_gp = LocalMatrixType::Zero(1, 1);
		P_G_a_gp = LocalMatrixType::Zero(1, 1);
		P_G_w_gp = LocalMatrixType::Zero(1, 1);
		massfractionX_gp = MathLib::LocalMatrix::Zero(1, 1);
		Xvap_gp = MathLib::LocalMatrix::Zero(1, 1);
		h_G_gp = MathLib::LocalMatrix::Zero(1, 1);
		dP_gw_dT_gp = LocalMatrixType::Zero(1, 1);
		dP_gw_dPC_gp = LocalMatrixType::Zero(1, 1);
		drho_mol_GdT_gp = LocalMatrixType::Zero(1, 1);
		drho_mol_LdT_gp = LocalMatrixType::Zero(1, 1);
		Lam_pm_gp = LocalMatrixType::Zero(1, 1);

		PC_node = MathLib::LocalVector::Zero(n_nodes);//
		vel_L_gp = MathLib::LocalVector::Zero(n_dim);
		vel_G_gp = MathLib::LocalVector::Zero(n_dim);
		vel_X_air_gp = MathLib::LocalVector::Zero(n_dim);
		vel_X_vap_gp = MathLib::LocalVector::Zero(n_dim);
		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively

		M = MathLib::LocalMatrix::Zero(3, 5);//method 2
		D = MathLib::LocalMatrix::Zero(3, 5);//method 2
		H = MathLib::LocalVector::Zero(3);//for gravity term
		A1 = MathLib::LocalMatrix::Zero(n_dof / n_nodes - 1, n_dof / n_nodes);//method 2
		A2 = MathLib::LocalMatrix::Zero(n_dof / n_nodes - 1, n_dof / n_nodes);//method 2
		tmp = MathLib::LocalMatrix::Zero(1, 1);
		tmp_adv = MathLib::LocalMatrix::Zero(1, n_dim);
		Input = MathLib::LocalVector::Zero(5);
		localMass_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localDispersion = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
		localAdvection_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		LocalJacb = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		localGravity_tmp = MathLib::LocalVector::Zero(n_nodes);//tmp vect for gravity 
		localRes = MathLib::LocalVector::Zero(n_nodes);//tmp vect for gravity 
		for (i = 0; i < n_nodes; i++)
		{
			node_id = e.getNodeID(i); // e is the current element

			PC_node[i] = _function_data->get_PC()->getValue(node_id);
		}
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		flag = 1.0;
		//-------------------Begin assembly on each gauss point---------------------------
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

			solid->density->eval(gp_pos, rho_S_std);
			solid->specific_heat->eval(gp_pos, C_solid);

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);


			// evaluation of the shape function
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();


			P_gp = Np*u1.head(n_nodes);//gas pressure
			X_G_a_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);//molar fraction of air in gas phase
			X_G_w_gp = Np*u1.block(2 * n_nodes, 0, n_nodes, 1);//molar fraction of air in gas phase
			S_gp = Np*u1.block(3 * n_nodes, 0, n_nodes, 1);// gas saturation
			T_gp = Np*u1.tail(n_nodes);

			/*if (X_gp(0, 0) < 0){
			X_gp(0, 0) = 0.05;
			flag = 0.0;
			}*/

			C_v = 1 / R / T_gp(0, 0);

			PC_gp(0, 0) = pm->getPc_bySat(S_gp(0, 0));
			dPC_dSg_gp(0, 0) = pm->Deriv_dPCdS(S_gp(0, 0));// (S_gp(0, 0));

			P_G_w_gp(0, 0) = _EOS->get_P_G_w(P_gp(0, 0), X_G_a_gp(0,0), PC_gp(0, 0), T_gp(0, 0));
			//Xvap_gp(0, 0) = 1 - X_gp(0, 0);// P_G_w_gp(0, 0) / P_gp(0, 0);// GIVE the molar fraction of water vapor in gas phase
			P_G_a_gp(0, 0) = P_gp(0, 0)*X_G_a_gp(0, 0);// 
			massfractionX_gp(0, 0) = _EOS->get_massfraction(X_G_a_gp(0, 0), X_G_w_gp(0, 0)); //(1-Xvap(0,0))indicates the molar fraction of the air in the gas phase 
			dMassfraction_air(0, 0) = _EOS->Deriv_massfraction_a(X_G_a_gp(0, 0), X_G_w_gp(0, 0));//!!!!
			dMassfraction_vap(0, 0) = _EOS->Deriv_massfraction_w(X_G_a_gp(0, 0), X_G_w_gp(0, 0));//!!!!
			h_G_gp(0, 0) = (C_pg*(T_gp(0, 0) - 273.15)*massfractionX_gp(0, 0) + (C_pl*(T_gp(0, 0) - 273.15) + delta_h_vap)*(1 - massfractionX_gp(0, 0)));
			//C_pg = 11.80 * X_gp(0, 0) + 19.96 * (1 - X_gp(0, 0));

			rho_mol_G_gp(0, 0) = P_gp(0, 0)*C_v;
			rho_mol_L_gp(0, 0) = rho_L_std / M_L;

			rho_mass_G_gp(0, 0) = rho_mol_G_gp(0, 0)*(X_G_a_gp(0, 0)*M_G + X_G_w_gp(0, 0)*M_L);
			//rho_mass_G_gp(0, 0) = rho_mol_G_gp(0, 0)*X_gp(0, 0)*M_G + P_G_w_gp(0, 0)* M_L*C_v;
			rho_mass_L_gp(0, 0) = rho_L_std;

			dP_gw_dPC_gp(0) = _EOS->Deriv_dPgw_dPC(P_gp(0, 0), PC_gp(0, 0), T_gp(0, 0));

			dP_gw_dT_gp(0) = _EOS->Deriv_dPgw_dT(P_gp(0, 0), PC_gp(0, 0), T_gp(0, 0));

			drho_mol_GdT_gp(0) = 0.0;// -P_gp(0, 0)*C_v / T_gp(0, 0);
			//drho_mol_LdT_gp(0) = -(P_gp(0, 0) - PC_gp(0, 0))*C_v / T_gp(0, 0);

			Lam_pm_gp(0) = _EOS->get_overall_Heat_Capacity(S_gp(0), lambda_n, lambda_w);

			isinf = _finite(dPC_dSg_gp(0, 0));
			if (isinf == 0)
			{
				dPC_dSg_gp(0, 0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0, 0) = dPC_dSg_gp(0, 0);
			}

			M(0, 0) = poro*X_G_a_gp(0, 0)*S_gp(0, 0)*C_v;//dPG
			M(0, 1) = poro*S_gp(0, 0)*rho_mol_G_gp(0, 0);//dX_G_a			
			M(0, 2) = 0.0;//dX_G_w
			M(0, 3) = poro*X_G_a_gp(0, 0)*rho_mol_G_gp(0, 0);//dSG
			M(0, 4) = poro*X_G_a_gp(0, 0)*S_gp(0, 0)*drho_mol_GdT_gp(0);//dT

			M(1, 0) = poro*S_gp(0, 0)*C_v*X_G_w_gp(0, 0);//dPG
			M(1, 1) = 0.0;//dX_G_a
			M(1, 2) = poro*S_gp(0, 0)*rho_mol_G_gp(0, 0);//dX_G_w	
			M(1, 3) = poro*(P_gp(0, 0)*C_v*X_G_w_gp(0, 0) - rho_mol_L_gp(0, 0));//dSG
			M(1, 4) = poro*S_gp(0, 0)*X_G_w_gp(0, 0)*drho_mol_GdT_gp(0); //dT dT

			M(2, 0) = poro*(X_G_a_gp(0, 0)*M_G + X_G_w_gp(0, 0)*M_L)*C_v*S_gp(0)*h_G_gp(0, 0) - poro*S_gp(0);//dPG
			M(2, 1) = poro*rho_mol_G_gp(0, 0)*(M_G)*S_gp(0)*h_G_gp(0, 0)
				+ poro*rho_mass_G_gp(0, 0)*S_gp(0)*dMassfraction_air(0, 0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap);//dX_G_a

			M(2, 2) = poro*rho_mol_G_gp(0, 0)*(M_L)*S_gp(0)*h_G_gp(0, 0)
				+ poro*rho_mass_G_gp(0, 0)*S_gp(0)*dMassfraction_vap(0, 0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap);//dX_G_w
			//+ poro*rho_mass_G_gp(0, 0)*S_gp(0)*dMassfraction(0, 0)*(C_pg*(T_gp(0, 0) - 273.15) - C_pl*(T_gp(0, 0) - 273.15) - delta_h_vap);
			M(2, 3) = poro*(rho_mass_G_gp(0, 0)*h_G_gp(0, 0) - rho_mass_L_gp(0, 0)*C_pl*(T_gp(0, 0) - 273.15))
				- poro*P_gp(0, 0);
			M(2, 4) = (1 - poro)*rho_S_std*C_solid + poro*rho_mass_G_gp(0, 0)*(C_pg*massfractionX_gp(0, 0) + C_pl*(1 - massfractionX_gp(0, 0)))*S_gp(0, 0)
				+ poro*rho_mass_L_gp(0, 0)*C_pl*(1 - S_gp(0, 0))
				+ poro*(X_G_a_gp(0, 0)*M_G + X_G_w_gp(0, 0)*M_L)*S_gp(0)*h_G_gp(0, 0)*drho_mol_GdT_gp(0);

			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 5; jj++){
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
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0, 0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0, 0));
			lambda_L = K_perm*Kr_L_gp / mu_L;
			lambda_G = K_perm*Kr_G_gp / mu_G;
			//+++++++++++++++++++++++++Calculate the vecocity+++++++++++++++++++++++++++++++++++++++
			vel_G_gp = -lambda_G*(dNp)*u1.head(n_nodes);

			vel_L_gp = -((lambda_L*(dNp)*u1.head(n_nodes)) - lambda_L*dPC_dSg_gp(0, 0)*(dNp)*u1.block(3 * n_nodes, 0, n_nodes, 1));// +lambda_L*(dNp)*PC_node;

			vel_X_air_gp = -(dNp)*u1.block(n_nodes, 0, n_nodes, 1);
			vel_X_vap_gp = -(dNp)*u1.block(2 * n_nodes, 0, n_nodes, 1);
			//vel_L_gp = -(lambda_L*(dNp)*u1.head(n_nodes) );//			

			//_function_data->get_vel_G()->setIntegrationPointValue(ele_id, j, vel_G_gp);
			//_function_data->get_vel_L()->setIntegrationPointValue(ele_id, j, vel_L_gp);
			//+++++++++++++++++++++++++ End Calculate+++++++++++++++++++++++++++++++++++++++
			A1(2, 4) = (C_pl*rho_l_std*vel_L_gp(0) + (C_pg*massfractionX_gp(0, 0) + C_pl*(1 - massfractionX_gp(0, 0)))*rho_mass_G_gp(0, 0)*vel_G_gp(0))
				+ C_pg*M_G * poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0)*vel_X_air_gp(0)
				+ C_pl*M_L* poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0)*vel_X_vap_gp(0);// 
			if (n_dim > 1){
				A2(2, 4) = (C_pl*rho_l_std*vel_L_gp(1) + (C_pg*massfractionX_gp(0, 0) + C_pl*(1 - massfractionX_gp(0, 0)))*rho_mass_G_gp(0, 0)*vel_G_gp(1))
					+ C_pg*M_G * poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0)*vel_X_air_gp(1)
					+ C_pl*M_L* poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0)*vel_X_vap_gp(1);// 
			}
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 5; jj++){
					tmp_adv(0, 0) = A1(ii, jj);
					if (n_dim > 1){
						tmp_adv(0, 1) = A2(ii, jj);
					}
					localAdvection_tmp.setZero();
					fe->integrateWxDN(j, tmp_adv, localAdvection_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localAdvection_tmp;
				}
			}
			//+++++++++++++++++++++++

			//Calc each entry of the Laplace Matrix
			D(0, 0) = X_G_a_gp(0, 0)*rho_mol_G_gp(0, 0)*lambda_G;
			D(0, 1) = poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0);//dX_G_a
			D(0, 2) = 0.0;//
			D(0, 3) = 0.0;
			D(0, 4) = 0.0;

			D(1, 0) = lambda_L*rho_mol_L_gp(0, 0) + lambda_G*rho_mol_G_gp(0, 0)*X_G_w_gp(0, 0);
			D(1, 1) = 0.0;//
			D(1, 2) = poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0);
			D(1, 3) = -lambda_L*rho_mol_L_gp(0, 0)*dPC_dSg_gp(0, 0);
			D(1, 4) = 0.0;

			D(2, 0) = lambda_G*rho_mass_G_gp(0, 0)*h_G_gp(0, 0) + lambda_L*rho_mass_L_gp(0, 0)*C_pl*(T_gp(0, 0) - 273.15);
			D(2, 1) = (C_pg*(T_gp(0, 0) - 273.15)*M_G)*poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0);//
			D(2, 2) = (C_pl*(T_gp(0, 0) - 273.15) + delta_h_vap)*M_L*poro*S_gp(0, 0)*D_G*rho_mol_G_gp(0, 0);
			D(2, 3) = -lambda_L*rho_mass_L_gp(0, 0)*C_pl*(T_gp(0, 0) - 273.15)*dPC_dSg_gp(0, 0);
			D(2, 4) = Lam_pm_gp(0, 0);

			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 5; jj++){
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
			//H(0) = -rho_L_h_gp(0, 0)*lambda_L*(rho_L_h_gp(0, 0) + rho_L_std) - rho_G_h_gp(0, 0)*rho_G_h_gp(0,0)*lambda_G; //-pow(RHO_G, 2)*lambda_G - pow(RHO_L, 2)*lambda_L;
			//H(1) = -rho_l_std*lambda_L*(rho_L_h_gp(0, 0) + rho_L_std);
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
			//----------------end assembly--------------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << localF << std::endl;

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
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];
	//LocalMatrixType* _M;
	//LocalMatrixType* _K;
	MathLib::LocalVector local_r_1;
	MathLib::LocalVector local_r_2;
	MathLib::LocalMatrix M1, M2, K1, K2, F1, F2, TMP_M_1, TMP_M_2;

	MathLib::LocalVector vec_g; // gravity term 
	MathLib::LocalVector H;//this is for local gravity vector 
	MathLib::LocalVector localGravity_tmp;
	MathLib::LocalVector localRes;

	MathLib::LocalVector PC;
	MathLib::LocalVector rho_L_h;
	MathLib::LocalVector rho_G_h;
	MathLib::LocalVector dPC_dSg;
	MathLib::LocalVector dMassfraction_air;
	MathLib::LocalVector dMassfraction_vap;
	MathLib::LocalVector Input;
	MathLib::LocalMatrix Xvap_gp;

	MathLib::LocalMatrix  P_G_w_gp;
	MathLib::LocalMatrix  P_G_a_gp;

	MathLib::LocalMatrix  rho_mol_G_gp;
	MathLib::LocalMatrix  rho_mol_L_gp;
	MathLib::LocalMatrix  rho_mass_G_gp;
	MathLib::LocalMatrix  rho_mass_L_gp;

	MathLib::LocalMatrix  PC_gp;
	MathLib::LocalMatrix  dPC_dSg_gp;
	MathLib::LocalMatrix  dP_gw_dT_gp;
	MathLib::LocalMatrix  dP_gw_dPC_gp;
	MathLib::LocalMatrix  massfractionX_gp;
	MathLib::LocalMatrix h_G_gp;//specific enthalpy of gas phase

	LocalVectorType drho_mol_GdT_gp;
	LocalVectorType drho_mol_LdT_gp;
	LocalVectorType drho_G_adX_gp;
	LocalVectorType drho_G_adS_gp;
	LocalVectorType drho_G_adT_gp;

	LocalVectorType Lam_pm_gp;
	LocalVectorType vel_L_gp;
	LocalVectorType vel_G_gp;
	LocalVectorType vel_X_air_gp;
	LocalVectorType vel_X_vap_gp;

	LocalVectorType PC_node;
	//LocalMatrixType PG_gp;
	//LocalMatrixType PL_gp;

	//primary variable on gauss point
	MathLib::LocalMatrix  P_gp;
	MathLib::LocalMatrix  X_G_a_gp;
	MathLib::LocalMatrix  X_G_w_gp;
	MathLib::LocalMatrix  S_gp;
	MathLib::LocalMatrix  T_gp;

	MathLib::LocalMatrix  M;//Local Mass Matrix
	MathLib::LocalMatrix  D;//Local Laplace Matrix
	MathLib::LocalMatrix  A1;//local advection matrix
	MathLib::LocalMatrix  A2;
	MathLib::LocalMatrix  tmp;//method 2
	MathLib::LocalMatrix  tmp_adv;

	MathLib::LocalMatrix  localMass_tmp; //for method2
	MathLib::LocalMatrix  localDispersion_tmp; //for method2
	MathLib::LocalMatrix  localDispersion;
	MathLib::LocalMatrix  localAdvection_tmp; //for method2
	MathLib::LocalMatrix  matN;
	MathLib::LocalMatrix  LocalJacb;


	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w

	const double R = 8.314;
	double C_v;
	double C_w;//M_L/RT
	double C_pg, C_pl;
	double lambda_n;
	double lambda_w;
	const double rho_l_std = 1000.0;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	const double delta_h_vap =  2258000;//J/kg
	//double delta_h_vap;//J/kg
	double rho_S_std;
	double C_solid;
	EOS_NonIsoGlobalNCPForm * _EOS;
	std::size_t isinf;
	double flag;
};

#endif  // end of ifndef
