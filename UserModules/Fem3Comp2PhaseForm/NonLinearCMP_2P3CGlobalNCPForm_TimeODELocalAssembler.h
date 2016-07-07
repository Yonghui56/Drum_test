/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2015-06-19 by Yonghui HUANG
*/

#ifndef NON_LINEAR_CMP_2P3CGLOBALNCPFORM_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_2P3CGLOBALNCPFORM_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"
#include "EOS_2P3CGlobalNCPForm.h"




/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;

	NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
	{
		std::vector<double> vec_x = { 0, 418.848, 837.696, 1256.54, 2094.24, 4188.48, 6701.57, 9633.51, 12146.6, 14240.8, 17172.8,
			20523.6, 23455.5, 26387.4, 30994.8, 34764.4, 39371.7, 43141.4, 46073.3, 49842.9, 54450.3, 64921.5, 74973.8, 79581.2 };
		std::vector<double> vec_y = { 0.0019, 0.00156027, 0.001435395, 0.001078989, 0.000846754, 0.000548841, 0.000387103
			, 0.000268217, 0.000141412, 4.01755e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0 };
		std::vector<double> vec_xx = { 0, 86.9565, 152.174, 195.652, 304.348, 413.043, 543.478, 717.391, 891.304, 1108.7, 1282.61
			, 1543.48, 2217.39, 2608.7, 3021.74, 3413.04, 3782.61, 4000 };
		std::vector<double> vec_yy = { 0.035, 0.02807498, 0.02331248, 0.020707208, 0.018903492, 0.01587534, 0.014197904, 0.010952603
			, 0.00747816, 0.006338911, 0.004480776, 0.001839415, 0.000335932, 0, 0, 0, 0, 0 };
		xy_curve = new //slow degradation curve
			MathLib::LinearInterpolation(vec_x, vec_y);
		curve_slow = new NumLib::FunctionLinear1D(xy_curve);
		xy_curve2 = new //slow degradation curve
			MathLib::LinearInterpolation(vec_xx, vec_yy);
	    curve_fast = new NumLib::FunctionLinear1D(xy_curve2);
		_EOS = new EOS_2P3CGlobalNCPForm();
	};
	
	virtual ~NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler() {
		_function_data = NULL;
		BaseLib::releaseObject(_EOS);
		BaseLib::releaseObject(curve_slow);
		BaseLib::releaseObject(curve_fast);
		BaseLib::releaseObject(xy_curve);
		BaseLib::releaseObject(xy_curve2);
	};
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}

protected:
	
	virtual void assembleODE(const NumLib::TimeStep & time, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
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
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
		const std::size_t n_gsp = q->getNumberOfSamplingPoints(); // number of Gauss points
		const std::size_t n_dof = u1.size();
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

		const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
		double geo_area = 1.0;
		pm->geo_area->eval(e_pos, geo_area);
		const bool hasGravityEffect =  _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES initialize
		*/
		
		PC_gp = LocalMatrixType::Zero(1, 1);
		PL_gp = LocalVectorType::Zero(1, 1);
		rho_mol_G_H_gp = LocalVectorType::Zero(1, 1);
		rho_mol_G_W_gp = LocalVectorType::Zero(1, 1);
		rho_mol_G_CH4_gp = LocalVectorType::Zero(1, 1);
		rho_mol_L_w_gp = LocalVectorType::Zero(1, 1);

		dPC_dSg_gp = LocalMatrixType::Zero(1, 1);
		X_L_h_gp = LocalMatrixType::Zero(1, 1);
		X_L_w_gp = LocalMatrixType::Zero(1, 1);
		X_L_c_gp = LocalMatrixType::Zero(1, 1);
		X_L_co2_gp = LocalMatrixType::Zero(1, 1);
		X_L_air_gp = LocalMatrixType::Zero(1, 1);

		rho_mol_G_gp = LocalVectorType::Zero(1, 1);
		rho_mol_L_gp = LocalVectorType::Zero(1, 1);
		rho_mass_G_gp = LocalVectorType::Zero(1, 1);
		rho_mass_L_gp = LocalVectorType::Zero(1, 1);

		P_sat_gp = LocalVectorType::Zero(1, 1);
		// Loop over each gauss point on this element e
		// calculate the secondary variables' values on each gauss point 
		// calculate the derivative of secondary variables based on P and X on each gauss point
		// store the value in the vector respectively

		M = LocalMatrixType::Zero(n_dof / n_nodes-2, n_dof / n_nodes);//method 2
		D = LocalMatrixType::Zero(n_dof / n_nodes - 2, n_dof / n_nodes);//method 2
		H = LocalVectorType::Zero(n_dof / n_nodes - 2);//for gravity term
		F = LocalVectorType::Zero(n_dof / n_nodes - 2);//for source term
		tmp = LocalMatrixType::Zero(1, 1);
		Input = LocalVectorType::Zero(7);
		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(n_nodes, n_nodes);
		//LocalJacb = LocalMatrixType::Zero(3*n_nodes, 3 * n_nodes);
		localGravity_tmp = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		localSource_tmp = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		localRes = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		double Q_concrete = 2.57635;//consume water
		double Q_steel = 0.28;//generate H2
		double Q_organic_fast_co2 = 0.035;//
		double Q_organic_fast_ch4 = 0.035;
		double Q_organic_slow_co2 = 0.0019;
		double Q_organic_slow_ch4 = 0.0019;
		
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

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			// evaluation of the shape function
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();
			F.setZero();
			P_gp = Np*u1.head(n_nodes);
			X1_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);//H2
			X2_gp = Np*u1.block(2 * n_nodes, 0, n_nodes, 1);//ch4
			X3_gp = Np*u1.block(3 * n_nodes, 0, n_nodes, 1);//h2o
			X4_gp = Np*u1.block(4 * n_nodes, 0, n_nodes, 1);//air
			X5_gp = Np*u1.block(5 * n_nodes, 0, n_nodes, 1);//co2
			S_gp = Np*u1.tail(n_nodes);
			
			P_sat_gp(0, 0) = _EOS->get_P_sat(T_0);
			//X_L_w_gp(0, 0) = P_gp(0, 0)*X3_gp(0, 0) / P_sat_gp(0, 0);
			X_L_h_gp(0, 0) = P_gp(0, 0)*X1_gp(0, 0) / Hen_L_h;
			X_L_c_gp(0, 0) = P_gp(0, 0)*X2_gp(0, 0) / Hen_L_c;
			X_L_co2_gp(0, 0) = P_gp(0, 0)*X5_gp(0, 0) / Hen_L_co2;
			X_L_air_gp(0, 0) = P_gp(0, 0)*X4_gp(0, 0) / Hen_L_air;
			X_L_w_gp(0, 0) = P_gp(0, 0)*X3_gp(0, 0) / P_sat_gp(0, 0);
			//S_gp(0, 0) = 0.3;
			PC_gp(0,0) = pm->getPc_bySat(S_gp(0, 0));
			PL_gp(0, 0) = P_gp(0, 0) - PC_gp(0, 0);
			dPC_dSg_gp(0,0) = pm->Deriv_dPCdS(S_gp(0, 0));
			//store sacondary variables value
			_function_data->get_PL()->setIntegrationPointValue(ele_id, j, PL_gp);
			

			rho_mol_G_gp(0, 0) = P_gp(0, 0) / R / T_0;
			rho_mol_L_w_gp(0, 0) = rho_L_std / M_L;// / (P_gp(0, 0)*X3_gp(0, 0) / P_sat_gp(0, 0));//N_
			rho_mol_L_gp(0, 0) = rho_mol_L_w_gp(0, 0) / X_L_w_gp(0, 0);
			///rho_mol_G_H_gp(0, 0) = rho_mol_L_gp(0, 0)*X_L_h_gp(0, 0);
			///rho_mol_G_CH4_gp(0, 0) = rho_mol_L_gp(0, 0)*X_L_c_gp(0, 0);
			_function_data->get_rhoGh()->setIntegrationPointValue(ele_id, j, rho_mol_G_H_gp);
			_function_data->get_rhoGch4()->setIntegrationPointValue(ele_id, j, rho_mol_G_CH4_gp);
			_function_data->get_rhoGw()->setIntegrationPointValue(ele_id, j, rho_mol_G_W_gp);


			rho_mass_G_gp(0, 0) = rho_mol_G_gp(0, 0)*(X1_gp(0, 0)*M_H + X2_gp(0, 0)*M_C + X3_gp(0, 0)*M_L +
				X4_gp(0, 0)*M_AIR + X5_gp(0, 0)*M_CO2);
			rho_mass_L_gp(0, 0) = rho_L_std;// rho_mol_L_gp(0, 0)*(X_L_h_gp(0, 0)*M_H + X_L_c_gp(0, 0)*M_C + X_L_w_gp(0, 0)*M_L);

			isinf = std::isfinite(dPC_dSg_gp(0, 0));
			if (isinf == 0)
			{
				dPC_dSg_gp(0, 0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0, 0) = dPC_dSg_gp(0, 0);
			}

			
			
			M(0, 0) = poro*(S_gp(0, 0)*X1_gp(0, 0) / R / T_0 + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*X1_gp(0, 0) / Hen_L_h);//dPG
			M(0, 1) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_h);//dX1 h2
			M(0, 2) = 0.0;//Dx2 ch4
			M(0, 3) = -poro*(1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0)*X1_gp(0, 0) / Hen_L_h/X3_gp(0,0);//dx3 h2o
			M(0, 4) = 0.0;//dx4 air
			M(0, 5) = 0.0;//dx5 co2
			M(0, 6) = poro*(rho_mol_G_gp(0, 0)*X1_gp(0, 0) - rho_mol_L_gp(0, 0)*P_gp(0, 0)*X1_gp(0, 0) / Hen_L_h);//dSat

			M(1, 0) = poro*(S_gp(0, 0)*X2_gp(0, 0) / R / T_0 + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / Hen_L_c);
			M(1, 1) = 0.0;
			M(1, 2) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_c);
			M(1, 3) = -poro*(1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0)*X2_gp(0, 0) / Hen_L_c / X3_gp(0, 0);//dx3 h2o
			M(1, 4) = 0.0;
			M(1, 5) = 0.0;
			M(1, 6) = poro*(rho_mol_G_gp(0, 0)*X2_gp(0, 0) - rho_mol_L_gp(0, 0)*P_gp(0, 0)*X2_gp(0, 0) / Hen_L_c);

			//M(2, 0) = poro*(S_gp(0, 0)*X3_gp(0, 0) / R / T_0 + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*X3_gp(0, 0) / P_sat_gp(0, 0));//dPG poro*(N_G*S_G*x_g_w+N_L*S_L*X_L_W)//dpg
			M(2, 0) = poro*(S_gp(0, 0)*X3_gp(0, 0) / R / T_0);//dPG poro*(N_G*S_G*x_g_w+N_L*S_L*X_L_W)//dpg
			M(2, 1) = 0.0;//x1
			M(2, 2) = 0.0;//x2
			M(2, 3) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) );
			M(2, 4) = 0.0;//x4
			M(2, 5) =0.0;//X5
			M(2, 6) = poro*(rho_mol_G_gp(0, 0)*X3_gp(0, 0) - rho_mol_L_w_gp(0, 0));//S

			M(3, 0) = poro*(S_gp(0, 0)*X4_gp(0, 0) / R / T_0 + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*X4_gp(0, 0) / Hen_L_air);//dPG
			M(3, 1) = 0.0;
			M(3, 2) = 0.0;
			M(3, 3) = -poro*(1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0)*X4_gp(0, 0) / Hen_L_air / X3_gp(0, 0);//dx3 h2o
			M(3, 4) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_air);//
			M(3, 5) = 0.0;
			M(3, 6) = poro*(rho_mol_G_gp(0, 0)*X4_gp(0, 0) - rho_mol_L_gp(0, 0)*P_gp(0, 0)*X4_gp(0, 0) / Hen_L_air);

			M(4, 0) = poro*(S_gp(0, 0)*X5_gp(0, 0) / R / T_0 + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*X5_gp(0, 0) / Hen_L_co2);//dPG
			M(4, 1) = 0.0;
			M(4, 2) = 0.0;
			M(4, 3) = -poro*(1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0)*X5_gp(0, 0) / Hen_L_co2 / X3_gp(0, 0);//dx3 h2o
			M(4, 4) = 0.0;//dX1 h2
			M(4, 5) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_co2);
			M(4, 6) = poro*(rho_mol_G_gp(0, 0)*X5_gp(0, 0) - rho_mol_L_gp(0, 0)*P_gp(0, 0)*X5_gp(0, 0) / Hen_L_co2);
			//-------------debugging------------------------
			 //std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 5; ii++){
				for (jj = 0; jj < 7; jj++){
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

			
			//Calc each entry of the Laplace Matrix
			D(0, 0) = lambda_G*rho_mol_G_gp(0, 0)*X1_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X1_gp(0, 0)*P_gp(0, 0) / Hen_L_h;
			D(0, 1) = poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(0, 2) = 0.0;
			D(0, 3) = 0.0;
			D(0, 4) = 0.0;
			D(0, 5) = 0.0;
			D(0, 6) = -lambda_L*rho_mol_L_gp(0, 0)*X1_gp(0, 0)*P_gp(0, 0)*dPC_dSg_gp(0, 0) / Hen_L_h;
			
			D(1, 0) = lambda_G*rho_mol_G_gp(0, 0)*X2_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0)*P_gp(0, 0) / Hen_L_c;
			D(1, 1) = 0.0;
			D(1, 2) = poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(1, 3) = 0.0;
			D(1, 4) = 0.0;
			D(1, 5) = 0.0;
			D(1, 6) = -lambda_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0)*P_gp(0, 0)*dPC_dSg_gp(0, 0) / Hen_L_c;

			D(2, 0) = lambda_G*rho_mol_G_gp(0, 0)*X3_gp(0, 0) + lambda_L*rho_mol_L_w_gp(0, 0);
				//+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / P_sat_gp(0,0);
			D(2, 1) = 0.0;
			D(2, 2) = 0.0;
			D(2, 3) = poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);// poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);// poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / P_sat_gp(0, 0)
			D(2, 4) = 0.0;
			D(2, 5) = 0.0;
			D(2, 6) = -lambda_L*rho_mol_L_w_gp(0, 0)*dPC_dSg_gp(0, 0);

			D(3, 0) = lambda_G*rho_mol_G_gp(0, 0)*X4_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X4_gp(0, 0)*P_gp(0, 0) / Hen_L_air;
			//+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / P_sat_gp(0,0);
			D(3, 1) = 0.0;
			D(3, 2) = 0.0;
			D(3, 3) = 0.0;// poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);// poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / P_sat_gp(0, 0)
			D(3, 4) = poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(3, 5) = 0.0;
			D(3, 6) = -lambda_L*rho_mol_L_gp(0, 0)*dPC_dSg_gp(0, 0)*X4_gp(0, 0)*P_gp(0, 0) / Hen_L_air;

			D(4, 0) = lambda_G*rho_mol_G_gp(0, 0)*X5_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X5_gp(0, 0)*P_gp(0, 0) / Hen_L_co2;
			//+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / P_sat_gp(0,0);
			D(4, 1) = 0.0;
			D(4, 2) = 0.0;
			D(4, 3) = 0.0;// poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);// poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / P_sat_gp(0, 0)
			D(4, 4) = 0.0;
			D(4, 5) = poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(4, 6) = -lambda_L*rho_mol_L_gp(0, 0)*dPC_dSg_gp(0, 0)*X5_gp(0, 0)*P_gp(0, 0) / Hen_L_co2;
			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			// std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 5; ii++){
				for (jj = 0; jj < 7; jj++){
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
			H(0) = -lambda_G*rho_mol_G_gp(0, 0)*X1_gp(0, 0)*rho_mass_G_gp(0, 0) 
				- lambda_L*rho_mol_L_gp(0, 0)*X1_gp(0, 0)*P_gp(0, 0)*rho_mass_L_gp(0, 0) / Hen_L_h;
			H(1) = -lambda_G*rho_mol_G_gp(0, 0)*X2_gp(0, 0)*rho_mass_G_gp(0, 0)
				- lambda_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0)*P_gp(0, 0)*rho_mass_L_gp(0, 0) / Hen_L_c;
			H(2) = -lambda_G*rho_mol_G_gp(0, 0)*X3_gp(0, 0)*rho_mass_G_gp(0, 0)
				- lambda_L*rho_mol_L_gp(0, 0)*X3_gp(0, 0)*P_gp(0, 0)*rho_mass_L_gp(0, 0) / P_sat_gp(0, 0);
			H(3) = -lambda_G*rho_mol_G_gp(0, 0)*X4_gp(0, 0)*rho_mass_G_gp(0, 0)
				- lambda_L*rho_mol_L_gp(0, 0)*X4_gp(0, 0)*P_gp(0, 0)*rho_mass_L_gp(0, 0) / Hen_L_air;
			H(4) = -lambda_G*rho_mol_G_gp(0, 0)*X5_gp(0, 0)*rho_mass_G_gp(0, 0)
				- lambda_L*rho_mol_L_gp(0, 0)*X5_gp(0, 0)*P_gp(0, 0)*rho_mass_L_gp(0, 0) / Hen_L_co2;
			
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			if (hasGravityEffect) {
				// F += dNp^T * H* gz				
				for (int idx = 0; idx < 5; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0,0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;

				}
			}
			
			//----------------end assembly--------------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << localF << std::endl;
			//================load source term===========================
			if (mat_id == 1){
				F(0) += Q_steel;
				F(2) -= Q_steel;
				//curve_slow->eval(time.getTime(), Q_organic_slow_co2);
				//curve_fast->eval(time.getTime(), Q_organic_fast_co2);
				F(4) += Q_organic_slow_co2/ 0.342895;//co2 organic degradation slow
				F(1) += Q_organic_slow_co2*31 / 0.342895 / 19;//ch4 organic degradation slow
				F(2) -= 2 * Q_organic_slow_co2 / 0.129; //h2o consume slow
				F(4) += Q_organic_fast_co2/ 0.342895;//co2 organic degradation fast
				F(1) += Q_organic_fast_co2/ 0.342895;//ch4 organic degradation fast
				F(2) -= Q_organic_fast_co2 / 0.342895 / 3; //h2o consume fast
			}
			else if (mat_id == 0){
				if (X5_gp(0, 0) > 0.01){
					F(4) -= Q_organic_slow_co2 / 0.342895; //carbonation consume co2
					F(2) += Q_organic_slow_co2 / 0.342895;
				}
				F(2) -= 0.007 / 0.342895;//concrete degradation
			}
			for (int idx = 0; idx < 5; idx++){
				tmp(0, 0) = F(idx);
				localSource_tmp.setZero();
				//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
				localSource_tmp = fac * Np.transpose()*tmp(0, 0);
				localF.block(n_nodes*idx, 0, n_nodes, 1) += localSource_tmp;

			}
			//----------------end assembly--------------------------------
			//std::cout << "F=" << std::endl;
			//std::cout << localF << std::endl;

		}
		//----------------end assemby gaus---------------------------------
		
		//std::cout << "H=" << std::endl;
		//std::cout << LocalJacb<< std::endl;
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;
		
		//_function_data->get_elem_J_matrix()[ele_id] = LocalJacb;
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

	LocalVectorType vec_g; // gravity term 
	LocalVectorType H;//this is for local gravity vector 
	LocalVectorType F;// THIS is for source term
	LocalVectorType localGravity_tmp;
	LocalVectorType localRes;
	LocalVectorType localSource_tmp;


	LocalVectorType Input;

	
	LocalMatrixType rho_mol_L_gp;
	LocalMatrixType rho_mol_G_gp;
	LocalMatrixType PC_gp;
	LocalMatrixType dPC_dSg_gp;
	//LocalMatrixType PG_gp;
	LocalVectorType PL_gp;
	LocalVectorType P_sat_gp;
	LocalVectorType rho_mol_G_H_gp;
	LocalVectorType rho_mol_G_W_gp;
	LocalVectorType rho_mol_G_CH4_gp;
	LocalVectorType rho_mol_L_w_gp;

	//primary variable on gauss point
	LocalMatrixType P_gp;
	LocalMatrixType X1_gp;
	LocalMatrixType X2_gp;
	LocalMatrixType X3_gp;
	LocalMatrixType X4_gp;
	LocalMatrixType X5_gp;
	LocalMatrixType S_gp;


	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType tmp;//method 2
	LocalMatrixType X_L_w_gp;
	LocalMatrixType X_L_h_gp;
	LocalMatrixType X_L_c_gp;
	LocalMatrixType X_L_air_gp;
	LocalMatrixType X_L_co2_gp;

	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localDispersion;
	LocalMatrixType matN;
	LocalMatrixType LocalJacb;

	LocalVectorType rho_mass_L_gp;
	LocalVectorType rho_mass_G_gp;

	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;
	const double M_H = 0.002;
	const double M_L = 0.018;
	const double M_C = 0.016;
	const double M_AIR = 0.029;
	const double M_CO2 = 0.044;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w
	
	const double R = 8.314;
	const double Hen_L_h = 7.26e+9; //Henry constant in [Pa]
	const double Hen_L_c = 4.13e+9; //Henry constant in [Pa]
	const double Hen_L_air = 9.077e+9; //Henry constant in [Pa]
	const double Hen_L_co2 = 0.163e+9; //Henry constant in [Pa]
	const double T_0 = 303.15; //in [K]
	const double rho_l_std = 1000.0;
	EOS_2P3CGlobalNCPForm * _EOS;
	std::size_t isinf;
	NumLib::FunctionLinear1D* curve_slow = NULL;
	NumLib::FunctionLinear1D* curve_fast = NULL;
	MathLib::LinearInterpolation* xy_curve = NULL;
	MathLib::LinearInterpolation* xy_curve2 = NULL;
};



#endif  // end of ifndef