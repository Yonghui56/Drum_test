/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.hpp
*
* Created on 2014-05-12 by Yonghui HUANG
*/

#include "logog.hpp"

#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXFunctionDirect.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "MathLib/ODE/RungeKutta4.h"
#include "EOS_NonIso_HeatPipe_CapPressureForm.h"

template <class T1, class T2>
bool FunctionCMP_NonIso_CapPressureForm<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
	NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

	//mesh and FE objects
	MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
	// get the number of elements 
	size_t n_ele = msh->getNumberOfElements();
	size_t n_nodes = msh->getNumberOfNodes();
	MyDiscreteSystem* dis = 0;
	dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);
	//FemLib::IFemNumericalIntegration *integral = _feObjects->getIntegrationMethod();
	// set names of the output parameters
	this->setOutputParameterName(0, "GAS_PRESSURE");
	this->setOutputParameterName(1, "CAPILLARY_PRESSURE");
	this->setOutputParameterName(2, "TEMPERATURE");
	this->setOutputParameterName(3, "Saturation");
	this->setOutputParameterName(4, "Liquid_Pressure");
	this->setOutputParameterName(5, "Mass_Density_L_H");
	this->setOutputParameterName(6, "Mass_Density_G_H");
	this->setOutputParameterName(7, "Mass_Density_G_W");
	this->setOutputParameterName(8, "Velocity_G");
	this->setOutputParameterName(9, "Velocity_L");
	this->setOutputParameterName(10, "Molar_Fraction_G_air");
	this->setOutputParameterName(11, "Molar_Fraction_G_vapor");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMPCapPressureFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPCapPressureForm* mean_pressure = _problem->addVariable("GAS_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("GAS_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, mean_pressure);
	SolutionLib::FemIC* femIC = _problem->getVariable(0)->getIC();
	_PG = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_PG->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
		femIC->setup(*_PG);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_PG->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction

	MyVariableCMPCapPressureForm* total_mass_density = _problem->addVariable("CAPILLARY_PRESSURE");
	var_builder.doit("CAPILLARY_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, total_mass_density);
	femIC = _problem->getVariable(1)->getIC();
	_PC = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_PC->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
		femIC->setup(*_PC);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_PC->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
	}


	MyVariableCMPCapPressureForm* temperature = _problem->addVariable("TEMPERATURE");
	var_builder.doit("TEMPERATURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, temperature);
	femIC = _problem->getVariable(2)->getIC();
	_T = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_T->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
		femIC->setup(*_T);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_T->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
	}

	// initialize the local EOS problem
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);
	_S = new MyIntegrationPointFunctionVector();
	_S->initialize(dis);
	_PL = new MyIntegrationPointFunctionVector();
	_PL->initialize(dis);
	_rho_L_h = new MyIntegrationPointFunctionVector();
	_rho_L_h->initialize(dis);
	_rho_G_h = new MyIntegrationPointFunctionVector();
	_rho_G_h->initialize(dis);
	_rho_G_w = new  MyIntegrationPointFunctionVector();
	_rho_G_w->initialize(dis);
	_vel_G = new MyIntegrationPointFunctionVector();
	_vel_G->initialize(dis);
	_vel_L = new MyIntegrationPointFunctionVector();
	_vel_L->initialize(dis);
	_x_G_air = new MyIntegrationPointFunctionVector();
	_x_G_air->initialize(dis);
	_x_G_vapor = new MyIntegrationPointFunctionVector();
	_x_G_vapor->initialize(dis);
	MathLib::LocalVector tmp = MathLib::LocalVector::Zero(1);
	for (size_t ele_id = 0; ele_id < n_ele; ele_id++)
	{
		
		_S->setNumberOfIntegationPoints(ele_id, 3);	
		_PL->setNumberOfIntegationPoints(ele_id, 3);
		_rho_G_h->setNumberOfIntegationPoints(ele_id, 3);
		_rho_G_w->setNumberOfIntegationPoints(ele_id,3);
		_rho_L_h->setNumberOfIntegationPoints(ele_id, 3);
		_x_G_air->setNumberOfIntegationPoints(ele_id, 3);
		_x_G_vapor->setNumberOfIntegationPoints(ele_id, 3);
		for (size_t jj = 0; jj < 3; jj++)
		{
			_S->setIntegrationPointValue(ele_id, jj, tmp);
			_PL->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_h->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_L_h->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_w->setIntegrationPointValue(ele_id, jj, tmp);
			_x_G_air->setIntegrationPointValue(ele_id, jj, tmp);
			_x_G_vapor->setIntegrationPointValue(ele_id, jj, tmp);
		}
	}

	_vel_G_3D = new My3DIntegrationPointFunctionVector(_vel_G, msh->getGeometricProperty()->getCoordinateSystem());
	_vel_L_3D = new My3DIntegrationPointFunctionVector(_vel_L, msh->getGeometricProperty()->getCoordinateSystem());
	
	_mat_secDer = new MyNodalFunctionMatrix();// define a matrix to store all the derivatives of the secondary variables 
	MathLib::LocalMatrix tmp_mat = MathLib::LocalMatrix::Zero(3, 2); //the size I will modify later
	_mat_secDer->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_mat);
	/**
	* initialize the vector of tempVar
	*/
	
	_vec_tempVar = new MyNodalFunctionVector();
	MathLib::LocalVector tmp_vec = MathLib::LocalVector::Zero(2);
	_vec_tempVar->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_vec);


	//initialize the local M matrix and K matrix
	for (size_t index = 0; index < n_ele; index++)
	{
		MathLib::LocalMatrix test1 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_M_matrix.push_back(test1);
		MathLib::LocalMatrix test2 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_K_matrix.push_back(test2);
	}
	
	// linear assemblers
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMPCapPressureFormProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("GAS_PRESSURE");
	_non_linear_problem->addVariable("CAPILLARY_PRESSURE");
	_non_linear_problem->addVariable("TEMPERATURE");
	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_PG->getDiscreteData()));
	var_builder.doit("GAS_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X_ic = new SolutionLib::FemIC(msh);
	X_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_PC->getDiscreteData()));
	var_builder.doit("CAPILLARY_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(X_ic);

	SolutionLib::FemIC* T_ic = new SolutionLib::FemIC(msh);
	T_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_T->getDiscreteData()));
	var_builder.doit("TEMPERATURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(2));
	_non_linear_problem->getVariable(2)->setIC(T_ic);

	// set up non-linear solution
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory(myNRIterator);
	this->_non_linear_solution = new MyNonLinearSolutionType(dis, this->_non_linear_problem, myNSolverFactory);
	this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
	this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

	// linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // set theta
    if (!optNum->hasOption("TimeTheta"))
    {
        ERR("Time theta setting not found!!!");
        exit(1);
    }
    else
    {
        double tmp_theta(optNum->getOptionAsNum<double>("TimeTheta"));
        non_linear_j_assembler->setTheta(tmp_theta);
        non_linear_r_assembler->setTheta(tmp_theta);
        _non_linear_eqs->getLinearAssembler()->setTheta(tmp_theta);
    }



	// get the nonlinear solution dof manager
	this->_nl_sol_dofManager = this->_non_linear_solution->getDofEquationIdTable();
	// set up solution
	_solution = new MyCMPCapPressureFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _T);
	femData->outController.setOutput(var3.name, var3);
	// 
	OutputVariableInfo var_Sec_1(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	
	OutputVariableInfo var_Sec_3(this->getOutputParameterName(5), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	OutputVariableInfo var_Sec_4(this->getOutputParameterName(6), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_Sec_4.name, var_Sec_4);

	OutputVariableInfo var_Sec_5(this->getOutputParameterName(7), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_w);
	femData->outController.setOutput(var_Sec_5.name, var_Sec_5);

	OutputVariableInfo var_Sec_6(this->getOutputParameterName(8), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vel_G_3D);
	femData->outController.setOutput(var_Sec_6.name, var_Sec_6);

	OutputVariableInfo var_Sec_7(this->getOutputParameterName(9), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vel_L_3D);
	femData->outController.setOutput(var_Sec_7.name, var_Sec_7);

	OutputVariableInfo var_Sec_8(this->getOutputParameterName(10), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _x_G_air);
	femData->outController.setOutput(var_Sec_8.name, var_Sec_8);
	OutputVariableInfo var_Sec_9(this->getOutputParameterName(11), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _x_G_vapor);
	femData->outController.setOutput(var_Sec_9.name, var_Sec_9);
	return true;
}

template <class T1, class T2>
void FunctionCMP_NonIso_CapPressureForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	/*
	size_t i;
	//const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for (i = 0; i < _linear_problems.size(); i++) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
	*/
}

template <class T1, class T2>
void FunctionCMP_NonIso_CapPressureForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_NonIso_CapPressureForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("GAS_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("CAPILLARY_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3("TEMPERATURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _T);
	femData->outController.setOutput(var3.name, var3);
	// and 8 secondary variables
	// add all seconary variables as well. 
	//we have several secondary variables the values of which are defined on each elements
	
	OutputVariableInfo var_Sec_1("Saturation", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Liquid_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	OutputVariableInfo var_Sec_3("Mass_Density_L_H", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	OutputVariableInfo var_Sec_4("Mass_Density_G_H", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_Sec_4.name, var_Sec_4);

	OutputVariableInfo var_Sec_5("Mass_Density_G_W", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_w);
	femData->outController.setOutput(var_Sec_5.name, var_Sec_5);

	OutputVariableInfo var_Sec_10("Molar_Fraction_G_air", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _x_G_air);
	femData->outController.setOutput(var_Sec_10.name, var_Sec_10);
	OutputVariableInfo var_Sec_11("Molar_Fraction_G_vapor", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _x_G_vapor);
	femData->outController.setOutput(var_Sec_11.name, var_Sec_11);
	

	_vel_G_3D->resetVectorFunction(_vel_G);
	OutputVariableInfo var_Sec_6("Velocity_G", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vel_G_3D);
	setOutput(0, _vel_G_3D);
	_vel_L_3D->resetVectorFunction(_vel_L);
	OutputVariableInfo var_Sec_7("Velocity_L", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vel_L_3D);
	setOutput(0, _vel_L_3D);

}


template <class T1, class T2>
void FunctionCMP_NonIso_CapPressureForm<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{

}


