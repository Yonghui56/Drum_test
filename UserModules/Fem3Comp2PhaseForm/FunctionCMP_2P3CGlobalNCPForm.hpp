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
#include "EOS_2P3CGlobalNCPForm.h"

template <class T1, class T2>
bool FunctionCMP_2P3CGlobalNCPForm<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<int>("TimeGroupID");
	NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

	//mesh and FE objects
	MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
	// get the number of elements 
	size_t n_ele = msh->getNumberOfElements();
	size_t n_nodes = msh->getNumberOfNodes();
	MyDiscreteSystem* dis = 0;
	dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// set names of the output parameters
	this->setOutputParameterName(0, "GAS_PRESSURE");
	this->setOutputParameterName(1, "MOLAR_FRACTION_GAS_H");
	this->setOutputParameterName(2, "MOLAR_FRACTION_GAS_C");
	this->setOutputParameterName(3, "MOLAR_FRACTION_GAS_W");
	this->setOutputParameterName(4, "MOLAR_FRACTION_GAS_AIR");
	this->setOutputParameterName(5, "MOLAR_FRACTION_GAS_CO2");
	this->setOutputParameterName(6, "SATURATION");
	this->setOutputParameterName(7, "Liquid_Pressure");
	this->setOutputParameterName(8, "Capillary_Pressure");
	this->setOutputParameterName(9, "Molar_Density_L");
	this->setOutputParameterName(10, "Molar_Density_G");

	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMP2P3CGlobalNCPFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating gas pressure vector
	MyVariableCMP2P3CGlobalNCPForm* Gas_Pressure = _problem->addVariable("GAS_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("GAS_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, Gas_Pressure);
	SolutionLib::FemIC* femIC = _problem->getVariable(0)->getIC();
	_P = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
		femIC->setup(*_P);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction for hydrogen in gas phase

	MyVariableCMP2P3CGlobalNCPForm* X_G_h = _problem->addVariable("MOLAR_FRACTION_GAS_H");
	var_builder.doit("MOLAR_FRACTION_GAS_H", option, msh, femData->geo, femData->geo_unique_name, _feObjects, X_G_h);
	femIC = _problem->getVariable(1)->getIC();
	_X1= new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X1->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
		femIC->setup(*_X1);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X1->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction for CH4 in gas phase

	MyVariableCMP2P3CGlobalNCPForm* X_G_c = _problem->addVariable("MOLAR_FRACTION_GAS_C");
	var_builder.doit("MOLAR_FRACTION_GAS_C", option, msh, femData->geo, femData->geo_unique_name, _feObjects, X_G_c);
	femIC = _problem->getVariable(2)->getIC();
	_X2 = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X2->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
		femIC->setup(*_X2);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X2->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction for vapor in gas phase

	MyVariableCMP2P3CGlobalNCPForm* X_G_w = _problem->addVariable("MOLAR_FRACTION_GAS_W");
	var_builder.doit("MOLAR_FRACTION_GAS_W", option, msh, femData->geo, femData->geo_unique_name, _feObjects, X_G_w);
	femIC = _problem->getVariable(3)->getIC();
	_X3 = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X3->initialize(*dis, _problem->getVariable(3)->getCurrentOrder(), 0.0);
		femIC->setup(*_X3);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X3->initialize(*dis, _problem->getVariable(3)->getCurrentOrder(), 0.0);
	}
	// creating molar fraction for AIR in gas phase
	MyVariableCMP2P3CGlobalNCPForm* X_G_air = _problem->addVariable("MOLAR_FRACTION_GAS_AIR");
	var_builder.doit("MOLAR_FRACTION_GAS_AIR", option, msh, femData->geo, femData->geo_unique_name, _feObjects, X_G_air);
	femIC = _problem->getVariable(4)->getIC();
	_X4 = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X4->initialize(*dis, _problem->getVariable(4)->getCurrentOrder(), 0.0);
		femIC->setup(*_X4);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X4->initialize(*dis, _problem->getVariable(4)->getCurrentOrder(), 0.0);
	}
	//CREATING molar fraction for CO2 in gas phase
	// creating molar fraction for AIR in gas phase
	MyVariableCMP2P3CGlobalNCPForm* X_G_co2 = _problem->addVariable("MOLAR_FRACTION_GAS_CO2");
	var_builder.doit("MOLAR_FRACTION_GAS_CO2", option, msh, femData->geo, femData->geo_unique_name, _feObjects, X_G_co2);
	femIC = _problem->getVariable(5)->getIC();
	_X5 = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X5->initialize(*dis, _problem->getVariable(5)->getCurrentOrder(), 0.0);
		femIC->setup(*_X5);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X5->initialize(*dis, _problem->getVariable(5)->getCurrentOrder(), 0.0);
	}
	// creating primary variable of molar fraction of water vapor in 
	MyVariableCMP2P3CGlobalNCPForm* saturation = _problem->addVariable("SATURATION");
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, saturation);
	femIC = _problem->getVariable(6)->getIC();
	_S = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_S->initialize(*dis, _problem->getVariable(6)->getCurrentOrder(), 0.0);
		femIC->setup(*_S);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_S->initialize(*dis, _problem->getVariable(6)->getCurrentOrder(), 0.0);
	}

	// initialize the local EOS problem
	
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);

	_PC = new MyNodalFunctionScalar();	
	_PC->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_rho_L_h = new MyNodalFunctionScalar();//MASS density of hydrogen in liquid phase
	_rho_L_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	//_rho_G_h = new MyNodalFunctionScalar();//mass density of hydrogen in gas phase
	//_rho_G_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_dPcdSg = new MyNodalFunctionScalar();
	_dPcdSg->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_X_L_c = new MyNodalFunctionScalar();//molar fraction of hydrogen in gas phase
	_X_L_c->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_X_L_h = new MyNodalFunctionScalar(); // molar fraction of hydrogen in gas phase
	_X_L_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_X_L_w = new MyNodalFunctionScalar();//molar fraction of hydrogen in gas phase
	_X_L_w->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_PL = new MyIntegrationPointFunctionVector();
	_rho_G_h = new MyIntegrationPointFunctionVector();
	_rho_G_w = new MyIntegrationPointFunctionVector();
	_rho_G_ch4 = new MyIntegrationPointFunctionVector();
	_PL->initialize(dis);
	_rho_G_h->initialize(dis);
	_rho_G_w->initialize(dis);
	_rho_G_ch4->initialize(dis);
	MathLib::LocalVector tmp = MathLib::LocalVector::Zero(1);
	for (size_t ele_id = 0; ele_id < n_ele; ele_id++)
	{
		MeshLib::IElement *e = msh->getElement(ele_id);
		const std::size_t node_ele = e->getNumberOfNodes();//
		_PL->setNumberOfIntegationPoints(ele_id, node_ele);
		_rho_G_h->setNumberOfIntegationPoints(ele_id, node_ele);
		_rho_G_w->setNumberOfIntegationPoints(ele_id, node_ele);
		_rho_G_ch4->setNumberOfIntegationPoints(ele_id, node_ele);
		for (size_t jj = 0; jj < node_ele; jj++)
		{
			_PL->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_h->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_w->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_ch4->setIntegrationPointValue(ele_id, jj, tmp);
		}
	}
	/**
	* initialize the vector of tempVar
	*/
	_mat_Jacob.resize(
	    7 * n_nodes,
	    7 * n_nodes);  
	// Resizing also initializes the matrix with zero
	_vec_Res = LocalVector::Zero(7 * n_nodes);

	
	//initialize the local M matrix and K matrix and partially Jacobian matrix
	for (size_t index = 0; index < n_ele; index++)
	{
		MathLib::LocalMatrix test1 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_M_matrix.push_back(test1);
		MathLib::LocalMatrix test2 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_K_matrix.push_back(test2);
		MathLib::LocalMatrix test3 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_J_matrix.push_back(test3);
	}
	

	// linear assemblers
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMP2P3CGlobalNCPFormProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	//_non_linear_eqs->initialize(nonlin_assembler,non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("GAS_PRESSIRE");
	_non_linear_problem->addVariable("MOLAR_FRACTION_GAS_H");
	_non_linear_problem->addVariable("MOLAR_FRACTION_GAS_C");
	_non_linear_problem->addVariable("MOLAR_FRACTION_GAS_W");
	_non_linear_problem->addVariable("MOLAR_FRACTION_GAS_AIR");
	_non_linear_problem->addVariable("MOLAR_FRACTION_GAS_CO2");
	_non_linear_problem->addVariable("SATURATION");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("GAS_PRESSIRE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X1_ic = new SolutionLib::FemIC(msh);
	X1_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X1->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION_GAS_H", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(X1_ic);

	SolutionLib::FemIC* X2_ic = new SolutionLib::FemIC(msh);
	X2_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X2->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION_GAS_C", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(2));
	_non_linear_problem->getVariable(2)->setIC(X2_ic);

	SolutionLib::FemIC* X3_ic = new SolutionLib::FemIC(msh);
	X3_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X3->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION_GAS_W", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(3));
	_non_linear_problem->getVariable(3)->setIC(X3_ic);

	SolutionLib::FemIC* X4_ic = new SolutionLib::FemIC(msh);
	X4_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X4->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION_GAS_AIR", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(4));
	_non_linear_problem->getVariable(4)->setIC(X4_ic);

	SolutionLib::FemIC* X5_ic = new SolutionLib::FemIC(msh);
	X5_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X5->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION_GAS_CO2", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(5));
	_non_linear_problem->getVariable(5)->setIC(X5_ic);

	SolutionLib::FemIC* S_ic = new SolutionLib::FemIC(msh);
	S_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_S->getDiscreteData()));
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(6));
	_non_linear_problem->getVariable(6)->setIC(S_ic);

	// set up non-linear solution
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory(myNRIterator);
	this->_non_linear_solution = new MyNonLinearSolutionType(dis, this->_non_linear_problem, this,myNSolverFactory);
	this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
	this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

	// linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // set theta----------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	_solution = new MyCMP2P3CGlobalNCPFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X1);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X2);
	femData->outController.setOutput(var3.name, var3);
	OutputVariableInfo var4(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X3);
	femData->outController.setOutput(var4.name, var4);
	OutputVariableInfo var5(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X4);
	femData->outController.setOutput(var5.name, var5);
	OutputVariableInfo var6(this->getOutputParameterName(5), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X5);
	femData->outController.setOutput(var6.name, var6);
	OutputVariableInfo var7(this->getOutputParameterName(6), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var7.name, var7);

	OutputVariableInfo var_sec_1(this->getOutputParameterName(7), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_sec_1.name, var_sec_1);

	return true;
}

template <class T1, class T2>
void FunctionCMP_2P3CGlobalNCPForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
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
void FunctionCMP_2P3CGlobalNCPForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_2P3CGlobalNCPForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("GAS_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("MOLAR_FRACTION_GAS_H", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X1);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3("MOLAR_FRACTION_GAS_C", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X2);
	femData->outController.setOutput(var3.name, var3);
	OutputVariableInfo var4("MOLAR_FRACTION_GAS_W", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X3);
	femData->outController.setOutput(var4.name, var4);
	OutputVariableInfo var5("MOLAR_FRACTION_GAS_AIR", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X4);
	femData->outController.setOutput(var5.name, var5);
	OutputVariableInfo var6("MOLAR_FRACTION_GAS_CO2", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X5);
	femData->outController.setOutput(var6.name, var6);
	OutputVariableInfo var7("SATURATION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var7.name, var7);
	// and 8 secondary variables
	// add all seconary variables as well. 
	OutputVariableInfo var_Sec_1("Liquid_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Capillary_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);


}


template <class T1, class T2>
void FunctionCMP_2P3CGlobalNCPForm<T1, T2>::calc_nodal_eos_sys(double dt)
{
	// solve the EOS system on each node
	std::size_t node_id(0);
	std::size_t num_nodes(0);
	num_nodes = _P->getDiscreteData()->getRangeEnd();
	double EPS = 1E-7;
	double _P_sat(0);
	double _X_L_h_np(0.0);
	double _X_L_c_np(0.0);
	double _X_L_air_np(0.0);
	double _X_L_co2_np(0.0);
	double _X_L_w_np(0.0);
	MathLib::LocalMatrix test = MathLib::LocalMatrix::Zero(1, 1);

	ogsChem::LocalVector input = ogsChem::LocalVector::Zero(7);
	
	//_mat_Jacob.setZero();
	INFO("Calculating EOS on each node...");

	// loop for each node
	for (node_id = _P->getDiscreteData()->getRangeBegin();
		node_id < _P->getDiscreteData()->getRangeEnd();
		node_id++)
	{
		input(0) = _P->getValue(node_id);
		input(1) = _X1->getValue(node_id);
		input(2) = _X2->getValue(node_id);
		input(3) = _X3->getValue(node_id);
		input(4) = _X4->getValue(node_id);
		input(5) = _X5->getValue(node_id);//SG gas phase
		input(6) = _S->getValue(node_id);//SG gas phase
		//This part is for standard newton iteration of local problem
		//solve EOS 
		//This part is for complementary condition
		_EOS->set_env_condition(input);//apply the primary variables as the input variables for Eos
		_P_sat= _EOS->get_P_sat(T_0);
		_X_L_h_np = input(0)*input(1) / Hen_L_h;
		_X_L_c_np = input(0)*input(2) / Hen_L_c;
		
		_X_L_air_np = input(0)*input(4) / Hen_L_air;
		_X_L_co2_np = input(0)*input(5) / Hen_L_co2;
		
		_X_L_w_np = input(0)*input(3) / _P_sat;
		double  check= 1 - (_X_L_h_np + _X_L_c_np  + _X_L_w_np + _X_L_air_np + _X_L_co2_np);
		//For minimum function 1
		//_vec_Res(node_id ) = std::min( input(2), C_h*(input(0) + output(2)) - rho_L_std*M_G*input(1) / M_L );//+ 2 * num_nodes
		_vec_Res(5 * num_nodes + node_id) = std::min(input(6), 1 - (input(1) + input(2) + input(3) + input(4) + input(5)));
		if (input(6) <= 1 - (input(1) + input(2) + input(3) + input(4)+input(5))){//then rho_L^h=C_h*PG
			//Calc each entry of the mass matrix
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 1) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 2) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 3) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 4) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 5) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 6) = 1.0;
		}
		else {//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id) = 0.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 1) = -1.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 2) = -1.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 3) = -1.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 4) = -1.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 5) = -1.0;
			_mat_Jacob.coeffRef(5 * num_nodes + node_id, 7 * node_id + 6) = 0.0;
		}
		//For minimum function 2

		_vec_Res(6 * num_nodes + node_id) = std::min(1 - input(6), 1 - (_X_L_h_np + _X_L_c_np + _X_L_w_np + _X_L_air_np + _X_L_co2_np));
		if (1 - input(6) <= 1 - (_X_L_h_np + _X_L_c_np + _X_L_w_np + _X_L_air_np + _X_L_co2_np)){//then rho_L^h=C_h*PG
			//Calc each entry of the mass matrix
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 1) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 2) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 3) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 4) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 5) = 0.0;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 6) = -1.0;
		}
		else {//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id ) = 
				-(input(1) / Hen_L_h + input(2) / Hen_L_c+input(3) / _P_sat + input(4) / Hen_L_air + input(5) / Hen_L_co2);
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 1) = -input(0) / Hen_L_h;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 2) = -input(0) / Hen_L_c;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 3) = -input(0) / _P_sat;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 4) = -input(0) / Hen_L_air;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 5) = -input(0) / Hen_L_co2;
			_mat_Jacob.coeffRef(6 * num_nodes + node_id, 7 * node_id + 6) = 0.0;
		}




		// end
	}

		
	

}
