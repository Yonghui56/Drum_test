/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.h
*
* Created on    2014-05-12 by Yonghui HUANG
*/

#ifndef FUNCTIONCMP_NONISO_CAPPRESSUREFORM_H
#define FUNCTIONCMP_NONISO_CAPPRESSUREFORM_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "MathLib/DataType.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "UserModules/FemCompNonIsoCapPressureForm/NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler.h"
#include "UserModules/FemCompNonIsoCapPressureForm/NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler.h"
#include "UserModules/FemCompNonIsoCapPressureForm/NonIso_CapPressureForm_Nested_EOS_NRIterationStepInitializer.h"
#include "Fem_CMP_NonIso_CapPressure_Solution.h"
#include "SingleStepCMP_NonIso_CapPressureForm.h"

#include "EOS_NonIso_HeatPipe_CapPressureForm.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionCMP_NonIso_CapPressureForm
	: public ProcessLib::Process
{
public:
	// input variable is velocity
	enum In {};//Velocity = 0 
	// no output variable
	enum Out { Mean_Pressure = 0, Total_Mass_Density = 0 };

	enum EOS_EVL_METHOD { 
		EOS_EVL_FIN_DIFF, 
		EOS_EVL_ANALYTICAL
	};

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionCMP_NonIso_CapPressureForm       MyFunctionData;    // FEM problem class
	typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
	typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver

	// memory for discretized concentration vector
	typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename FemLib::FemNodalFunctionMatrix<MyDiscreteSystem>::type MyNodalFunctionMatrix; 
	typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
	typedef typename NumLib::TXWrapped3DVectorFunction<MyIntegrationPointFunctionVector> My3DIntegrationPointFunctionVector;

	// local assembler
	typedef NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar, FunctionCMP_NonIso_CapPressureForm>      MyNonLinearAssemblerType;
	typedef NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar, FunctionCMP_NonIso_CapPressureForm> MyNonLinearResidualAssemblerType;
	typedef NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType;
	typedef NonIso_CapPressureForm_Nested_EOS_NRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer> MyDiscreteNonlinearSolverFactory;

	/**
	* nonlinear coupled part solving xi
	* Equation definition
	*/
	typedef SolutionLib::TemplateFemEquation<
		MyDiscreteSystem,
		MyLinearSolver,
		MyNonLinearAssemblerType,
		MyNonLinearResidualAssemblerType,
		MyNonLinearJacobianAssemblerType
	> MyNonLinearEquationType;
	/**
	* FEM IVBV problem definition
	*/
	typedef SolutionLib::FemIVBVProblem<
		MyDiscreteSystem,
		MyNonLinearEquationType
	> MyNonLinearCMPCapPressureFormProblemType;

	/**
	* Solution algorithm definition
	*/
	typedef SolutionLib::SingleStepFEM<
		MyNonLinearCMPCapPressureFormProblemType,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;

	/**
	  * the general CompMultiPhase problem part
	  */
	typedef SolutionLib::Fem_CMP_NonIso_CapPressure_Solution<MyDiscreteSystem> MyCMPCapPressureFormProblemType;
	typedef SolutionLib::SingleStepCMP_NonIso_CapPressureForm<MyFunctionData,
		MyCMPCapPressureFormProblemType,
		MyNonLinearCMPCapPressureFormProblemType,
		MyNonLinearSolutionType> MyCMPCapPressureFormSolution;
	typedef typename MyCMPCapPressureFormProblemType::MyVariable MyVariableCMPCapPressureForm;

	/**
	  * constructor
	  */
	FunctionCMP_NonIso_CapPressureForm()
		: Process("CMP_NonIso_CapPressureForm", 0, 3),
		_feObjects(0), _n_Comp(2), _n_Phases(2)
	{
		_EOS = new EOS_NonIso_HeatPipe_CapPressureForm();
		m_EOS_EVL_METHOD = EOS_EVL_ANALYTICAL;
	};

	/**
	  * destructor, reclaim the memory
	  */
	virtual ~FunctionCMP_NonIso_CapPressureForm()
	{
		BaseLib::releaseObject(_PG);
		BaseLib::releaseObject(_PC);

		
		BaseLib::releaseObject(_S);
		BaseLib::releaseObject(_mat_secDer); 
		BaseLib::releaseObject(_EOS);

		/*
		BaseLib::releaseObject(_problem);
		BaseLib::releaseObjectsInStdVector(_concentrations);
		*/
	};

	/**
	  * initialization of the problem class
	  */
	virtual bool initialize(const BaseLib::Options &option);

	/**
	  * finalize but nothing to do here
	  */
	virtual void finalize() {};

	/**
	  * returns the convergence checker
	  */
	virtual NumLib::IConvergenceCheck* getConvergenceChecker()
	{
		return &_checker;
	}

	/**
	  * function to solve the current time step
	  */
	virtual int solveTimeStep(const NumLib::TimeStep &time)
	{
		INFO("Solving %s...", getProcessName().c_str());
		initializeTimeStep(time);
		getSolution()->solveTimeStep(time);
		updateOutputParameter(time);
		return 0;
	}

	/**
	  * function to suggest the next time step
	  */
	virtual double suggestNext(const NumLib::TimeStep &time_current)
	{
		return getSolution()->suggestNext(time_current);
	}

	/**
	  * called when this problem is awake
	  */
	virtual bool isAwake(const NumLib::TimeStep &time)
	{
		return getSolution()->isAwake(time);
	}

	virtual bool accept(const NumLib::TimeStep &time)
	{
		return getSolution()->accept(time);
	}

	/**
	  * called when this time step is accepted
	  */
	virtual void finalizeTimeStep(const NumLib::TimeStep &time)
	{
		output(time);
		getSolution()->finalizeTimeStep(time);
	};

	/**
	  * update the value of primary variable P
	  */
    void set_PG(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
        double nodal_PG;

        for (node_idx = _PG->getDiscreteData()->getRangeBegin();
             node_idx < _PG->getDiscreteData()->getRangeEnd();
			 node_idx++)
		{
            nodal_PG = new_nodal_values->getValue(node_idx);
            this->_PG->setValue(node_idx, nodal_PG);
		}
	};

    /**
      * update the value of primary variable X
      */
    void set_PC(MyNodalFunctionScalar* new_nodal_values)
    {
        std::size_t node_idx;
        double nodal_PC;

        for (node_idx = _PC->getDiscreteData()->getRangeBegin();
             node_idx < _PC->getDiscreteData()->getRangeEnd();
             node_idx++)
        {
			nodal_PC = new_nodal_values->getValue(node_idx);
			this->_PC->setValue(node_idx, nodal_PC);
        }
    };

	/**
	* update the value of primary variable X
	*/
	void set_T(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_T;

		for (node_idx = _T->getDiscreteData()->getRangeBegin();
			node_idx < _T->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_T = new_nodal_values->getValue(node_idx);
			this->_T->setValue(node_idx, nodal_T);
		}
	};

	/*
	MyNodalFunctionScalar* get_concentrations(size_t idx_conc)
	{
		return _concentrations[idx_conc];
	};
	*/

	/**
	  * set mean pressure nodal value
	  */
	void set_PG_node_values(std::size_t node_idx, double node_value){ _PG->setValue(node_idx, node_value);  };

	/**
	  * set molar fraction nodal value
	  */
	void set_PC_node_values(std::size_t node_idx, double node_value){ _PC->setValue(node_idx, node_value); };

	void set_T_node_values(std::size_t node_idx, double node_value){ _T->setValue(node_idx, node_value); };

	/**
	  * calculate nodal Equation-of-State (EOS) system
	  */
	void calc_nodal_eos_sys(double dt);
	/**
	* Here define two functions 
	* return the matrix of M and K
	* would be helpful for constructing the Jacobian matrix
	*/
	std::vector<LocalMatrix>& get_elem_M_matrix(void) { return _elem_M_matrix; };
	std::vector<LocalMatrix>& get_elem_K_matrix(void) { return _elem_K_matrix; };

	/**
	  * return the number of chemical components
	  */
	std::size_t get_n_Comp(void) { return _n_Comp; };


	/**
	* return the value of gas phase saturation
	*/
	MyIntegrationPointFunctionVector* getS(void) { return _S; };
	
	/**
	* return the value of Liquid  pressure
	*/
	MyIntegrationPointFunctionVector* get_PL(void) { return _PL; };
	/**
	* return the dissolved gas mass density
	*/
	MyIntegrationPointFunctionVector* get_rhoLh(void) { return _rho_L_h; };
	/**
	* return the  gas mass density IN GAS PHASE
	*/
	MyIntegrationPointFunctionVector* get_rhoGh(void) { return _rho_G_h; };
	/**
	* return the  water vapor mass density IN GAS PHASE
	*/
	MyIntegrationPointFunctionVector* get_rhoGw(void) { return _rho_G_w; };

	MyIntegrationPointFunctionVector* get_dPGh_dPg(void) { return _dPGh_dPg; };
	MyIntegrationPointFunctionVector* get_dPc_dSg(void) { return _dPcdSg; };

	MyIntegrationPointFunctionVector* get_Func_C(void)  { return _Func_C; };

	MyIntegrationPointFunctionVector* get_vel_L(void)  { return _vel_L; };
	MyIntegrationPointFunctionVector* get_vel_G(void)  { return _vel_G; };

	MyIntegrationPointFunctionVector* get_x_G_air(void)  { return _x_G_air; };
	MyIntegrationPointFunctionVector* get_x_G_vapor(void)  { return _x_G_vapor; };
	
	My3DIntegrationPointFunctionVector* get_vel_3D_L(void) { return _vel_L_3D };
	My3DIntegrationPointFunctionVector* get_vel_3D_G(void) { return _vel_G_3D };
	/**
	* return the matrix of derivative of secondary variable based on P and X
	*/
	MyNodalFunctionMatrix* get_mat_secDer(void) { return _mat_secDer; };
	MyNodalFunctionVector* get_vec_tempVar(void) { return _vec_tempVar; };
protected:
	virtual void initializeTimeStep(const NumLib::TimeStep &time);

	/**
	  * this function is called to exchange output parameters
	  */
	virtual void updateOutputParameter(const NumLib::TimeStep &time);

	/**
	  * get the pointer of solution class for current problem
	  */
	virtual MyCMPCapPressureFormSolution* getSolution() { return _solution; };

	/**
	  * output the result of current solution
	  */
	virtual void output(const NumLib::TimeStep &time);

private:
	DISALLOW_COPY_AND_ASSIGN(FunctionCMP_NonIso_CapPressureForm);

private:
	/**
	* Nonlinear iterator
	*/
	MyNRIterationStepInitializer*              myNRIterator;

	/**
	* Nonlinear solver factory
	*/
	MyDiscreteNonlinearSolverFactory*          myNSolverFactory;

	/**
	* nonlinear equation
	*/
	MyNonLinearEquationType*                  _non_linear_eqs;

	/**
	* non-linear problem
	*/
	MyNonLinearCMPCapPressureFormProblemType*     _non_linear_problem;

	/**
	* non-linear solution
	*/
	MyNonLinearSolutionType*                  _non_linear_solution; 

	/**
	* degree of freedom equation ID talbe for the nonlinear problem
	*/
	DiscreteLib::DofEquationIdTable * _nl_sol_dofManager;

	/**
	  * Component based multiphase flow problem
	  */
	MyCMPCapPressureFormProblemType* _problem;

	/**
	  * solution class for the component based multiphase flow problem
	  */
	MyCMPCapPressureFormSolution*    _solution;

	/**
	  * FEM object
	  */
	FemLib::LagrangeFeObjectContainer* _feObjects;

	/**
	  * convergence checker
	  */
	NumLib::DiscreteDataConvergenceCheck _checker;

	/**
	  * Primary variable 1): vector of mean pressure values on each node
	  */
	MyNodalFunctionScalar* _PG;

	/**
	  * Primary variable 2): vector of molar fraction values of the light component on each node
	  */
	MyNodalFunctionScalar* _PC;
	/**
	* Primary variable 3): vector of Temperature
	*/
	MyNodalFunctionScalar* _T;

	/**
	*		Vector _output
	*		store the eight values of the secondary variables
	*/
	
	/** 
	  * secondary variable --saturation
	  */

	MyIntegrationPointFunctionVector* _S;
	MyIntegrationPointFunctionVector* _PL;
	MyIntegrationPointFunctionVector* _rho_L_h;
	MyIntegrationPointFunctionVector* _rho_G_h;
	MyIntegrationPointFunctionVector* _rho_G_w;
	MyIntegrationPointFunctionVector* _dPGh_dPg;
	MyIntegrationPointFunctionVector* _dPcdSg;//
	MyIntegrationPointFunctionVector* _Func_C;//
	MyIntegrationPointFunctionVector* _vel_L;//
	MyIntegrationPointFunctionVector* _vel_G;//
	MyIntegrationPointFunctionVector* _x_G_air;
	MyIntegrationPointFunctionVector* _x_G_vapor;

	My3DIntegrationPointFunctionVector*  _vel_L_3D;
	My3DIntegrationPointFunctionVector*  _vel_G_3D;
	
	/**
	* store the local matrix M and K 
	*/
	std::vector<LocalMatrix> _elem_M_matrix; 
	std::vector<LocalMatrix> _elem_K_matrix;
		
	/**
	  * derivative 
	  * for each node 8 by 2 matrix. 
	  */
	MyNodalFunctionMatrix* _mat_secDer;
	/**
	* define a vectior to store the value of Omega M and Characteristic function
	* therefore the size should be 3.
	*/
	MyNodalFunctionVector* _vec_tempVar;

	/**
	  * the id of the msh applied in this process
	  */
	size_t _msh_id;

	/**
	  * number of chemical components, 
	  * in the current UserModule, it is fixed to two. 
	  */
	const size_t _n_Comp;

	/**
	  * number of fluid phases,
	  * in the current UserModule, it is also fixed to two.
	  */
	const size_t _n_Phases;

	EOS_NonIso_HeatPipe_CapPressureForm* _EOS;
	EOS_EVL_METHOD m_EOS_EVL_METHOD;
	double real_x[3], gp_x[3];
};

#include "FunctionCMP_NonIso_CapPressureForm.hpp"

#endif  // end of ifndef
