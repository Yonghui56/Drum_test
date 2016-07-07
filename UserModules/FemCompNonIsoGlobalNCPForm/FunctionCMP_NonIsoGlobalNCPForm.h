/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.h
*
* Created on    2014-09-12 by Yonghui HUANG
*/

#ifndef FUNCTIONCMP_NONISOGLOBALNCP_H
#define FUNCTIONCMP_NONISOGLOBALNCP_H

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
#include "UserModules/FemCompNonIsoGlobalNCPForm/NonLinearCMP_NonIsoGlobalNCPForm_TimeODELocalAssembler.h"
#include "UserModules/FemCompNonIsoGlobalNCPForm/NonLinearCMP_NonIsoGlobalNCPForm_JacobianLocalAssembler.h"
#include "UserModules/FemCompNonIsoGlobalNCPForm/NonIsoGlobalNCPForm_Nested_EOS_NRIterationStepInitializer.h"
#include "Fem_CMP_NonIsoGlobalNCP_Solution.h"
#include "SingleStepNonIsoGlobalNCPFormNonlinear.h"
#include "SingleStepCMP_NonIsoGlobalNCPForm.h"
#include "TemplateFemEquation_NonIsoGlobalNCPForm.h"
#include "EOS_NonIsoGlobalNCPForm.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionCMP_NonIsoGlobalNCPForm
	: public ProcessLib::Process
{
public:
	// input variable is velocity
	enum In {};//Velocity = 0 
	// no output variable
	//enum Out {  = 0, Total_Mass_Density = 0 };

	enum EOS_EVL_METHOD { 
		EOS_EVL_FIN_DIFF, 
		EOS_EVL_ANALYTICAL
	};

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionCMP_NonIsoGlobalNCPForm       MyFunctionData;    // FEM problem class
	typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
	typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver

	// memory for discretized concentration vector
	typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename FemLib::FemNodalFunctionMatrix<MyDiscreteSystem>::type MyNodalFunctionMatrix; 

	// local assembler
	typedef NonLinearCMP_NonIsoGlobalNCPForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar, FunctionCMP_NonIsoGlobalNCPForm>      MyNonLinearAssemblerType;
	typedef NonLinearCMP_NonIsoGlobalNCPForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar, FunctionCMP_NonIsoGlobalNCPForm> MyNonLinearResidualAssemblerType;
	typedef NonLinearCMP_NonIsoGlobalNCPForm_JacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType;
	typedef NonIsoGlobalNCPForm_Nested_EOS_NRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer> MyDiscreteNonlinearSolverFactory;

	/**
	* nonlinear coupled part solving xi
	* Equation definition
	*/
	/*typedef SolutionLib::TemplateFemEquation<
		MyDiscreteSystem,
		MyLinearSolver,
		MyNonLinearAssemblerType,
		MyNonLinearResidualAssemblerType,
		MyNonLinearJacobianAssemblerType
	> MyNonLinearEquationType;
	*/
	typedef TemplateFemEquation_NonIsoGlobalNCPForm<
		MyDiscreteSystem,
		MyLinearSolver,
		MyNonLinearAssemblerType,
		MyFunctionData,
		MyNonLinearResidualAssemblerType,
		MyNonLinearJacobianAssemblerType
	> MyNonLinearEquationType;
	
	/**
	* FEM IVBV problem definition
	*/
	typedef SolutionLib::FemIVBVProblem<
		MyDiscreteSystem,
		MyNonLinearEquationType
	> MyNonLinearCMPNonIsoGlobalNCPFormProblemType;

	/**
	* Solution algorithm definition
	*/
	/*
	typedef SolutionLib::SingleStepFEM<
		MyNonLinearCMPGlobalComplementaryFORMProblemType,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;
	*/
	typedef SingleStepNonIsoGlobalNCPFormNonlinear<
		MyNonLinearCMPNonIsoGlobalNCPFormProblemType,
		MyFunctionData,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;

	
	/**
	  * the general CompMultiPhase problem part
	  */
	typedef SolutionLib::Fem_CMP_NonIsoGlobalNCP_Solution<MyDiscreteSystem> MyCMPNonIsoGlobalNCPFormProblemType;
	typedef SolutionLib::SingleStepCMP_NonIsoGlobalNCPForm<MyFunctionData,
		MyCMPNonIsoGlobalNCPFormProblemType,
		MyNonLinearCMPNonIsoGlobalNCPFormProblemType,
		MyNonLinearSolutionType> MyCMPNonIsoGlobalNCPFormSolution;
	typedef typename MyCMPNonIsoGlobalNCPFormProblemType::MyVariable MyVariableCMPNonIsoGlobalNCPForm;

	/**
	  * constructor
	  */
	FunctionCMP_NonIsoGlobalNCPForm()
		: Process("CMP_NonIso_GlobalNCPForm", 0, 5),//5 primary variable
		_feObjects(0), _n_Comp(2), _n_Phases(2)
	{
		_EOS = new EOS_NonIsoGlobalNCPForm();
		m_EOS_EVL_METHOD = EOS_EVL_ANALYTICAL;
	};

	/**
	  * destructor, reclaim the memory
	  */
	virtual ~FunctionCMP_NonIsoGlobalNCPForm()
	{
		BaseLib::releaseObject(_P);
		BaseLib::releaseObject(_X_G_a);
		BaseLib::releaseObject(_X_G_w);
		BaseLib::releaseObject(_T);
		
		BaseLib::releaseObject(_S);
		BaseLib::releaseObject(_mat_Jacob_C);
		BaseLib::releaseObject(_EOS);
		BaseLib::releaseObject(_PC);
		BaseLib::releaseObject(_PG);
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
    void set_P(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
        double nodal_P;

        for (node_idx = _P->getDiscreteData()->getRangeBegin();
             node_idx < _P->getDiscreteData()->getRangeEnd();
			 node_idx++)
		{
            nodal_P = new_nodal_values->getValue(node_idx);
            this->_P->setValue(node_idx, nodal_P);
		}
	};

    /**
      * update the value of primary variable X
      */
    void set_X_G_a(MyNodalFunctionScalar* new_nodal_values)
    {
        std::size_t node_idx;
        double nodal_X;

		for (node_idx = _X_G_a->getDiscreteData()->getRangeBegin();
			node_idx < _X_G_a->getDiscreteData()->getRangeEnd();
             node_idx++)
        {
            nodal_X = new_nodal_values->getValue(node_idx);
			this->_X_G_a->setValue(node_idx, nodal_X);
        }
    };
	/**
	* update the value of primary variable X
	*/
	void set_X_G_w(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X_G_w;

		for (node_idx = _X_G_w->getDiscreteData()->getRangeBegin();
			node_idx < _X_G_w->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X_G_w = new_nodal_values->getValue(node_idx);
			this->_X_G_w->setValue(node_idx, nodal_X_G_w);
		}
	};

	/**
	* update the value of primary variable S
	*/
	void set_S(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_S;

		for (node_idx = _S->getDiscreteData()->getRangeBegin();
			node_idx < _S->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_S = new_nodal_values->getValue(node_idx);
			this->_S->setValue(node_idx, nodal_S);
		}
	};
	/**
	* update the value of primary variable T
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
	void set_P_node_values(std::size_t node_idx, double node_value){ _P->setValue(node_idx, node_value);  };

	/**
	  * set molar fraction nodal value
	  */
	void set_X_G_a_node_values(std::size_t node_idx, double node_value){ _X_G_a->setValue(node_idx, node_value); };

	/**
	* set molar fraction of water vapor in Gas phase ==nodal value
	*/
	void set_X_G_w_node_values(std::size_t node_idx, double node_value){ _X_G_w->setValue(node_idx, node_value); };

	/**
	* set saturation nodal value
	*/
	void set_S_node_values(std::size_t node_idx, double node_value){ _S->setValue(node_idx, node_value); };
	/**
	* set Temperature nodal value
	*/
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
	std::vector<LocalMatrix>& get_elem_J_matrix(void) { return _elem_J_matrix; };
	/**
	  * return the number of chemical components
	  */
	std::size_t get_n_Comp(void) { return _n_Comp; };


	/**
	* return the value of gas phase saturation
	*/
	
	/**
	* return the value of capillary pressure
	*/
	MyNodalFunctionScalar* get_PC(void) { return _PC; };
	/**
	* return the value of Liquid  pressure
	*/
	MyNodalFunctionScalar* get_PL(void) { return _PL; };
	/**
	* return the value of Gas phase pressure
	*/
	MyNodalFunctionScalar* get_PG(void) { return _PG; };
	/**
	* return the dissolved gas mass density
	*/
	MyNodalFunctionScalar* get_rhoLh(void) { return _rho_L_h; };
	/**
	* return the  gas mass density IN GAS PHASE
	*/
	MyNodalFunctionScalar* get_rhoGh(void) { return _rho_G_h; };
	/**
	* return the  water vapor mass density IN GAS PHASE
	*/
	MyNodalFunctionScalar* get_rhoGw(void) { return _rho_G_w; };

	MyNodalFunctionScalar* get_dPGh_dPg(void) { return _dPGh_dPg; };
	MyNodalFunctionScalar* get_dPc_dSg(void) { return _dPcdSg; };

	MyNodalFunctionScalar* get_Func_C(void)  { return _Func_C; };
	/**
	* return the matrix of derivative of secondary variable based on P and X
	*/
	MyNodalFunctionMatrix* get_mat_secDer(void) { return _mat_secDer; };
	MyNodalFunctionVector* get_vec_tempVar(void) { return _vec_tempVar; };

	LocalVector get_vec_Res(void) { return _vec_Res; };
	//LocalMatrix get_mat_Jacob(void){ return _mat_Jacob; };
	MathLib::SparseMatrix& get_mat_Jacob(void){ return _mat_Jacob; };
protected:
	virtual void initializeTimeStep(const NumLib::TimeStep &time);

	/**
	  * this function is called to exchange output parameters
	  */
	virtual void updateOutputParameter(const NumLib::TimeStep &time);

	/**
	  * get the pointer of solution class for current problem
	  */
	virtual MyCMPNonIsoGlobalNCPFormSolution* getSolution() { return _solution; };

	/**
	  * output the result of current solution
	  */
	virtual void output(const NumLib::TimeStep &time);

private:
	DISALLOW_COPY_AND_ASSIGN(FunctionCMP_NonIsoGlobalNCPForm);

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
	MyNonLinearCMPNonIsoGlobalNCPFormProblemType*     _non_linear_problem;

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
	MyCMPNonIsoGlobalNCPFormProblemType* _problem;

	/**
	  * solution class for the component based multiphase flow problem
	  */
	MyCMPNonIsoGlobalNCPFormSolution*    _solution;

	/**
	  * FEM object
	  */
	FemLib::LagrangeFeObjectContainer* _feObjects;

	/**
	  * convergence checker
	  */
	NumLib::DiscreteDataConvergenceCheck _checker;

	/**
	  * Primary variable 1): vector of Gas Pressure values on each node
	  */
	MyNodalFunctionScalar* _P;

	/**
	* Primary variable 2): vector of molar fraction values of the light component on each node
	*/
	
	MyNodalFunctionScalar* _X_G_a;

	/**
	* Primary variable 3): vector of molar fraction values of the light component on each node
	*/

	MyNodalFunctionScalar* _X_G_w;
	
	/**
	  * Primary variable 4): vector of Saturation values on each noda
	  */
	MyNodalFunctionScalar* _S;
	/**
	* Primary variable 5): vector of Temperature values on each noda
	*/
	MyNodalFunctionScalar* _T;

	/**
	*		Vector _output
	*		store the eight values of the secondary variables
	*/
	/** 
	  * secondary variable --saturation
	  */

	
	MyNodalFunctionScalar* _PG;
	MyNodalFunctionScalar* _PL;
	MyNodalFunctionScalar* _PC;
	MyNodalFunctionScalar* _rho_L_h;
	MyNodalFunctionScalar* _rho_G_h;
	MyNodalFunctionScalar* _rho_G_w;
	MyNodalFunctionScalar* _dPGh_dPg;
	MyNodalFunctionScalar* _dPcdSg;//
	MyNodalFunctionScalar* _Func_C;//

	/**
	* store the local matrix M and K 
	*/
	std::vector<LocalMatrix> _elem_M_matrix; 
	std::vector<LocalMatrix> _elem_K_matrix;
	std::vector<LocalMatrix> _elem_J_matrix;
		
	//LocalMatrix _mat_Jacob;
	MathLib::SparseMatrix _mat_Jacob;
	LocalVector _vec_Res;
	/**
	  * derivative 
	  * for each node 8 by 2 matrix. 
	  */
	MyNodalFunctionMatrix* _mat_Jacob_C;
	/**
	* define a vectior to store the value of Omega M and Characteristic function
	* therefore the size should be 3.
	*/
	MyNodalFunctionVector* _vec_Res_C;

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

	EOS_NonIsoGlobalNCPForm* _EOS;
	EOS_EVL_METHOD m_EOS_EVL_METHOD;
	//const double C_h = 1.530e-8;
	const double rho_L_std = 1000;
	const double M_G = 0.02896;
	const double M_L = 0.018;

};

#include "FunctionCMP_NonIsoGlobalNCPForm.hpp"

#endif  // end of ifndef
