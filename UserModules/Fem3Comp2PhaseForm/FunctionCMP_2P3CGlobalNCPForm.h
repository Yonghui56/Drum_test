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

#ifndef FUNCTIONCMP_2P3CGLOBALNCPFORM_H
#define FUNCTIONCMP_2P3CGLOBALNCPFORM_H

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
#include "UserModules/Fem3Comp2PhaseForm/NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler.h"
#include "UserModules/Fem3Comp2PhaseForm/NonLinearCMP_2P3CGlobalNCPForm_JacobianLocalAssembler.h"
#include "UserModules/Fem3Comp2PhaseForm/NonLinearCMP_2P3CGlobalNCPForm_Nested_EOS_NRIterationStepInitializer.h"
#include "Fem_CMP_2P3CGlobalNCPForm_Solution.h"
#include "SingleStep2P3CGlobalNCPFormNonlinear.h"
#include "SingleStepCMP_2P3CGlobalNCPForm.h"
#include "TemplateFemEquation_2P3CGlobalNCPForm.h"
#include "EOS_2P3CGlobalNCPForm.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionCMP_2P3CGlobalNCPForm
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

	typedef FunctionCMP_2P3CGlobalNCPForm       MyFunctionData;    // FEM problem class
	typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
	typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver

	// memory for discretized concentration vector
	typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename FemLib::FemNodalFunctionMatrix<MyDiscreteSystem>::type MyNodalFunctionMatrix; 
	typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
	typedef typename NumLib::TXWrapped3DVectorFunction<MyIntegrationPointFunctionVector> My3DIntegrationPointFunctionVector;
	// local assembler
	typedef NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar, FunctionCMP_2P3CGlobalNCPForm>      MyNonLinearAssemblerType;
	typedef NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar, FunctionCMP_2P3CGlobalNCPForm> MyNonLinearResidualAssemblerType;
	typedef NonLinearCMP_2P3CGlobalNCPForm_JacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType;
	typedef NonLinearCMP_2P3CGlobalNCPForm_Nested_EOS_NRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
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
	typedef TemplateFemEquation_2P3CGlobalNCPForm<
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
	> MyNonLinearCMP2P3CGlobalNCPFormProblemType;

	/**
	* Solution algorithm definition
	*/
	
	typedef SingleStep2P3CGlobalNCPFormNonlinear<
		MyNonLinearCMP2P3CGlobalNCPFormProblemType,
		MyFunctionData,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;

	
	/**
	  * the general CompMultiPhase problem part
	  */
	typedef SolutionLib::Fem_CMP_2P3CGlobalNCPForm_Solution<MyDiscreteSystem> MyCMP2P3CGlobalNCPFormProblemType;
	typedef SolutionLib::SingleStepCMP_2P3CGlobalNCPForm<MyFunctionData,
		MyCMP2P3CGlobalNCPFormProblemType,
		MyNonLinearCMP2P3CGlobalNCPFormProblemType,
		MyNonLinearSolutionType> MyCMP2P3CGlobalNCPFormSolution;
	typedef typename MyCMP2P3CGlobalNCPFormProblemType::MyVariable MyVariableCMP2P3CGlobalNCPForm;

	/**
	  * constructor
	  */
	FunctionCMP_2P3CGlobalNCPForm()
		: Process("CMP_2P3CGlobalNCPForm", 0, 7),
		_feObjects(0), _n_Comp(5), _n_Phases(2)
	{
		_EOS = new EOS_2P3CGlobalNCPForm();
		m_EOS_EVL_METHOD = EOS_EVL_ANALYTICAL;
	};

	/**
	  * destructor, reclaim the memory
	  */
	virtual ~FunctionCMP_2P3CGlobalNCPForm()
	{
		BaseLib::releaseObject(_P);
		BaseLib::releaseObject(_X1);
		BaseLib::releaseObject(_X2);
		BaseLib::releaseObject(_X3);
		BaseLib::releaseObject(_X4);
		BaseLib::releaseObject(_X5);
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
      * update the value of primary variable X1
      */
    void set_X1(MyNodalFunctionScalar* new_nodal_values)
    {
        std::size_t node_idx;
        double nodal_X;

        for (node_idx = _X1->getDiscreteData()->getRangeBegin();
             node_idx < _X1->getDiscreteData()->getRangeEnd();
             node_idx++)
        {
            nodal_X = new_nodal_values->getValue(node_idx);
            this->_X1->setValue(node_idx, nodal_X);
        }
    };
	/**
	* update the value of primary variable X2
	*/
	void set_X2(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X;

		for (node_idx = _X2->getDiscreteData()->getRangeBegin();
			node_idx < _X2->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X = new_nodal_values->getValue(node_idx);
			this->_X2->setValue(node_idx, nodal_X);
		}
	};
	/**
	* update the value of primary variable X1
	*/
	void set_X3(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X;

		for (node_idx = _X3->getDiscreteData()->getRangeBegin();
			node_idx < _X3->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X = new_nodal_values->getValue(node_idx);
			this->_X3->setValue(node_idx, nodal_X);
		}
	};

	/**
	* update the value of primary variable X1
	*/
	void set_X4(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X;

		for (node_idx = _X4->getDiscreteData()->getRangeBegin();
			node_idx < _X4->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X = new_nodal_values->getValue(node_idx);
			this->_X4->setValue(node_idx, nodal_X);
		}
	};
	/**
	* update the value of primary variable X1
	*/
	void set_X5(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X;

		for (node_idx = _X5->getDiscreteData()->getRangeBegin();
			node_idx < _X5->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X = new_nodal_values->getValue(node_idx);
			this->_X5->setValue(node_idx, nodal_X);
		}
	};

	/**
	* update the value of primary variable X
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
	void set_X1_node_values(std::size_t node_idx, double node_value){ _X1->setValue(node_idx, node_value); };
	/**
	* set molar fraction nodal value
	*/
	void set_X2_node_values(std::size_t node_idx, double node_value){ _X2->setValue(node_idx, node_value); };
	/**
	* set molar fraction nodal value
	*/
	void set_X3_node_values(std::size_t node_idx, double node_value){ _X3->setValue(node_idx, node_value); };
	/**
	* set molar fraction nodal value
	*/
	void set_X4_node_values(std::size_t node_idx, double node_value){ _X4->setValue(node_idx, node_value); }
	/**
	* set molar fraction nodal value
	*/
	void set_X5_node_values(std::size_t node_idx, double node_value){ _X5->setValue(node_idx, node_value); };
	/**
	* set SATURATION nodal value
	*/
	void set_S_node_values(std::size_t node_idx, double node_value){ _S->setValue(node_idx, node_value); };
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
	MyIntegrationPointFunctionVector* get_PL(void) { return _PL; };
	/**
	* return the molar density of H2
	*/
	MyIntegrationPointFunctionVector* get_rhoGh(void) { return _rho_G_h; };

	/**
	* return the molar density of H2O
	*/
	MyIntegrationPointFunctionVector* get_rhoGw(void) { return _rho_G_w; };

	/**
	* return the molar density of CH4
	*/
	MyIntegrationPointFunctionVector* get_rhoGch4(void) { return _rho_G_ch4; };



	/**
	* return the dissolved gas mass density
	*/
	MyNodalFunctionScalar* get_rhoLh(void) { return _rho_L_h; };
	

	MyNodalFunctionScalar* get_dPGh_dPg(void) { return _dPGh_dPg; };
	MyNodalFunctionScalar* get_dPc_dSg(void) { return _dPcdSg; };

	MyNodalFunctionScalar* get_Func_C(void)  { return _Func_C; };
	/**
	* return the matrix of derivative of secondary variable based on P and X
	*/
	MyNodalFunctionMatrix* get_mat_secDer(void) { return _mat_secDer; };
	MyNodalFunctionVector* get_vec_tempVar(void) { return _vec_tempVar; };

	LocalVector get_vec_Res(void) { return _vec_Res; };
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
	virtual MyCMP2P3CGlobalNCPFormSolution* getSolution() { return _solution; };

	/**
	  * output the result of current solution
	  */
	virtual void output(const NumLib::TimeStep &time);

private:
	DISALLOW_COPY_AND_ASSIGN(FunctionCMP_2P3CGlobalNCPForm);

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
	MyNonLinearCMP2P3CGlobalNCPFormProblemType*     _non_linear_problem;

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
	MyCMP2P3CGlobalNCPFormProblemType* _problem;

	/**
	  * solution class for the component based multiphase flow problem
	  */
	MyCMP2P3CGlobalNCPFormSolution*    _solution;

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
	MyNodalFunctionScalar* _P;

	/**
	* Primary variable 2): vector of molar fraction values of the light component on each node
	*/
	
	MyNodalFunctionScalar* _X1;//H2
	/**
	* Primary variable 3): vector of molar fraction values of the light component on each node
	*/
	MyNodalFunctionScalar* _X2;//CH4
	/**
	* Primary variable 4): vector of molar fraction values of the light component on each node
	*/
	MyNodalFunctionScalar* _X3;//VAPOR
	/**
	* Primary variable 4): vector of molar fraction values of the light component on each node
	*/
	MyNodalFunctionScalar* _X4;//AIR
	/**
	* Primary variable 5): vector of molar fraction values of the light component on each node
	*/
	MyNodalFunctionScalar* _X5;//CO2
	/**
	* Primary variable 6): vector of saturation values on each noda
	*/
	MyNodalFunctionScalar* _S;

	/**
	*		Vector _output
	*		store the eight values of the secondary variables
	*/
	/** 
	  * secondary variable --saturation
	  */

	MyIntegrationPointFunctionVector* _PL;
	MyIntegrationPointFunctionVector* _rho_G_h;
	MyIntegrationPointFunctionVector* _rho_G_w;
	MyIntegrationPointFunctionVector* _rho_G_ch4;

	MyNodalFunctionScalar* _PG;
	
	MyNodalFunctionScalar* _PC;
	MyNodalFunctionScalar* _rho_L_h;
	MyNodalFunctionScalar* _dPGh_dPg;
	MyNodalFunctionScalar* _dPcdSg;//
	MyNodalFunctionScalar* _Func_C;//
	MyNodalFunctionScalar* _X_L_h;
	MyNodalFunctionScalar* _X_L_c;
	MyNodalFunctionScalar* _X_L_w;
	
	/**
	* store the local matrix M and K 
	*/
	std::vector<LocalMatrix> _elem_M_matrix; 
	std::vector<LocalMatrix> _elem_K_matrix;
	std::vector<LocalMatrix> _elem_J_matrix;
	/**
	* store the Jacobian matrix for complementary function 1~~~~~~~
	*/
	MathLib::SparseMatrix _mat_Jacob;
	LocalVector _vec_Res;

	MathLib::SparseMatrix _mat_Jacob_L;
	LocalVector _vec_Res_L;
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

	EOS_2P3CGlobalNCPForm* _EOS;
	EOS_EVL_METHOD m_EOS_EVL_METHOD;
	const double Hen_L_h = 7.26e+9; //Henry constant in [Pa]
	const double Hen_L_c = 4.13e+9; //Henry constant in [Pa]
	const double Hen_L_air = 9.077e+9; //Henry constant in [Pa]
	const double Hen_L_co2 = 0.163e+9; //Henry constant in [Pa]
	const double rho_L_std = 1000;
	const double M_G = 0.002;
	const double M_L = 0.018;
	const double T_0 = 303.15;
};

#include "FunctionCMP_2P3CGlobalNCPForm.hpp"

#endif  // end of ifndef