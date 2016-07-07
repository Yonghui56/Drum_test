/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EOS_H2_H2O_ISOT.h
 *
 * Created on 2014-05-12 by Yonghui HUANG
 */

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_NonIsoTotalMassEXPForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_NonIsoTotalMassEXPForm : public AbstractEOS_NonIsoTotalMassEXPForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_NonIsoTotalMassEXPForm() : AbstractEOS_NonIsoTotalMassEXPForm(2)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_NonIsoTotalMassEXPForm()
    {

    };

	
	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */


	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
		//
		T_L = env_conditon(2);
	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		// getting unknowns
		
		double Sg = vec_unknowns(0);
		double rho_G_w = vec_unknowns(1);
		double PGW = get_P_sat(P_L,T_L);
		// calculating residual
		res(0) = X_L - rho_G_w*Sg - rho_l_std*(1-Sg);
		res(1) = std::min(1 - Sg, fac*(PGW*M_L / R / T_L - rho_G_w));
	};
	/*
	* Calculate the Jacobian Matrix
	*/
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_G_w = vec_unknowns(1);
		double PG_w (0.0);
		PG_w = get_P_sat(P_L,T_L);
		double F1 = 1 - Sg;
		double G1 = fac*(PG_w*M_L / R / T_L - rho_G_w);//C_v*PG_h; /)/
		// evaluate J
		J.setZero();
		J(0, 0) = rho_l_std - rho_G_w;//-----dF(1)/dSg
		J(0, 1) = -Sg;

		if (F1 <= G1) {
			//std::cout << "phase transition happens" << std::endl;
			//std::cout << X_L << std::endl;
			J(1, 0) = -1.0;
			J(1, 1) = 0.0;
		}
		else{
			J(1, 0) = 0.0;
			J(1, 1) = -fac*1.0;
		}
	};
	virtual double LocalF1(double Sg, double rho_G_w)
	{
		return X_L - Sg*rho_G_w - rho_l_std*(1 - Sg);
	}
	virtual double LocalF2(double Sg, double rho_G_w)
	{
		double PG_w = get_P_sat(P_L, T_L);
		return std::min(1 - Sg, PG_w*M_L / R / T_L - rho_G_w);
	}

	/*
	* Calculate the Jacobian Matrix
	*/
	virtual void calc_Jacobian_FD(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());
		double Sg = vec_unknowns(0);
		double rho_G_w = vec_unknowns(1);
		//double PG_w(0.0);
		//PG_w = get_P_sat(P_L, T_L);	
		// evaluate J
		J.setZero();
		J(0, 0) = LocalF1(Sg + eps, rho_G_w) - LocalF1(Sg - eps, rho_G_w) / 2 / eps;
		J(0, 1) = LocalF1(Sg , rho_G_w + eps) - LocalF1(Sg , rho_G_w - eps) / 2 / eps;
		J(1, 0) = LocalF2(Sg + eps, rho_G_w) - LocalF2(Sg - eps, rho_G_w) / 2 / eps;
		J(1, 1) = LocalF2(Sg, rho_G_w + eps) - LocalF2(Sg, rho_G_w  - eps) / 2 / eps;
		
	};
	
	virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
	{
		double Pc = 0.0;
		double S_eff = getEffectSat_lbySg(Sg);
		double p0 = 6.324555320336758e+05;
		Pc = p0*0.05878*((1.263*(1 - S_eff) - 2.120)*(1 - S_eff) + 1.417) * (1 - S_eff);
		if (Sg > 1 - S_lr) {
			double y = p0*0.05878*((1.263*1.0 - 2.120)*1.0 + 1.417)*1.0;
			double m = p0*0.05878*((3 * 1.263*1.0 - 2 * 2.120)*1.0 + 1.417);// *(1 / (1 - S_lr));
			return -S_eff*m + y;
		}
		else if (Sg < 0.0) {
			double y = 0.0;
			double m = p0*0.05878*1.417;// *(1 / (1 - S_lr));
			return (1 - S_eff)*m + y;
		}
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{
		double dPcdSg(0.0);
		double S_eff = getEffectSat_lbySg(Sg);

		double p0 = 6.324555320336758e+05;

		if (Sg > 1.0 - S_lr)
			S_eff = 0.0;
		dPcdSg = p0*0.05878*((3 * 1.263*(1 - S_eff) - 2 * 2.120)*(1 - S_eff) + 1.417)*(1 / (1 - S_lr));
		if (Sg < 0.0){
			dPcdSg = p0*0.05878*1.417*(1 / (1 - S_lr));
		}

		return dPcdSg;
	}
	/**
	* calculate the effective saturation of Liquid phase by Gas saturation
	*/
	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		/*if (Sg < S_gr)
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		else if (Sg>(1 - S_lr))
		EffectSat_l = 0.0;
		else*/
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	/**
	* based on primary variables to mean_pressure get the value of PL ---Liquid phase pressure
	*/
	virtual double getPL(double PG, double PC)
	{
		
		return PG-PC;
	}


	virtual double getPGH(double PG, double T)
	{
		double P_sat = get_P_sat(PG,T);//MODIFY
		double PGH = PG - P_sat;
		return PGH;
	}

	virtual double get_P_sat(double PG,double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		//double PC = getPcbySg(0.01);
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//P_sat *= exp(-(PC - 0.7 * 101325)*M_L / R / T / 1000);
		//return std::min(P_sat,PG);
		return P_sat;
	}

	

	
	/**
	* based on Sg  get the value of hydrogen mass density in the GAS phase
	*/
	virtual double get_RHO_G_H(double PG, double T, double rho_G_w)
	{
		
		double RHO_G_H= PG*M_G / R / T - M_G*rho_G_w / M_L;
		return RHO_G_H;
		/*
		if (RHO_G_H < 0){
			return 0.0;
		}
		else
			return RHO_G_H;
		*/
	}


	virtual double Deriv_dPsat_dT(double PG, double T)
	{

		// Here unit of T is Celsius;
		double dPsat_dT(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		dPsat_dT = P_0*(M_L*h_wg / R)*(1. / T / T)*exp(((1. / T_0) - (1. / T))*M_L*h_wg / R);
		//double P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//if (T > T_0+2)
		//dPsat_dT = 0.0;
		return dPsat_dT;
	}
	
	/**
	* based on pg  get the value of Deriv_dSgdP ---
	*/
	
	virtual double Deriv_dSgdP(double rho_G_w, double drhoGWdP, double Sg)
	{
		
		return -drhoGWdP*Sg / (rho_G_w-rho_l_std);
	}
	virtual double Deriv_dSgdX(double rho_G_w, double drhoGWdX, double Sg)
	{
		
		
		return (1 - drhoGWdX*Sg) / (rho_G_w - rho_l_std);
	}

	virtual double Deriv_dSgdT(double rho_G_w, double drhoGWdT, double Sg)
	{


		return (- drhoGWdT*Sg) / (rho_G_w - rho_l_std);
	}
	
	/**
	* Deriv rho_G_w based on P mean pressure
	*/
	virtual double Deriv_drhoGw_dP(double dPGWdPG, size_t chara)
	{
		return dPGWdPG*M_L*chara / R / T_L;
	}

	/**
	* Deriv rho_G_w based on X --total mass density
	*/
	virtual double Deriv_drhoGw_dX(size_t chara)
	{
		return 1 - chara;
	}

	/**
	* Deriv rho_G_w based on T --TEMPERATURE
	*/
	virtual double Deriv_drhoGw_dT(double PGW, double dPGWdT, size_t chara)
	{
		return chara*(dPGWdT*M_L / R / T_L - PGW*M_L / R / T_L / T_L);
	}

	/**
	* Deriv rho_G_w based on P mean pressure
	*/
	virtual double Deriv_drhoGh_dP(double drhoGw_dP, double T)
	{
		
		return (M_G / R / T - M_G*drhoGw_dP / M_L);
	}

	/**
	* Deriv rho_G_w based on X --total mass density
	*/
	virtual double Deriv_drhoGh_dX(double drhoGw_dX)
	{
		return -M_G*drhoGw_dX / M_L;
	}

	/**
	* Deriv rho_G_w based on T --TEMPERATURE
	*/
	virtual double Deriv_drhoGh_dT(double drhoGw_dT, double PG, double T)
	{
		return -PG*M_G / R / T / T - M_G*drhoGw_dT / M_L;
	}

	

	/**
	* water vapor pressure of pure water
	* Function of Temperature
	* comes from Amaziane's paper
	*/
	virtual double get_Vapor_Pressure(double T)
	{
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		//double PC = getPcbySg(0.01);
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		return P_sat;
	}

	/**
	* OVERALL HEAT CAPACITY
	*/
	virtual double get_overall_Heat_Capacity(double Sg, double lambda_pm_dry = 0.582, double lambda_pm_wet = 1.14)
	{
		double lambda_pm(0);
		lambda_pm = lambda_pm_dry + pow(1 - Sg, 1 / 2)*(lambda_pm_wet - lambda_pm_dry);
		if (Sg < 0)
			lambda_pm = lambda_pm_wet;
		else if (Sg > 1)
			lambda_pm = lambda_pm_dry;
		return lambda_pm;

	}

private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double R = 8.314;
	const double S_gr = 0.0;//residual saturation of the gas phase
	const double S_lr = 0.15;//residual saturation of the liquid phase
	//const double P_vapor = 0.0;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	const double eps = 1e-8;
	
	const double fac = 1.0;

	const double rho_l_std = 1000;
	double P_L;
	double X_L;
	double T_L;

	double _alpha;
	double _beta;

	double xi = 1e-5;
};


