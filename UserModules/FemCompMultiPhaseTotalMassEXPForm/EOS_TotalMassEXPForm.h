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
#include "AbstractEOS_TotalMassEXPForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_TotalMassEXPForm : public AbstractEOS_TotalMassEXPForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_TotalMassEXPForm() : AbstractEOS_TotalMassEXPForm(3)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_TotalMassEXPForm()
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
	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	
	

	virtual double getS_bar(double Sg)
	{
		double S_bar = 0.0;
		
		S_bar = S_gr + (1 - xi)*(Sg - S_gr) + 0.5*xi*(1 - S_gr - S_lr);
		return S_bar;
	}
	/**
	*  van Genuchten capillary pressure-saturation Model
	*/
	virtual double getPc_vG_Sg(double Sg)
	{
		double Pc_vG = 0.0;
		double S_le = 0.0;
	    //effective saturation
		S_le =  (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		//Pc_vG = P_r*(S_le. ^ (-1 / m) - 1). ^ (1 / n);
		Pc_vG = P_r*pow(pow(S_le, -1 / m) - 1, 1 / n);
		return Pc_vG;
	}
	/**
	* regularized van Genuchten capillary pressure-saturation Model
	*/
	virtual double getPc_bar_vG_Sg(double Sg)
	{
		double Pc_bar_vG = 0.0;
		double S_bar;
		S_bar = getS_bar(Sg);
		Pc_bar_vG = getPc_vG_Sg(S_bar) - getPc_vG_Sg(S_gr + (1 - S_gr - S_lr)*xi / 2);
		return Pc_bar_vG;
	}
	/**
	* derivative dPCdS based on standard van Genuchten capillary pressure-saturation Model
	*/
	virtual double get_dPCdS_vG(double Sg)
	{
		double dPcdSg = 0.0;
		double S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		dPcdSg = P_r*(1 / (m*n))*(1 / (1 - S_lr - S_gr))*pow(pow(S_le, (-1 / m)) - 1, (1 / n) - 1)*pow(S_le, (-1 / m)) / S_le;
		return dPcdSg;
	}
	/**
	* derivative dPCdS based on regularized van Genuchten capillary pressure-saturation Model
	*/
	virtual double get_dPCdS_vG_bar(double Sg)
	{
		double dPCdS(0.0);
		double S_bar = 0.0;
		S_bar = getS_bar(Sg);
		dPCdS = get_dPCdS_vG(S_bar)*(1 - xi);
		return dPCdS;
	}

	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/

	virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
	{
		double Pc = 0.0;
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			Pc = getPc_bar_vG_Sg(Sg);
		}
		else if (Sg < S_gr)
		{
			Pc = getPc_bar_vG_Sg(S_gr) + get_dPCdS_vG_bar(S_gr)*(Sg - S_gr);
		}
		else
		{
			Pc = getPc_bar_vG_Sg(1 - S_lr) + get_dPCdS_vG_bar(1 - S_lr)*(Sg - 1 + S_lr);
		}
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{

		double dPcdSg(0.0);
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			dPcdSg = get_dPCdS_vG_bar(Sg);
		}
		else if (Sg < S_gr)
		{
			dPcdSg = get_dPCdS_vG_bar(S_gr);
		}
		else
		{
			dPcdSg = get_dPCdS_vG_bar(1 - S_lr);
		}
		return dPcdSg;
	}
	/**
	* calculate the effective saturation of Liquid phase by Gas saturation
	*/
	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;
		
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	
	/**
	* based on primary variable mean_pressure to get the value of PG ---Gas phase pressure
	*/
	virtual double getPG(double Sg)
	{
		double PG = 0.0;
		double PC = 0.0;
		PC= getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PG = P_L + PC;// (1 - gamma)*PC;
		return PG;
	}

	/**
	* based on primary variables to mean_pressure get the value of PL ---Liquid phase pressure
	*/
	virtual double getPL(double Sg)
	{
		double PL = 0.0;
		double PC = 0.0;
		PC = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PL = P_L - gamma*PC;
		return PL;
	}


	virtual double getPGH(double PG, double T)
	{
		double P_sat = get_P_sat(T);
		double PGH(0.0);
		//PGH = (Hen*M_L* PG - rho_l_std + std::sqrt(pow(Hen* M_L* PG + rho_l_std, 2) - 4 * Hen *M_L* P_sat* rho_l_std)) / (2 * Hen* M_L);
		PGH = PG;
		return PGH;
	}

	virtual double get_P_sat(double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		double PC = getPcbySg(0.01);
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//P_sat *= exp(-(PC - 0.7 * 101325)*M_L / R / T / 1000);
		return P_sat;
	}

	virtual double get_Sat(double PG, double X)
	{
		double PGH = getPGH(PG, T0);
		double C_L_h = Hen*M_G;
		double C_G_h = M_G/R/T0;
		if (X <= PGH*Hen*M_G)
			return 0.0;
		else if (X >= PGH*M_G / R / T0)
			return 1.0;
		else
		{
			return (X / PGH - Hen*M_G) / (C_G_h - C_L_h);
		}
	}

	virtual double get_alpha(double Sg, double dPGHdPG)
	{
		return ((C_G_h - C_L_h)*Sg + C_L_h)*dPGHdPG;
	}
	 
	virtual double get_beta(double PG)
	{
		double PGH = getPGH(PG, T0);
		return (C_G_h - C_L_h)*PGH;
	}

	/**
	* weighted function 
	*/
	virtual double getweightedFunc(double Sg)
	{
		double gamma=0.0;
		
		/*if (Sg >= 1)
			gamma = 1;
		else if (Sg <= 0)
			gamma = 0;
		else
			gamma = Sg*(2 - Sg);
		*/
		return gamma;
	}
	/**
	* Derivation of weighted function based on Sg
	*/
	virtual double get_deriv_weight_Func(double Sg)
	{
		double dgamma_dSg(0.0);
		/*
		dgamma_dSg = 2 - 2 * Sg;
		if (Sg >= 1)
			dgamma_dSg = 0;
		else if (Sg <= 0)
			dgamma_dSg = 0;
			*/
		return dgamma_dSg;
	}

	/**
	* based on Sg  get the value of dissoved hydrogen mass density in the LIQUID phase
	*/
	virtual double get_RHO_L_H(double PG, double X, double Sg)
	{
		double rho_L_h(0.0);
		//double PG_H(0.0);
		//PG_H = Function_PG_h(PG);
		double PGH = getPGH(PG,T0);
		rho_L_h = std::min(C_L_h*PGH, X);//Henry law
		return rho_L_h;
	}
	/**
	* based on Sg  get the value of hydrogen mass density in the GAS phase
	*/
	virtual double get_RHO_G_H(double PG, double X, double sg)
	{
		double rho_G_h(0.0);
		double PGH = getPGH(PG, T0);
		rho_G_h = std::max(C_G_h*PGH, X);//ideal gas law
		return rho_G_h;
	}
	/**
	* based on Sg  get the value of water vapor mass density in the GAS phase
	*/
	virtual double get_RHO_G_W(double PG, double rho_G_h)
	{
		
		
		return C_G_w*PG - C_G_w*rho_G_h / C_G_h;
	}
	/**
	* based on pg  get the value of P_G^h ---PARTIAL GAS PRESSURE OF HYDROGEN
	*/

	virtual double Function_PG_h(double PG)
	{
		double PG_H(0.0);
		double P_vapor(0.0);
		P_vapor = get_Vapor_Pressure(T0-273);
		//PG_H = (1 / (2 * C_h*M_L))* (C_h*M_L*PG - M_G* rho_L_std + sqrt(pow((-C_h*M_L*PG + M_G*rho_L_std), 2) - 4 * C_h*M_L*(-M_G*PG*rho_L_std + M_G* P_vapor* rho_L_std)));
		PG_H = PG;
		return PG_H;

	}


	/**
	* derivative dP_G^HdPG
	*/
	virtual double Deriv_PGH_dPG(double PG, double PGH, double T=303 )
	{
		//double test=(getPGH(PG + EPS*PG, T) - getPGH(PG - EPS*PG, T)) / 2 / EPS/PG;
		double P_sat = get_P_sat(T);
		return 1.0;// 1 / (1 - (Hen * M_L * rho_l_std *P_sat) / pow((rho_l_std + Hen * M_L * PGH), 2));
	}
	
	

	/**
	* Derivation of PG based on saturation
	*/
	virtual double Deriv_dPGdSg(double Sg, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double dPcdSg(0.0);
		dPcdSg = Deriv_dPCdS(Sg);
		double PC(0.0);
		PC = getPcbySg(Sg);
		double Gamma(0.0);
		Gamma = getweightedFunc(Sg);
		double dPGdSg(0.0);
		dPGdSg = dPcdSg*(1 - Gamma) - PC*get_deriv_weight_Func(Sg);
		//if (Sg < 1e-14)
			//dPGdSg = 0.0;
		return dPGdSg;
	}
	
	virtual double Deriv_dSgdP(double alpha, double beta, size_t chara)
	{
		
		
		double dPG_dSg(0.0);
		double dSgdP(0.0);
		dSgdP = -alpha*chara/ (beta + alpha*dPG_dSg);
		return dSgdP;
	}
	virtual double Deriv_dSgdX(double alpha, double beta, size_t chara)
	{
		
		double dPG_dSg(0.0);
		double dSgdX(0.0);
		dSgdX = chara / (beta + alpha*dPG_dSg);
		return dSgdX;
	}
	
	/**
	* Deriv rho_L_h based on mean pressure
	*/
	virtual double Deriv_drhoLh_dP(double dPGHdPG, size_t chara)
	{
		return C_L_h*dPGHdPG*chara;
	}

	/**
	* Deriv rho_L_h based on total mass density
	*/
	virtual double Deriv_drhoLh_dX(size_t chara)
	{
		return 1 - chara;
	}
	/**
	* Deriv rho_G_h based on P mean pressure
	*/
	virtual double Deriv_drhoGh_dP(double dPGHdPG, size_t chara)
	{
		return C_G_h* dPGHdPG*chara;
	}

	/**
	* Deriv rho_G_h based on X --total mass density
	*/
	virtual double Deriv_drhoGh_dX(size_t chara)
	{
		return 1 - chara;
	}

	

	/**
	* water vapor pressure of pure water
	* Function of Temperature
	* comes from Amaziane's paper
	*/
	virtual double get_Vapor_Pressure(double T)
	{
		// Here unit of T is Celsius;
		double P_vapor(0.0);
		/*double tmp(0.0);
		tmp = 2.786 + 0.031514*T - 1.2373e-4*pow(T, 2) + 4.2267e-7*pow(T, 3) - 8.1308e-10*pow(T, 4);
		P_vapor = pow(10, tmp);*/
		return P_vapor;
	}

private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double R = 8.314;
	const double T0 = 303;// [K]
	const double Hen = 7.65e-6; //Henry constant
	//const double P_vapor = 0.0;
	const double M_G = 0.002;
	const double M_L = 0.01;
	const double C_L_h = Hen*M_G;
	const double C_G_h = M_G / R / T0;
	const double C_G_w = 0.0;// M_L / R / T0;
	const double eps = 1e-12;
	const double rho_L_std = 1000;
	
	/**
	* parameters for van-Genuchten capillary pressure-saturation model
	*/
	const double S_gr = 0.0;//residual saturation of the gas phase
	const double S_lr = 0.4;//residual saturation of the liquid phase
	const double P_r = 2e+6;//pa
	const double n = 1.49;
	const double m = 1 - 1 / n;
	const double rho_l_std = 1000;
	double P_L;
	double X_L;

	double _alpha;
	double _beta;

	double xi = 1e-5;
};


