/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file EOS_H2_H2O_ISOT_3x_Complementary.h
*
* Created on 2014-05-30 by Haibing SHAO
*/

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_NonIso_CapPressureForm.h"
#include <cmath>


/**
*  EOS formulation for calculating each secondary variable
*  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg, with only 3 unknows! 
*/
class EOS_NonIso_HeatPipe_CapPressureForm : public AbstractEOS_NonIso_CapPressureForm
{
public:

	/**
	  * Constructor
	  */
	EOS_NonIso_HeatPipe_CapPressureForm() : AbstractEOS_NonIso_CapPressureForm(2)
	{
	};

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_NonIso_HeatPipe_CapPressureForm()
	{

	};

	/**
	  * realization of the eval function.
	  * this function will evaluate the resiudal of the governing equation,
	  * based on the unknown values given.
	  */
	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		PG_L = env_conditon(0);
		// molar fraction of the light component
		PC_L = env_conditon(1);

		T_L = env_conditon(2);

	};	
	
	virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
	{
		double Pc = 0.0;
		double EffectSat_l = getEffectSat_lbySg(Sg);
		Pc = P_r*Lambda_Leverett*(1.417*(1 - EffectSat_l) - 2.12*pow(1 - EffectSat_l, 2.) + 1.263*pow(1 - EffectSat_l, 3.));
		//Pc = P_r*Lambda_Leverett*(1.263*pow(Sg, 3) - 2.12*pow(Sg, 2) + 1.417*Sg);// from lauser's paper
		if (Sg < 0)
			Pc = 0.0;// P_r*Lambda_Leverett*(1.417*Sg);
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{

		double dPcdSg(0.0);
		double EffectSat_l = getEffectSat_lbySg(Sg);
		double dSedSg = 1 / (1 - S_lr);
		dPcdSg = P_r*Lambda_Leverett*(1.263 * 3 * pow(1 - EffectSat_l, 2.)*dSedSg - 2.12 * 2 * (1 - EffectSat_l)* dSedSg + 1.417*dSedSg);
		//dPcdSg = P_r*Lambda_Leverett*(1.263 * 3 * pow(Sg, 2)- 2.12*2*Sg + 1.417 );
		if (Sg < 0)
			dPcdSg = 0.0;// P_r*Lambda_Leverett*1.417;
		return dPcdSg;
	}

	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;

		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	virtual double getPL(double Sg)
	{
		double PL = 0.0;
		PL = PG_L - PC_L;// (1 - gamma)*PC;
		return PL;
	}
	virtual double getrho_G(double PG, double P_vap, double T)
	{
		return (M_G*PG + (M_L - M_G)*P_vap) / R / T;
	}


	virtual double Deriv_dSgdT(double Sg)
	{
		double dSgdT(0);
		
		return dSgdT;
		
	}
	virtual double Deriv_drho_GdP(double T)
	{
		
		return  M_G / R / T;
	}
	virtual double Deriv_drho_GdPC(double T, double dP_GwdPC)
	{

		return  (M_L - M_G)* dP_GwdPC / R / T;
	}
	virtual double Deriv_drho_GdT(double PG, double T, double P_vap, double dP_vapdT)
	{
		//drho_GdT = -((PG - P_vap)*M_G + P_vap*M_L) / R / T + (M_L - M_G)*dP_vapdT / R / T;
		return  -((PG - P_vap)*M_G + P_vap*M_L) / R / T / T +(M_L - M_G)*dP_vapdT / R / T;
	}




	virtual double get_RHO_G_a(double P_G_a, double T)
	{
		double rho_G_a(0.0);
		double C_v = M_G / R / T;
		rho_G_a = P_G_a*C_v;
		if (P_G_a < 0)
			rho_G_a = 0.0;
		return rho_G_a;
	}

	virtual double get_RHO_G_W(double P_gw,double T)
	{
		double rho_G_w(0.0);
		double C_w = M_L / R / T;
		rho_G_w = P_gw*C_w;
		return rho_G_w;
	}	


	/**
	* Deriv rho_G_h based on P mean pressure
	*/
	virtual double Deriv_drhoGa_dPG(double T)
	{
		double C_v = M_G / R / T;
		return C_v;
	}

	/**
	* Deriv rho_G_a based on PC --
	*/
	virtual double Deriv_drhoGa_dPC(double PC, double T)
	{
		double drhoGa_dPC(0.0);
		double dP_gwdPC(0.0);
		double C_v = M_G / R / T;
		dP_gwdPC = Deriv_dPgw_dPC(PC, T);
		drhoGa_dPC = -C_v*dP_gwdPC;
		return drhoGa_dPC;
	}

	virtual double Deriv_drhoGa_dT(double P_G_a,double PC, double T)
	{
		double drhoGa_dT(0.0);
		double C_v = M_G / R / T;
		double dP_gwdT = Deriv_dPgw_dT(PC,T);
		drhoGa_dT = -P_G_a*C_v / T - C_v*dP_gwdT;
		return drhoGa_dT;
	}

	virtual double Deriv_drhoGW_dT(double P_gw,double PC,double T)
	{
		double C_w(0.0);
		C_w = M_L / R / T;
		double drhoGW_dT(0.0);
		double dP_gwdT(0.0);
		dP_gwdT = Deriv_dPgw_dT(PC, T);
		drhoGW_dT = -P_gw*C_w / T + C_w*dP_gwdT;
		return  drhoGW_dT;

	}
	virtual double Deriv_drhoGW_dPC(double PC, double T)
	{
		double drhoGw_dPC(0.0);
		double C_w = M_L / R / T;
		double dP_gwdPC = Deriv_dPgw_dPC(PC, T);
		drhoGw_dPC = C_w*dP_gwdPC;
		return drhoGw_dPC;
	}

	/*
	*water vapor saturation pressure
	*/
	virtual double get_P_sat_S(double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325.0;
		double h_wg = 2258000.0;
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//double P_sat1 = 611 * exp(17.27*T / (237.3 + T));
		//double P_sat1 = 1E-3*exp(19.84 - 4975 / T)*R*T / M_L;
		return P_sat;
	}
	
	/*
	*Kelvin equation
	*/
	virtual double get_P_G_w(double PC, double T)
	{
		double P_sat(0.0);
		double P_gw(0.0);
		P_sat = get_P_sat_S(T);
		double C_w = M_L / R / T;
		P_gw = P_sat*exp(-PC*C_w / rho_l_std);
		return P_gw;
	}
	virtual double Deriv_dPsat_dT(double T)
	{
		// Here unit of T is Celsius;
		double dPsat_dT(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		dPsat_dT = P_0*(M_L*h_wg / R)*(1./T/T)*exp(((1. / T_0) - (1. / T))*M_L*h_wg / R);
		return dPsat_dT;
	}
	virtual double Deriv_dPgw_dT(double PC, double T)
	{
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		double dPgw_dT(0.0);
		double C_w = M_L / R / T;
		double dPsat_dT(0.0);
		double P_sat(0);
		P_sat = get_P_sat_S(T);
		dPsat_dT = Deriv_dPsat_dT(T);
		dPgw_dT = dPsat_dT*exp(-PC*C_w / rho_l_std) + P_sat*exp(-PC*C_w / rho_l_std)*(PC*M_L / rho_l_std / R / T / T);
		//dPgw_dT = exp((-PC*M_L / (rho_l_std* R*T)) + (((1 / T_0) - (1 / T))*M_L*h_wg / R))*P_0*((h_wg*M_L / (R*T*T)) + (M_L*PC / (R*rho_l_std*T*T)));
		return dPgw_dT;
	}
	virtual double Deriv_dPgw_dPC(double PC, double T)
	{
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		double dPgw_dPC(0.0);
		double C_w = M_L / R / T;
		double P_sat(0.0);
		P_sat = get_P_sat_S(T);
		dPgw_dPC = P_sat*exp(-PC*C_w / rho_l_std)*(-C_w / rho_l_std);
		//dPgw_dPC = -exp((-PC*M_L / (rho_l_std* R*T)) + (((1 / T_0) - (1 / T))*M_L*h_wg / R)) * P_0*M_L / (R*rho_l_std*T);
		return dPgw_dPC;
	}

	/**
	* OVERALL HEAT CAPACITY
	*/
	virtual double get_overall_Heat_Capacity(double Sg, double lambda_pm_dry = 0.582, double lambda_pm_wet = 1.14)
	{
		double lambda_pm(0.0);
		lambda_pm = lambda_pm_dry + std::pow(1 - Sg, 0.5)*(lambda_pm_wet - lambda_pm_dry);
		if (Sg < 0)
			lambda_pm = lambda_pm_wet;
		else if (Sg > 1)
			lambda_pm = lambda_pm_dry;
		return lambda_pm;

	}


private:
	const double R = 8.314;
	const double S_gr = 0.0;
	const double S_lr = 0.15;
	const double P_r = 5.916e+5;
	const double Lambda_Leverett = 0.05878;
	const double C_h = 0.0;// Here I define a constant value C_h which is equal to Hen*M_G(Molar mass of gas)
	 //const double P_vapor =1.0 ;
	 const double rho_l_std = 1000.0;
	 const double M_G=0.02896;
	 const double M_L=0.018;
	double PG_L;
	double PC_L;
	double T_L;
};


