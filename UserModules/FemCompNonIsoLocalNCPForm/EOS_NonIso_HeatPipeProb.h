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
#include "AbstractEOS_NonIso_LocalNCPForm.h"
#include <cmath>


/**
*  EOS formulation for calculating each secondary variable
*  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg, with only 3 unknows! 
*/
class EOS_NonIso_HeatPipeProb : public AbstractEOS_NonIso_LocalNCPForm
{
public:

	/**
	  * Constructor
	  */
	EOS_NonIso_HeatPipeProb() : AbstractEOS_NonIso_LocalNCPForm(3)
	{
	};

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_NonIso_HeatPipeProb()
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
		P_L = env_conditon(0);//Liquid phase pressure
		// molar fraction of the light component
		X_L = env_conditon(1);// total mass density of water in each phase

		T_L = env_conditon(2);//temperature

	};

	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		// calculating PG_h partial pressure of light component
		double PGH = getPGH(P_L, T_L);//
		//calculate the residual
		res(0) =  X_L - rho_G_h*Sg - rho_L_h*(1 - Sg);
		res(1) = std::min(Sg, PGH*Hen*M_G - rho_L_h);
		res(2) = std::min(1 - Sg, rho_G_h - PGH*M_G / R / T_L);
	};

	virtual double NCP1(double Sg, double rho_L_h, double rho_G_h)
	{
		return X_L - rho_G_h*Sg - rho_L_h*(1 - Sg);
	}
	virtual double NCP2(double Sg, double rho_L_h, double rho_G_h)
	{
		double PGH = getPGH(P_L, T_L);//
		return std::min(Sg, PGH*Hen*M_G - rho_L_h);
	}
	virtual double NCP3(double Sg, double rho_L_h, double rho_G_h)
	{
		double PGH = getPGH(P_L, T_L);//
		return std::min(1 - Sg, rho_G_h - PGH*M_G / R / T_L);
	}

	/**
	* realization of the calc_Jacobian function.
	* this function will evaluate the Jacobian matrix,
	* based on the unknown values given.
	*/

	virtual void calc_Jacobian_fd(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		
		J.setZero();
		J(0, 0) = -rho_G_h + rho_L_h;//-----dF(1)/dSg
		J(0, 1) = -1 + Sg;//dRHO_L_H
		J(0, 2) = -Sg;
		J(1, 0) = (NCP2(Sg + EPS, rho_L_h, rho_G_h) - NCP2(Sg - EPS, rho_L_h, rho_G_h)) / 2 / EPS;
		J(1, 1) = (NCP2(Sg , rho_L_h + EPS, rho_G_h) - NCP2(Sg , rho_L_h - EPS, rho_G_h)) / 2 / EPS;
		J(1, 2) = (NCP2(Sg, rho_L_h , rho_G_h+ EPS) - NCP2(Sg, rho_L_h , rho_G_h- EPS)) / 2 / EPS;
		J(2, 0) = (NCP3(Sg + EPS, rho_L_h, rho_G_h) - NCP3(Sg - EPS, rho_L_h, rho_G_h)) / 2 / EPS;
		J(2, 1) = (NCP3(Sg, rho_L_h + EPS, rho_G_h) - NCP3(Sg, rho_L_h - EPS, rho_G_h)) / 2 / EPS;
		J(2, 2) = (NCP3(Sg, rho_L_h, rho_G_h + EPS) - NCP3(Sg, rho_L_h, rho_G_h - EPS)) / 2 / EPS;
	}
	/**
	* realization of the calc_Jacobian function.
	* this function will evaluate the Jacobian matrix,
	* based on the unknown values given.
	*/

	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		
		
		double PGH = getPGH(P_L,T_L);
		double F1 = Sg;
		double G1 = PGH*Hen*M_G - rho_L_h;
		double F2 = 1 - Sg;
		double G2 = rho_G_h - PGH*M_G / R / T_L;
		// evaluate J
		J.setZero();
		J(0, 0) = -rho_G_h+rho_L_h;//-----dF(1)/dSg
		J(0, 1) = -1 + Sg;//dRHO_L_H
		J(0, 2) = -Sg;
		if (F1 <= G1) {
			J(1, 0) = 1.0;
			J(1, 1) = 0.0;
			J(1, 2) = 0.0;
		}
		else{
			J(1, 0) = 0.0;// -C_v*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
			J(1, 1) = -1.0;
			J(1, 2) = 0.0;	
		}
		if (F2 <= G2) {
			J(2, 0) = -1.0;
			J(2, 1) = 0.0;
			J(2, 2) = 0.0;
		}
		else{
			J(2, 0) = 0.0;// -C_v*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
			J(2, 1) = 0.0;
			J(2, 2) = 1.0;

		}
	};
	


	virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var)
	{
		
	};
	virtual double Calc_Res_Sg(double Sg, double rho_L_h, double rho_G_h)
	{
		double Res_Sg(0.0);
		Res_Sg = X_L - (Sg*rho_G_h + (1 - Sg)*rho_L_h);
		return Res_Sg;
	}

	
	/**
	*  return mass fraction of air in gas phase
	*/
	virtual double get_massfraction(double x_G_a)
	{
		double massfractionX(0.0);
		massfractionX = x_G_a / (x_G_a + (1 - x_G_a)*M_L / M_G);
		return massfractionX;
	}
	//calculate \rho_G_w
	virtual double get_RHO_G_W(double sg, double T)
	{
		double rho_G_w(0.0);
		double P_G_w(0.0);
		P_G_w = get_P_sat(T);
		rho_G_w = M_L*P_G_w / R / T;// std::max(, X_L);
		return rho_G_w;
	}
	//calculate \rho_G_h
	virtual double get_RHO_G_H(double PG)
	{
		return P_L*M_G / R / T_L;
		//return std::max(P_L*M_G / R / T_L, X_L);
	}


	virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
	{
		double p0 = 6.324555320336758e+05;
		double Pc = p0*0.0588*((1.263*Sg - 2.120)*Sg + 1.417) * Sg;
		if (Sg >= 1.0) {
			double y = p0*0.0588*((1.263*1.0 - 2.120)*1.0 + 1.417)*1.0;
			double m = p0*0.0588*((3 * 1.263*1.0 - 2 * 2.120)*1.0 + 1.417);
			return (Sg - 1)*m + y;
		}
		else if (Sg <= 0.0) {
			double y = 0.0;
			double m = p0*0.0588*1.417;
			return Sg*m + y;
		}
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{

		double dPcdSg(0.0);
		double p0 = 6.324555320336758e+05;
		if (Sg > 1.0)
			Sg = 1.0;
		else if (Sg <= 0.0) {
			double dPcdSg = p0*0.0588*1.417;
			return dPcdSg;
		}
		dPcdSg = p0*0.0588*((3 * 1.263*Sg - 2 * 2.120)*Sg + 1.417);
		return dPcdSg;
	}

	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;

		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	virtual double getPGH(double PG, double T)
	{
		double P_sat = get_P_sat(T);
		double PGH(0.0);
		PGH = (Hen*M_L* PG - rho_l_std + std::sqrt(pow(Hen* M_L* PG + rho_l_std, 2) - 4 * Hen *M_L* P_sat* rho_l_std)) / (2 * Hen* M_L);
		if (pow(Hen* M_L* PG + rho_l_std, 2) < 4 * Hen *M_L* P_sat* rho_l_std){
			//WARN("Solving Partial Pressure Comes an Error!");
			PGH = (Hen*M_L* PG - rho_l_std ) / (2 * Hen* M_L);
		}
		
		return PGH;
	}

	virtual double getPGW(double PG, double PGH, double T)
	{
		double P_sat = get_P_sat(T);
		double PGW(0.0);
		PGW = P_sat*(rho_l_std / (rho_l_std + (M_L / M_G)*Hen*M_G*PGH));

		return PGW;
	}
	virtual double Deriv_PGH_dPG(double PG, double PGH, double T)
	{
		//double test=(getPGH(PG + EPS*PG, T) - getPGH(PG - EPS*PG, T)) / 2 / EPS/PG;
		double P_sat = get_P_sat(T);
		return 1 / (1 - (Hen * M_L * rho_l_std *P_sat) / pow((rho_l_std + Hen * M_L * PGH), 2));
	}

	virtual double Deriv_PGH_dT(double PG, double PGH, double T)
	{
		double P_sat = get_P_sat(T);
		double Derivative_Psat = Deriv_dPsat_dT(T);
		//double test= (getPGH(PG, T + EPS*T) - getPGH(PG, T - EPS*T)) / 2 / EPS/T;
		return -((rho_l_std* (rho_l_std + Hen * M_L * PGH) * Derivative_Psat) / (pow((rho_l_std + Hen * M_L * PGH), 2) - Hen * M_L *rho_l_std*  P_sat));
	}

	virtual double Deriv_dSgdX(double Sg, double rho_L_h, double rho_G_h, double drho_L_h_dX, double drho_G_hdX)
	{
		
		return 1 / (rho_G_h - rho_L_h);
	}
	
	virtual double Deriv_dSgdP(double Sg, double rho_L_h, double rho_G_h, double drho_L_h_dP, double drho_G_hdP)
	{

		return -(drho_L_h_dP*(1 - Sg) + drho_G_hdP*Sg) / (rho_G_h - rho_L_h);
	}
	
	virtual double Deriv_dSgdT(double Sg, double rho_L_h, double rho_G_h, double drho_L_h_dT, double drho_G_hdT)
	{

		return -(drho_L_h_dT*(1 - Sg) + drho_G_hdT*Sg) / (rho_G_h - rho_L_h);
	}

	/*
	*calculate the partial vapor pressure 
	* regulated by water saturation pressure 
	*/
	virtual double get_P_G_w(double PG, double rho_L_h, double PC, double T)
	{
		double P_sat(0.0);
		double P_gw(0.0);
		P_sat = get_P_sat_S(T);
		double C_w = M_L / R / T;
		P_gw = P_sat*exp(-(PC)*C_w / rho_l_std)*(rho_l_std / (rho_l_std + M_L*rho_L_h/M_G));
		return P_gw;
	}

	virtual double get_P_sat_S(double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325.0;
		double h_wg = 2258000.0;
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//if (T>T_0+2)
		//P_sat = P_0;
		return P_sat;
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
	/*
	* alternative  
	*/
	virtual double get_Vapor_Pressure()
	{
		// Here unit of T is Celsius;
		double P_vapor(0.0);
		double tmp(0.0);
		double T_tmp = T_L - 273;
		tmp = 2.786 + 0.031514*T_tmp - 1.2373e-4*pow(T_tmp, 2) + 4.2267e-7*pow(T_tmp, 3) - 8.1308e-10*pow(T_tmp, 4);
		P_vapor = pow(10, tmp);
		return P_vapor;
	}
	
	virtual double Deriv_dPsat_dT(double T)
	{
		// Here unit of T is Celsius;
		double dPsat_dT(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;

		dPsat_dT = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R)*(M_L*h_wg / R)*(1/T/T);
		return dPsat_dT;
	}
	
	virtual double Deriv_dPGw_dP(double rho_L_h,double drho_L_hdP, double T)
	{
		double A = M_L / M_G;
		double P_sat = get_P_sat(T);
		return P_sat*(-A*rho_l_std / pow((rho_l_std + A*rho_L_h), 2))*drho_L_hdP;
		
	}
	virtual double Deriv_dPGw_dT(double rho_L_h, double drho_L_hdT, double T)
	{
		double A = M_L / M_G;
		double P_sat = get_P_sat(T);
		double dPsat_dT = Deriv_dPsat_dT(T);
		return dPsat_dT*(rho_l_std / (rho_l_std + A*rho_L_h)) + P_sat*(-A*rho_l_std / pow((rho_l_std + A*rho_L_h), 2))*drho_L_hdT;

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
	 const double Hen = 1e-4;
	 const double EPS = 1e-8;
	double P_L;
	double X_L;
	double T_L;
};


