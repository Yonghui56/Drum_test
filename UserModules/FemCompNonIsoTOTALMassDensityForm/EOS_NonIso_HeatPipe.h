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
#include "AbstractEOS_NonIso_TotalDensityForm.h"
#include <cmath>


/**
*  EOS formulation for calculating each secondary variable
*  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg, with only 3 unknows! 
*/
class EOS_NonIso_HeatPipe : public AbstractEOS_NonIso_TotalDensityForm
{
public:

	/**
	  * Constructor
	  */
	EOS_NonIso_HeatPipe() : AbstractEOS_NonIso_TotalDensityForm(2)
	{
	};

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_NonIso_HeatPipe()
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
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);

		T_L = env_conditon(2);

	};

	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_G_h = vec_unknowns(1);

		// calculating residual
		res(0) = Calc_Res_Sg(Sg, rho_G_h);
		res(1) = Calc_Res_rho_G_h(Sg, rho_G_h);
	};

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
		double rho_G_h = vec_unknowns(1);
		double PG(0.0);
		PG = getPG(Sg);
		double C_v = M_G / R / T_L;
		double PG_h = Function_PG_h(PG);
		double F1 = 1 - Sg;
		double G1 = rho_G_h -  std::max(C_v*PG_h, X_L);//C_v*PG_h; //
		// evaluate J
		J.setZero();
		J(0, 0) =  - rho_G_h;//-----dF(1)/dSg
		J(0, 1) = -Sg;

		if (F1 <= G1) {
			J(1, 0) = -1.0;
			J(1, 1) = 0.0;
		}
		else{
			J(1, 0) = -C_v*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
			J(1, 1) = 1.0;
			/*if (C_v*PG_h >= X_L){
			
			}
			else{
				INFO("something happens");
				J(1, 0) = 0.0;
				J(1, 1) = 1.0;
			}
			*/
			

		}
	};
	virtual void calc_Jacobian_fd(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_G_h = vec_unknowns(1);
		double eps = 1e-12;
		// evaluate J
		J.setZero();
		J(0, 0) = (Calc_Res_Sg(Sg + eps,rho_G_h) - Calc_Res_Sg(Sg - eps, rho_G_h)) / eps / 2;
		J(0, 1) = (Calc_Res_Sg(Sg, rho_G_h + eps) - Calc_Res_Sg(Sg, rho_G_h - eps)) / eps / 2;
		J(1, 0) = (Calc_Res_rho_G_h(Sg + eps, rho_G_h) - Calc_Res_rho_G_h(Sg - eps, rho_G_h)) / eps / 2;
		J(1, 1) = (Calc_Res_rho_G_h(Sg, rho_G_h + eps) - Calc_Res_rho_G_h(Sg, rho_G_h - eps)) / eps / 2;
		
	};


	virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var)
	{
		
	};
	virtual double Calc_Res_Sg(double Sg, double rho_G_h)
	{
		double Res_Sg(0.0);
		Res_Sg = X_L - (Sg*rho_G_h);
		return Res_Sg;
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 2
	* rho_L_h=min(C_h*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_L_h(double Sg, double rho_L_h)
	{
		double Res_rho_L_h(0.0);
		Res_rho_L_h = 0.0;// std::min(Sg, -rho_L_h + get_RHO_L_H(Sg));//
		return Res_rho_L_h;
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 3
	* rho_G_h=MAX(C_V*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_G_h(double Sg, double rho_G_h)
	{
		double Res_rho_G_h(0.0);
		Res_rho_G_h = std::min((1 - Sg), rho_G_h - get_RHO_G_H(Sg));
		return Res_rho_G_h;
	}
	virtual double get_RHO_G_H(double sg)
	{
		double rho_G_h(0.0);
		double PG(0.0);
		PG = getPG(sg);
		double PG_H(0.0);
		PG_H = Function_PG_h(PG);
		rho_G_h = M_G*PG_H / R / T_L;// std::max(, X_L);
		return rho_G_h;
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

	virtual double getPG(double Sg)
	{
		double PG = 0.0;
		double PC = 0.0;
		PC = getPcbySg(Sg);
		PG = P_L + PC;// (1 - gamma)*PC;
		return PG;
	}

	virtual double Deriv_dPGH_dPG(double Sg)
	{
		double PG_h(0.0);
		double PG(0.0);
		double dPGH_dPG(0.0);
		PG = getPG(Sg);
		PG_h = Function_PG_h(PG);
		double P_vapor(0.0);
		P_vapor = get_P_sat(T_L);
		double A = 0.0;
		A = M_L / M_G;
		double B = 0.0;
		//B = P_vapor*rho_l_std*A*C_h / (rho_l_std + A*C_h*PG_h);
		dPGH_dPG = 1 / (1 - B);
		return dPGH_dPG;
	}

	virtual double Function_Alpha(double Sg)
	{
		double dPGH_dPG(0.0);
		dPGH_dPG = Deriv_dPGH_dPG(Sg);
		double Alpha(0.0);
		double C_v = M_G / R / T_L;
		Alpha = (C_v*Sg )*dPGH_dPG;
		return Alpha;
	}
	virtual double Function_Beta(double Sg)
	{
		double PGh(0.0);
		double PG(0.0);
		double Beta(0.0);
		PG = getPG(Sg);
		PGh = Function_PG_h(PG);
		double C_v = M_G / R / T_L;
		Beta = C_v*PGh;
		return Beta;
	}
	virtual double Deriv_dPGdSg(double Sg, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double dPcdSg(0.0);
		dPcdSg = Deriv_dPCdS(Sg);
		double PC(0.0);
		PC = getPcbySg(Sg);
		double Gamma(0.0);
		//Gamma = getweightedFunc(Sg);
		double dPGdSg(0.0);
		dPGdSg = dPcdSg*(1 - Gamma);// -PC*get_deriv_weight_Func(Sg);
		//if (Sg < 1e-14)
		//dPGdSg = 0.0;
		return dPGdSg;
	}

	virtual double Deriv_dSgdP(double Sg)
	{

		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double dSgdP(0.0);
		dSgdP = -alpha / (beta + alpha*dPG_dSg);
		return dSgdP;
	}

	virtual double Deriv_dSgdX(double Sg)
	{
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double dSgdX(0.0);
		dSgdX = 1 / (beta + alpha*dPG_dSg);
		return dSgdX;
	}

	virtual double Deriv_dSgdT(double Sg,double rho_G_h, double drho_G_hdT)
	{
		
		double dSgdT(0.0);
		dSgdT = -drho_G_hdT*Sg / rho_G_h;
		return dSgdT;
	}


	virtual double get_RHO_G_W(double T)
	{
		double rho_G_w(0.0);
		
		double PG_w = get_P_sat(T);
		double C_w = M_L / R / T_L;
		rho_G_w = PG_w*C_w;
		return rho_G_w;
	}	
	
	virtual double Function_PG_h(double PG)
	{
		double PG_H(0.0);
		double P_vapor(0.0);
		P_vapor = get_P_sat(T_L);
		double C_h = 0.0;
		PG_H = PG - P_vapor;
		//if (PG_H < 0)
			//PG_H = 0.0;
		//PG_H = (1 / (2 * C_h*M_L))* (C_h*M_L*PG - M_G* rho_l_std + sqrt(pow((-C_h*M_L*PG + M_G*rho_l_std), 2) - 4 * C_h*M_L*(-M_G*PG*rho_l_std + M_G* P_vapor* rho_l_std)));
		return PG_H;

	}

	/**
	* The intermedia function
	* Omega=beta(p,Sg)*Func_characteristic(P,X)/(beta(P,Sg)+alpha(P,Sg)*dPgdSg(P,Sg))
	* M=Func_characteristic(P,X)*dPgdSg(P,Sg)/(beta(P,Sg)+alpha(P,Sg)*dPgdSg(P,Sg))
	*/
	virtual double Func_Omega(double Sg)
	{
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double Omega(0.0);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		Omega = beta / (beta + alpha*dPG_dSg);
		return Omega;

	}
	virtual double Func_InterM(double Sg)
	{
		double M(0.0);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		M = 1 / (beta / dPG_dSg + alpha);
		return M;

	}
	/**
	* Deriv rho_L_h based on mean pressure
	*/
	virtual double Deriv_drhoLh_dP(double Sg, double Omega)
	{
		double drho_L_h_dP(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_L_h_dP = C_h*dPG_H_dPG*Omega;
		return drho_L_h_dP;
	}

	/**
	* Deriv rho_L_h based on total mass density
	*/
	virtual double Deriv_drhoLh_dX(double Sg, double M)
	{
		double drho_L_h_dX(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_L_h_dX = C_h*dPG_H_dPG*M;
		return drho_L_h_dX;
	}

	/**
	* Deriv rho_G_h based on P mean pressure
	*/
	virtual double Deriv_drhoGh_dP(double Sg, double Omega)
	{
		double drho_G_h_dP(0.0);
		double dPG_H_dPG(0.0);
		double C_v = M_G / R / T_L;
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_G_h_dP = C_v*dPG_H_dPG*Omega;
		return drho_G_h_dP;
	}

	/**
	* Deriv rho_G_h based on X --total mass density
	*/
	virtual double Deriv_drhoGh_dX(double Sg, double M)
	{
		double drho_G_h_dX(0.0);
		double dPG_H_dPG(0.0);
		double C_v = M_G / R / T_L;
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_G_h_dX = C_v*dPG_H_dPG*M;
		return drho_G_h_dX;
	}

	virtual double Deriv_drhoGh_dT(double PG_h, double T)
	{
		
		double C_v = M_G / R / T;
		double dPG_hdT = -Deriv_dPsat_dT(T);
		return dPG_hdT*C_v - PG_h*C_v / T;
		
	}

	virtual double Deriv_drhoGW_dT(double PG_w, double T)
	{
		double drhoGW_dT(0.0);
		double C_w = M_L / R / T;
		double dPG_wdT = Deriv_dPsat_dT(T);
		drhoGW_dT = dPG_wdT*C_w - PG_w*C_w / T;
		return drhoGW_dT;
	}

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
	virtual double get_P_sat(double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//P_sat *= exp(-(PC - 0.7 * 101325)*M_L / R / T / 1000);
		return P_sat;
	}
	virtual double Deriv_dPsat_dT(double T)
	{
		// Here unit of T is Celsius;
		double dPsat_dT(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;

		dPsat_dT = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R)*(M_L*h_wg / R)*(1/T/T);
		//dPsat_dT = P_0*(M_L*h_wg / R)*(1. / T / T)*exp(((1. / T_0) - (1. / T))*M_L*h_wg / R);
		return dPsat_dT;
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
	const double P_r = 6.324555320336758e+05;
	const double Lambda_Leverett = 0.05878;
	 const double rho_l_std = 1000.0;
	 const double M_G=0.02896;
	 const double M_L=0.018;
	 const double C_h = 0.0;
	double P_L;
	double X_L;
	double T_L;
};


