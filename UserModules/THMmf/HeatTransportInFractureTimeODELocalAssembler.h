/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file HeatTransportInFractureTimeODELocalAssembler.h
 *
 * Created on 2012-10-23 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/Fluid.h"

#include "Ogs6FemData.h"

/**
 * \brief Local assembly of time ODE components for heat transport in fracture
 *
 * \tparam T_TIME_DIS   Time discretization scheme
 */
template <class T_TIME_DIS>
class HeatTransportInFractureTimeODELocalAssembler: public T_TIME_DIS
{
public:

    /**
     *
     * @param feObjects
     */
    HeatTransportInFractureTimeODELocalAssembler(const FemLib::IFeObjectContainer* feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
    : _feObjects(feObjects->clone()), _problem_coordinates(problem_coordinates), _vel(NULL)
    {
    };

    /**
     *
     */
    virtual ~HeatTransportInFractureTimeODELocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    };

    /**
     *
     * @param vel
     */
    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

protected:
    /**
     * assemble components of local time ODE
     *
     * @param time
     * @param e
     * @param u1
     * @param u0
     * @param localM
     * @param localK
     * @param localF
     */
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localM, MathLib::LocalMatrix &localK, MathLib::LocalVector &/*localF*/);

private:
    FemLib::IFeObjectContainer* _feObjects;
    const MeshLib::CoordinateSystem _problem_coordinates;
    NumLib::ITXFunction* _vel;
};

template <class T>
void HeatTransportInFractureTimeODELocalAssembler<T>::assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localM, MathLib::LocalMatrix &localK, MathLib::LocalVector &/*localF*/)
{
    // element information
    const size_t mat_id = e.getGroupID();
    const size_t n_dim = _problem_coordinates.getDimension();
    MeshLib::ElementCoordinatesMappingLocal* ele_local_coord;
    ele_local_coord = (MeshLib::ElementCoordinatesMappingLocal*)e.getMappedCoordinates();
    assert(ele_local_coord!=NULL);
    //const MathLib::LocalMatrix* matR = &ele_local_coord->getRotationMatrixToOriginal();

    // material data
    MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
    MaterialLib::Solid* solid = Ogs6FemData::getInstance()->list_solid[mat_id];
    MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];

    // numerical integration
    MathLib::LocalMatrix localDispersion = MathLib::LocalMatrix::Zero(localK.rows(), localK.cols());
    MathLib::LocalMatrix localAdvection = MathLib::LocalMatrix::Zero(localK.rows(), localK.cols());
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        // compute shape functions and parameters for integration
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        MathLib::LocalMatrix &Np = *fe->getBasisFunction();
        MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();
        const double fac = fe->getDetJ() * q->getWeight(j);

        // evaluate material properties
        fe->getRealCoordinates(real_x);
        NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
        //fluid
        double rho_f = .0;
        double cp_f = .0;
        double lambda_f = .0;
        fluid->density->eval(gp_pos, rho_f);
        fluid->specific_heat->eval(gp_pos, cp_f);
        fluid->thermal_conductivity->eval(gp_pos, lambda_f);
        //solid
        double rho_s = .0;
        double cp_s = .0;
        double lambda_s = .0;
        solid->density->eval(gp_pos, rho_s);
        solid->specific_heat->eval(gp_pos, cp_s);
        solid->thermal_conductivity->eval(gp_pos, lambda_s);
        //porous media
        double b = .0;
        double n = .0;
        double rho_cp_eff = n * rho_f * cp_f + (1.-n) * rho_s * cp_s;
        double lambda_eff = n * lambda_f + (1.-n) * lambda_s;
        pm->geo_area->eval(gp_pos, b);
        pm->porosity->eval(gp_pos, n);

        NumLib::ITXFunction::DataType v;
        _vel->eval(gp_pos, v);
        NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

        // M += Np^T * b * rho_cp_eff * Np
        localM.noalias() += fac * Np.transpose() * b* rho_cp_eff * Np;

        // Diff += dNp^T * b * lambda * dNp
        localDispersion.noalias() += fac * dNp.transpose() * b * lambda_eff * dNp;

        // Adv += Np^T * b*rho*cp*v * dNp
        localAdvection.noalias() += fac * Np.transpose() * b * rho_f * cp_f * v2 * dNp;
    }

    localK = localDispersion + localAdvection;

    //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
    //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
    //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
}
