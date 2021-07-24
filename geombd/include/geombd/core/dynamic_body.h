/*
 * MIT License
 *
 * Copyright (c) 2020 Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Carla Villanueva <carla.villanueva@cinvestav.edu.mx>,
 * Alvaro Paz <alvaro.paz@cinvestav.edu.mx>
 * CINVESTAV - Saltillo Campus
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions
 * of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
 * TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 *	\file include/geombd/core/dynamic_body.h
 *	\author Gustavo Arechavaleta, Carla Villanueva, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 *  DynamicBody class. It has functions to build a body from a xml file. It is a derived class from body.
 */

#ifndef GEOMBD_CORE_DYNAMIC_BODY_H
#define GEOMBD_CORE_DYNAMIC_BODY_H

#include "body.h"
#include "types.h"

namespace geo{

//! DynamicBody class, derived from Body class.
/*! \ingroup core_module
    \sa Body
  */
class DynamicBody : public Body
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:

    //! Custom constructor. Initializes the body with the given parameters.
    /*! \param id The id of the body.*/
    /*! \param name The name of the body.*/
    DynamicBody(const unsigned int id, const std::string &name);

    //! Custom constructor. Initializes the body with the given parameters.
    /*! \param id The id of the body.*/
    /*! \param name The name of the body.*/
    /*! \param pos The position of the body.*/
    /*! \param rot The rotation matrix (Rx, Ry, Rz).*/
    DynamicBody(const unsigned int id, const std::string &name,
                const Vector3r &pos,
                const Matrix3r &rot);

    //! Default destructor.
    ~DynamicBody();

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:

    //! Mass of the body.
    real_t mass;

    //! 3D position of the body CoM w.r.t the body reference frame.
    Vector3r centerOfMass;

    //! body inertia tensor w.r.t. the body reference frame.
    Matrix3r inertiaTensor;

    //! body configuration (parameterized by a homogeneous matrix) at home pose w.r.t. the body reference frame.
    Matrix4r homeConfig;

    //! body global configuration (parameterized by a homogeneous matrix) w.r.t. the inertial reference frame.
    Matrix4r globalConfiguration;

    //! spatial inertia matrix J w.r.t. the body reference frame.
    SpatialMatrix spatialInertia;

    //! OperationalHandle linear and angular velocity vector (vx,vy,vz,omegax,omegawy,omegaz) relative to the associated body.
    SpatialVector twist;

    //! OperationalHandle linear and angular acceleration vector (dvx,dvy,dvz,dwx,dwy,dwz) relative to the associated body.
    SpatialVector acceleration;

    //! OperationalHandle linear and angular acceleration vector (dvx,dvy,dvz,dwx,dwy,dwz) relative to the associated body.
    SpatialVector quasiAcceleration;

    //! OperationalHandle force and moments vector (fx,fy,fz,nx,ny,nz) relative to the associated body.
    SpatialVector wrench;

    //! Spatial screw axes relative to the associated body.
    SpatialVector screwAxes;

    //! Adjoint_{i-1}^{i} operator.
    SpatialMatrix Adjoint;

    //! AdjointDual_{i-1}^{i} operator.
    SpatialMatrix AdjointDual;

    //! adjoint_{i-1}^{i} operator.
    SpatialMatrix adjoint;

    //! adjointDual_{i-1}^{i} operator.
    SpatialMatrix adjointDual;

    //! Differentiation of Twist vector
    D_SpatialVector Diff_Twist;

    //! Differentiation of Acceleration vector
    D_SpatialVector Diff_Acceleration;

    //! Differentiation of Wrench vector
    D_SpatialVector Diff_Wrench;

    //! adjoint of ScrewAxes
    SpatialMatrix adjointScrew;

    //! adjoint dual of -ScrewAxes
    SpatialMatrix adjointDualScrew;

    //! cummulative spatial inertia matrix
    SpatialMatrix In;

    //! spatial acceleration bias
    SpatialVector c_bias;

    //! cummulative wrench bias
    SpatialVector Pn;

    //! cummulative wrench bias
    SpatialVector Pa;

    //! U <- I*S
    SpatialVector Ufd;

    //! inverse of inertial factor
    real_t invD;

    //! control variable
    real_t ufd;

    //! current spatial inertia matrix for back-projection
    SpatialMatrix Ia;

    //! Differentiation of Ufd vector
    D_SpatialVector Diff_Ufd;

    //! Differentiation of inverse of inertial vector
    D_real_t Diff_invD;

    //! Differentiation of control variable
    D_real_t Diff_ufd;



    /// Under development (Enhanced diff)--------------------

public:

    SpatialVector UinvD;

    //! body configuration (parameterized by a homogeneous matrix) w.r.t. the body reference frame.
    Matrix4r configuration;

    //! Dense terms of differentiated Twist vector
    D_SpatialVector D_q_V, D_dq_V;

    //! Dense terms of differentiated spatial acceleration bias
    D_SpatialVector D_q_c, D_dq_c;

    //! Dense terms of differentiated wrench vector bias
    D_SpatialVector D_q_p, D_dq_p;

    //! Dense term of differentiated vector U (horizontally arranged)
    D_SpatialVector D_U_h;

    //! Dense term of differentiated vector U (vertically arranged)
    VectorXr D_U_v;

    //! Dense term of differentiated D^{-1}
    RowVectorXr D_invD;

    //! Dense terms of differentiated scalar u
    RowVectorXr D_q_u, D_dq_u;

    //! Dense term of differentiated accumulated inertia matrix
    D_SpatialMatrix D_Ia;

    //! Dense terms of differentiated accumulated wrench bias
    D_SpatialVector D_q_Pa, D_dq_Pa;

    //! Dense term of differentiated inertia matrix back-projection
    D_SpatialMatrix D_In;

    //! Dense terms of differentiated wrench bias back-projection
    D_SpatialVector D_q_Pn, D_dq_Pn;

    //! Dense terms of differentiated pseudo spatial acceleration
    D_SpatialVector D_q_dVa, D_dq_dVa;

    //! Dense terms of differentiated joint acceleration
    RowVectorXr D_q_ddq, D_dq_ddq;

    //! Dense terms of differentiated spatial acceleration
    D_SpatialVector D_q_dV, D_dq_dV;

    //! Dense term of D_Ia*c (vertically arranged)
    VectorXr D_IaC_v;

    //! Dense term of D_Ia*c (horizontally arranged)
    D_SpatialVector D_IaC_h;



    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:

    //! Load the geometry of the body.
    /*! \param XML file name with the the geometry
    * \param body's relative position
    * \return status
             */
    bool loadGeometry( const std::string &xmlFileName, Vector3r &relativeBodyPosition);

    //! Load the geometry of the body.
    /*! \param XML file name with the the geometry
    * \return status
             */
    bool loadGeometryURDF( const std::string &xmlFileName);

    //! Set the mass of the body.
    /*! \param a real_t precision number
        \return void
             */
    void setMass( real_t mass ){ this->mass = mass; }

    //! Set the center of mass of the body.
    /*! \param a 3-dimensional real_t precision array.
        \return void
             */
    void setCoM( const Vector3r &centerOfMass );

    //! Set the inertia Tensor of the body
    /*! \param a real_t precision matrix (3x3 matrix)
        \return void
             */
    void setInertiaTensor( const Matrix3r &inertiaTensor );

    //! Set the body configuration
    /*! \param a real_t precision matrix (4x4 matrix)
        \return void
             */
    void setConfiguration( const Matrix4r &configuration );

    //! Set the body global configuration
    /*! \param a real_t precision matrix (4x4 matrix)
        \return void
             */
    void setGlobalConfiguration( const Matrix4r &globalConfiguration );

    //! Set the body home configuration
    /*! \param a real_t precision matrix (4x4 matrix)
        \return void
             */
    void setHomeConfig( const Matrix4r &homeConfig );

    //! Set the spatial inertia of the body
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setSpatialInertia( const SpatialMatrix &spatialInertia );

    //! Set the linear and angular velocity at opHandle
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setTwist( const SpatialVector &twist );

    //! Set the linear and angular acceleration at opHandle
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setAcceleration( const SpatialVector &acceleration );

    //! Set the linear and angular acceleration at opHandle
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setQuasiAcceleration( const SpatialVector &quasiAcceleration );

    //! Set the force and moment at opHandle
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setWrench( const SpatialVector &wrench );

    //! Set the screw axes
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setScrewAxes( const SpatialVector &screwAxes );

    //! Set the Adjoint_{i-1}^{i} operator
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setAdjoint( const SpatialMatrix &Adjoint );

    //! Set the AdjointDual_{i-1}^{i} operator
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setAdjointDual( const SpatialMatrix &AdjointDual );

    //! Set the adjoint_{i-1}^{i} operator
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setadjoint( const SpatialMatrix &adjoint );

    //! Set the adjointDual_{i-1}^{i} operator
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setadjointDual( const SpatialMatrix &adjointDual );

    //! Set the Differentiation of the Twist vector
    /*! \param a real_t precision matrix (6xX matrix)
        \return void
             */
    void setDiff_Twist( const D_SpatialVector &Diff_Twist );

    //! Set the Differentiation of the Acceleration vector
    /*! \param a real_t precision matrix (6xX matrix)
        \return void
             */
    void setDiff_Acceleration( const D_SpatialVector &Diff_Acceleration );

    //! Set the Differentiation of the Wrench vector
    /*! \param a real_t precision matrix (6xX matrix)
        \return void
             */
    void setDiff_Wrench( const D_SpatialVector &Diff_Wrench );

    //! Set the adjoint operator of ScrewAxes
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setadjointScrew( const SpatialMatrix &adjointScrew );

    //! Set the adjoint dual operator of -ScrewAxes
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setadjointDualScrew( const SpatialMatrix &adjointDualScrew );

    //! Set the cummulative spatial inertia of the body
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setIn( const SpatialMatrix &In );

    //! Set the spatial acceleration bias
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setC_bias( const SpatialVector &c_bias );

    //! Set the wrench bias
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setPn( const SpatialVector &Pn );

    //! Set the wrench bias
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setPa( const SpatialVector &Pa );

    //! Set inertial current variable
    /*! \param a 6D point of real_t precision numbers.
        \return void
             */
    void setUfd( const SpatialVector &Ufd );

    //! Set the inertial factor
    /*! \param a real_t precision number
        \return void
             */
    void setInvD( real_t invD ){ this->invD = invD; }

    //! Set the control factor
    /*! \param a real_t precision number
        \return void
             */
    void setufd( real_t ufd ){ this->ufd = ufd; }

    //! Set the current spatial inertia for back-projection
    /*! \param a real_t precision matrix (6x6 matrix)
        \return void
             */
    void setIa( const SpatialMatrix &Ia );

    //! Set the Differentiation of the Twist vector
    /*! \param a real_t precision matrix (6xX matrix)
        \return void
             */
    void setDiff_Ufd( const D_SpatialVector &Diff_Ufd );

    //! Set the Differentiation of the inverse of D
    /*! \param a real_t precision matrix (1xX matrix)
        \return void
             */
    void setDiff_invD( const D_real_t &Diff_invD );

    //! Set the Differentiation of the inertial factor u
    /*! \param a real_t precision matrix (1xX matrix)
        \return void
             */
    void setDiff_ufd( const D_real_t &Diff_ufd );


    //! Set the dense terms of differentiated Twist vector
    /*! \param two real_t spatial vector Jacobian
        \return void
             */
    void setDiff_x_V( const D_SpatialVector &D_q_V, const D_SpatialVector &D_dq_V );


    //! Get the mass of the body.
    /*! \param none
        \return a real_t precision number
             */
    real_t getMass(){ return mass; }

    //! Get the center of mass of the body.
    /*! \param none
        \return a 6-dimensional real_t precision array (3D position, orientation)
             */
    Vector3r getCoM(){ return centerOfMass; }

    //! Get the inertia tensor of the body.
    /*! \param none
        \return a real_t precision matrix (3x3 matrix)
             */
    Matrix3r getInertiaTensor(){ return inertiaTensor; }

    //! Get the body configuration.
    /*! \param none
        \return a real_t precision matrix (4x4 matrix)
             */
    Matrix4r getConfiguration(){ return configuration; }

    //! Get the body global configuration.
    /*! \param none
        \return a real_t precision matrix (4x4 matrix)
             */
    Matrix4r getGlobalConfiguration(){ return globalConfiguration; }

    //! Get the body home configuration.
    /*! \param none
        \return a real_t precision matrix (4x4 matrix)
             */
    Matrix4r getHomeConfig(){ return homeConfig; }

    //! Get the inertia tensor of the body.
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getSpatialInertia() { return spatialInertia; }

    //! Get the twist of the OperationalHandle.
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular velocity)
             */
    SpatialVector getTwist() { return twist; }

    //! Get the acceleration of the OperationalHandle.
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular acceleration)
             */
    SpatialVector getAcceleration() { return acceleration; }

    //! Get the acceleration of the OperationalHandle.
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular acceleration)
             */
    SpatialVector getQuasiAcceleration() { return quasiAcceleration; }

    //! Get the wrench of the OperationalHandle.
    /*! \param none
        \return a 6-dimensional real_t precision array (force and moment)
             */
    SpatialVector getWrench() { return wrench; }

    //! Get the screw axes spatial vector.
    /*! \param none
        \return a 6-dimensional real_t precision array
             */
    SpatialVector getScrewAxes() { return screwAxes; }

    //! Get the Adjoint_{i-1}^{i} operator
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getAdjoint() { return Adjoint; }

    //! Get the AdjointDual_{i-1}^{i} operator
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getAdjointDual() { return AdjointDual; }

    //! Get the adjoint_{i-1}^{i} operator
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getadjoint() { return adjoint; }

    //! Get the adjointDual_{i-1}^{i} operator
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getadjointDual() { return adjointDual; }

    //! Get the Differentiation of the Twist spatial vector
    /*! \param none
        \return a real_t precision matrix (6xX matrix)
             */
    D_SpatialVector getDiff_Twist() { return Diff_Twist; }

    //! Get the Differentiation of the Acceleration spatial vector
    /*! \param none
        \return a real_t precision matrix (6xX matrix)
             */
    D_SpatialVector getDiff_Acceleration() { return Diff_Acceleration; }

    //! Get the Differentiation of the Wrench spatial vector
    /*! \param none
        \return a real_t precision matrix (6xX matrix)
             */
    D_SpatialVector getDiff_Wrench() { return Diff_Wrench; }

    //! Get the adjoint element of ScrewAxes
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getadjointScrew() { return adjointScrew; }

    //! Get the adjoint dual element of -ScrewAxes
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getadjointDualScrew() { return adjointDualScrew; }

    //! Get the cummulative spatial inertia
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getIn() { return In; }

    //! Get the spatial acceleration bias
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular velocity)
             */
    SpatialVector getC_bias() { return c_bias; }

    //! Get the cummulative wrench bias
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular velocity)
             */
    SpatialVector getPn() { return Pn; }

    //! Get the cummulative wrench bias
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular velocity)
             */
    SpatialVector getPa() { return Pa; }

    //! Get the spatial inertia factor
    /*! \param none
         \return a 6-dimensional real_t precision array (linear and angular velocity)
             */
    SpatialVector getUfd() { return Ufd; }

    //! Get the inertia factor
    /*! \param none
        \return a real_t precision number
             */
    real_t getInvD(){ return invD; }

    //! Get the control factor
    /*! \param none
        \return a real_t precision number
             */
    real_t getufd(){ return ufd; }

    //! Get the current spatial inertia
    /*! \param none
        \return a real_t precision matrix (6x6 matrix)
             */
    SpatialMatrix getIa() { return Ia; }

    //! Get the Differentiation of the U vector
    /*! \param none
        \return a real_t precision matrix (6xX matrix)
             */
    D_SpatialVector getDiff_Ufd() { return Diff_Ufd; }

    //! Get the Differentiation of the invD factor
    /*! \param none
        \return a real_t precision matrix (1xX matrix)
             */
    D_real_t getDiff_invD() { return Diff_invD; }

    //! Get the Differentiation of the inertial factor u
    /*! \param none
        \return a real_t precision matrix (1xX matrix)
             */
    D_real_t getDiff_ufd() { return Diff_ufd; }

};
} // end of namespace geo
#endif // GEOMBD_CORE_DYNAMIC_BODY_H
