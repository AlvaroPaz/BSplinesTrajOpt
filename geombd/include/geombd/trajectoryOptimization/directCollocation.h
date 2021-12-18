/*
 * MIT License
 *
 * Copyright (c) 2020 Alvaro Paz <alvaro.paz@cinvestav.edu.mx>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx>
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
 *	\file include/geombd/trajectoryOptimization/directCollocation.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Class to implement the Trajectory Optimization -> Direct Collocation
 */

#ifndef GEOMBD_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H
#define GEOMBD_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H

#include <memory>
#include "geombd/core.h"
#include "geombd/dynamics.h"

namespace geo{

//! The type of the differentiation
/*! \ingroup core_module  */
enum DifferentiationType
{
    wrt_time,            /*!< Differentiate wrt the time. */
    wrt_state,           /*!< Differentiate wrt the robot state {q,dq}. */
    wrt_controlPoints    /*!< Differentiate wrt the control points. */
};


//! The type of constraint
/*! \ingroup core_module  */
enum ConstraintType
{
    constraint_initialConfiguration,                    /*!< constraint in q_i */
    constraint_finalConfiguration,                      /*!< constraint in q_f */
    constraint_initialGeneralizedVelocity,              /*!< constraint in dq_i */
    constraint_finalGeneralizedVelocity,                /*!< constraint in dq_f */
    constraint_pelvisSymmetry,                          /*!< constraint in pelvis for Nao robot */
    constraint_centerOfMass,                            /*!< constraint in CoM */
    constraint_centroidalMomentum                       /*!< constraint in centroidal Mu */
};


//! Defines a new eigen-based-matrix std::vector container type that is based on the ConstraintType.
/*!
  Defines a std::vector container of ConstraintType
 */
typedef std::vector< ConstraintType > ConstraintsStack;


struct robotSettingsTrajectoryOptimization {
    //! General Parameters
    short int n;                            // degrees of freedom
    short int numberControlPoints;          // number of control points
    int numberPartitions;                   // number of partitions
    int numberConstraints;                  // number of constraints
    real_t si, sf;                          // initial and final time
    VectorXr S;                             // vector time parameter
    VectorXr weights;                       // weights vector for cost function
    int numberFinalInterpolation;           // number of partitions of the final interpolation

    //! Differentiate with respect to
    DifferentiationType DifferentiationWRT;

    //! Boundaries
    VectorXr initialConfiguration;
    VectorXr finalConfiguration;
    VectorXr initialGeneralizedVelocity;
    VectorXr finalGeneralizedVelocity;

    //! Upper and lower tolerances for baundaries
    VectorXr qiBound_u;
    VectorXr qfBound_u;
    Vector2r comBound_u;
    SpatialVector muBound_u;

    VectorXr qiBound_l;
    VectorXr qfBound_l;
    Vector2r comBound_l;
    SpatialVector muBound_l;

    //! Constraints customization
    ConstraintsStack StackConstraints;
};


class DirectCollocation
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor
    DirectCollocation(std::shared_ptr< MultiBody > robot, std::shared_ptr< robotSettingsTrajectoryOptimization > robotSettings) : robotMotion{ robot } {

        this->robot = robot;
        this->n = robot->getDoF();

        this->f = 0;

        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->ddq = VectorXr::Zero(n,1);
        this->Tau = VectorXr::Zero(n,1);

        this->robotSettings = robotSettings;
        updateSettings();

    }


    ~DirectCollocation(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    protected:

    //! Dof
    short int n;

    //! Differentiation size
    int m;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;

    //! Robot dynamics object
    InverseDynamics robotMotion;

    //! Settings structure for NLP
    std::shared_ptr< robotSettingsTrajectoryOptimization > robotSettings;

    //! Stack of constraints
    ConstraintsStack StackConstraints;

    //! Index for storing constraints
    std::vector< Eigen::Vector2i > consIndex;

    //! Inner element of consIndex
    Eigen::Vector2i cIn;

    //! Time vector
    MatrixXr S;

    //! Cost function
    real_t f;

    //! Gradient of cost function
    MatrixXr g;

    //! Hessian of cost function
    MatrixXr H;

    //! number of constraints
    int numberConstraints;

    //! Constraints
    VectorXr constraintsValue;

    //! Jacobian of constraints
    MatrixXr JacobianConstraints;

    //! Jacobian of constraints
    MatrixXr HessianConstraints;

    //! basis function
    MatrixXr B;

    //! first derivative of B wrt s
    MatrixXr dB;

    //! second derivative of B wrt s
    MatrixXr ddB;

    //! Kronecker product of D_q and D_q
    MatrixXr DD_q;

    //! initial configuration
    VectorXr q_initial;

    //! final configuration
    VectorXr q_final;

    //! initial generalized velocity
    VectorXr dq_initial;

    //! final generalized velocity
    VectorXr dq_final;

    //! Generalized vectors
    VectorXr q, dq, ddq, Tau;

    //! Variables for differentiation
    MatrixXr Diff_Tau, DDiff_Tau;
    MatrixXr D_X, D_q, D_dq, D_ddq;


    public:



    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //! Update settings
         /*! \param void
         * \return void
         */
    void updateSettings(  );

    /*! Description: Builds the basis functions and its first two derivatives, i.e. B, dB, ddB
         * \parameters: Vector time S and the numer of control points per DoF N
         * \return: int 0
         */
    int buildBasisFunctions(  );    

    /*! Description: Builds the basis functions and its first two derivatives, i.e. B, dB, ddB
         * \parameters: Vector time S and the numer of control points per DoF N
         * \return: int 0
         */
    int buildBasisFunctions( short int N, VectorXr S );

    //! Compute the objective function
         /*! \param control points c and boolean flags for computing partial derivatives
         * \return void
         */
    void computeObjectiveFunction(const MatrixXr &c, const VectorXr &weights, const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute constraints of NLP
         /*! \param control points c and boolean flag for computing partials
         * \return void
         */
    void computeConstraints(const MatrixXr &c, const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Set initial configuration
         /*! \param initial configuration
         * \return void
         */
    void setInitialConfiguration( const VectorXr &q_initial ){ this->q_initial = q_initial; }

    //! Set final configuration
         /*! \param final configuration
         * \return void
         */
    void setFinalConfiguration( const VectorXr &q_final ){ this->q_final = q_final; }

    //! Set initial generalized velocity
         /*! \param initial generalized velocity
         * \return void
         */
    void setInitialGeneralizedVelocity( const VectorXr &dq_initial ){ this->dq_initial = dq_initial; }

    //! Set final generalized velocity
         /*! \param final generalized velocity
         * \return void
         */
    void setFinalGeneralizedVelocity( const VectorXr &dq_final ){ this->dq_final = dq_final; }

    //! Get the cost function value
    /*! \param none
        \return a real_t precision number
             */
    real_t getCost(){ return f; }

    //! Get the cost function gradient
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getCostGradient(){ return g; }

    //! Get the cost function Hessian
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getCostHessian(){ return H; }

    //! Get the constraints value
    /*! \param none
        \return a real vector
             */
    VectorXr getConstraints(){ return constraintsValue; }

    //! Get the constraints Jacobian
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getConstraintsJacobian(){ return JacobianConstraints; }

    //! Get the constraints Hessian
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getConstraintsHessian(){ return HessianConstraints; }

    //! Get the basis function
    /*! \param none
        \return Matrix real_t
             */
    MatrixXr getBasis(){ return B; }

    //! Get the basis function first derivative
    /*! \param none
         \return Matrix real_t
              */
    MatrixXr getDBasis(){ return dB; }

    //! Get the basis function second derivative
    /*! \param none
          \return Matrix real_t
               */
    MatrixXr getDDBasis(){ return ddB; }


protected:


};
} // end of namespace geo

#endif // GEOMBD_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H
