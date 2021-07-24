/*
 *
 * Copyright (C) 2018
 * Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Alvaro Paz
 * CINVESTAV - Saltillo Campus
 *
 * This file is part of OpenHRC
 * OpenHRC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * OpenHRC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/**
 *	\file
 *	\author Carla Villanueva, Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Class to implement the variational dynamics.
 */

#ifndef HR_VARIATIONAL_DYNAMIC_H
#define HR_VARIATIONAL_DYNAMIC_H

#include <memory>
#include "geombd/core.h"

namespace geo{
namespace variational{


class Dynamic : public Kinematics
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    //! Custom constructor.
    /*! \param id The id of the body.
        \param name The name of the body.*/
    Dynamic(std::shared_ptr< MultiBody > robot, float alpha, float dt): Kinematics(robot){

        this->robot = robot;
        this->n = robot->getDoF();
        this->alpha = alpha;
        this->dt = dt;

        this->D_X = MatrixXr::Identity(n,n);
        this->D_q = D_X;
        this->D_dq = D_X;

        this->D_Xx = MatrixXr::Identity(2*n,2*n);
        this->D_qx = D_Xx.topRows(n);
        this->D_dqx = D_Xx.bottomRows(n);


        DynamicBody* currentBody;
        typename std::vector< Body* >::iterator bodyIt;

        Diff_Twist.conservativeResize(Eigen::NoChange, 2*n);
        Diff_Acceleration.conservativeResize(Eigen::NoChange, 2*n);
        Diff_Wrench.conservativeResize(Eigen::NoChange, 2*n);


        //! Set root Twist and Accelerations
        bodyIt = bodies.begin();
        currentBody = (DynamicBody *)(*bodyIt);
        currentBody->setAcceleration(robot->getGravity());
        currentBody->setDiff_Twist(MatrixXr::Zero(6,2*n));
        currentBody->setDiff_Acceleration(MatrixXr::Zero(6,2*n));


        int body_id;
        int parent_id;
        int kx=0;


        for ( bodyIt = bodies.begin()+1; bodyIt != bodies.end(); bodyIt++ )
        {
            currentBody = (DynamicBody *)(*bodyIt);
            body_id = currentBody->getId();
            parent_id = robot->getParentId(currentBody->getId());
            kx= 0;
            VectorXr ancestors_aux(parent_id*parent_id);

            MatrixXr adjointScrew = ad(currentBody->getScrewAxes());
            currentBody->setadjointScrew(adjointScrew);

            MatrixXr adjointDualScrew = adDual(-currentBody->getScrewAxes());
            currentBody->setadjointDualScrew(adjointDualScrew);

            //! Pre-allocation to spped up computation
            currentBody->setDiff_Twist(Diff_Twist);
            currentBody->setDiff_Acceleration(Diff_Acceleration);
            currentBody->setDiff_Wrench(Diff_Wrench);

            for(int ix = 0; ix<parent_id; ix++){

                for(int jx = ix * body_id; jx<((ix*body_id)+parent_id); jx++){
                    ancestors_aux(kx) = jx;
                    kx++;
                }
            }
            if (body_id>1){
                aux_parent[body_id] = ancestors_aux;
            }
        }

        for ( bodyIt = bodies.begin()+1; bodyIt != bodies.end(); bodyIt++ )
        {
            currentBody = (DynamicBody *)(*bodyIt);
            body_id = currentBody->getId();
            parent_id = robot->getParentId(currentBody->getId());
            kx= 0;
            VectorXr ancestors_aux(4*parent_id*parent_id);
            for(int ix = 0; ix<parent_id; ix++){

                for(int jx = 4 * ix * body_id; jx<((4 * ix*body_id)+ 4 * parent_id); jx++){
                    ancestors_aux(kx) = jx;
                    kx++;
                    //std::cout<<jx<<std::endl;
                }
            }
            if (body_id>1){
                aux_parent_pos[body_id] = ancestors_aux;
            }
        }


    }


    ~Dynamic(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    public:

    //! Robot's Degrees of freedom
    short int n;

    //! Quadrature variable
    float alpha;

    //! Time step
    float dt;



    protected:

    //! MultiBody
    std::shared_ptr< MultiBody > robot;

    //! Lagragian derivative with respect to generalized position
    VectorXr dL_q;

    //! Lagrangian derivative with respecto to generalized velocity
    VectorXr dL_dq;

    //! Second Lagrangian derivative with respect to generalized position
    MatrixXr dL_qq;

    //! Second Lagrangian derivative with respect to generalized velocity and then with respect to generalized position
    MatrixXr dL_qdq;

    //! Second Lagrangian derivative with respect to generalized position and then with respect to generalized velocity
    MatrixXr dL_dqq;

    //! Second Lagrangian derivative with respect to generalized velocity
    MatrixXr dL_dqdq;

    //! This map contains the location of parent tensor information in the current tensor information
    std::map<int,geo::VectorXr> aux_parent;

    //! This map contains the location of parent tensor information in the current tensor information with respect a transformation matrices
    std::map<int,geo::VectorXr> aux_parent_pos;

    Vector3r multibodyCoM;
    MatrixXr D_multibodyCoM;
    MatrixXr D_multibodyCoMx;
    Matrix4r G_CoMi, G_aux;         // Local and global i-th center of mass
    //! Variables for centroidal momentum
    SpatialVector SpatialMomentum, CentroidalMomentum;
    D_SpatialVector D_SpatialMomentum, D_CentroidalMomentum;
    D_SpatialVector Diff_Twist, Diff_Acceleration, Diff_Wrench;

    MatrixXr D_Xx, D_qx, D_dqx, D_ddqx;




    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //! Compute the discrete configuration for the multibody system
         /*! \param A configuration vector in instant of time "k" and "k+1"
         * \return void
         */
    void setDiscreteConfiguration(const VectorXr &qk, const VectorXr &qk1);

    //! Compute the first Lagrangian derivatives with respect to generalized position and velocity
         /*! \param none
         * \return void
         */
    void computeFirstLagrangianDerivative();

    //! Compute the second Lagrangian derivatives with respect generalized position and velocity and their combination
         /*! \param none
         * \return void
         */
    void computeSecondLagrangianDerivative();

    //! Compute the Inertia Matrix of the multibody system
         /*! \param none
         * \return Inertia Matrix
         */
    MatrixXr computeInertiaMatrix();

    //! Compute the vector of Nonlinear terms of the multibody system
         /*! \param none
         * \return Nonlinear Term Vector
         */
    VectorXr computeNonLinearTerms();

    //! Compute the first slot Derivatives that are used to compute the dynamic equations of motion
         /*! \param configuration vectors in instant of time "k", "k+1" and "k-1"
         *   \param Slot Derivative Vector with respect to the first variable
         *   \param Slot Derivative Vector with respect to the second variable
         * \return void
         */
    void computeSlotDerivatives(const VectorXr &q0, const VectorXr &q1, const  VectorXr &q2, VectorXr &D1, VectorXr &D2);

    //! Compute the second slot Derivatives that are used to compute the derivative of  dynamic equations of motion
         /*! \param configuration vectors in instant of time "k", "k+1" and "k-1"
         * \return void
         */
    void computeSecondSlotDerivatives(const VectorXr &q0, const VectorXr &q1, const  VectorXr &q2, MatrixXr &D1D1, MatrixXr &D2D1, MatrixXr &D1D2, MatrixXr &D2D2);


    void computeCenterofMass(const VectorXr &qk, const bool &computeFirstDerivative);


    //! Get the multibody center of mass
    /*! \param none
        \return a Vector3r
             */
    Vector3r getMultibodyCoM(){ return multibodyCoM; }

    //! Get the multibody center of mass differentiation
    /*! \param none
        \return a MatrixXr
             */
    MatrixXr getMultibodyCoMDifferentiation(){ return D_multibodyCoM; }

    //! Compute the first slot Derivatives that are used to compute the dynamic equations of motion
         /*! \param configuration vectors in instant of time "k", "k+1" and "k-1"
         *   \param Slot Derivative Vector with respect to the first variable
         *   \param Slot Derivative Vector with respect to the second variable
         * \return void
         */
    void computeContinuousSlotDerivatives(const VectorXr &q, const VectorXr &dq, VectorXr &D1, VectorXr &D2);

    //! Compute the multibody spatial momentum at inertial frame
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeSpatialMomentum(const VectorXr &qk, const VectorXr &qk1, const bool &computePartialDerivatives);

    //! Compute the multibody centroidal momentum
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCentroidalMomentum(const VectorXr &qk, const VectorXr &qk1, const bool &computePartialDerivatives);

    //! Get the multibody spatial momentum
    /*! \param none
        \return a SpatialVector
             */
    SpatialVector getSpatialMomentum(){ return SpatialMomentum; }

    //! Get the multibody spatial momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getSpatialMomentumDifferentiation(){ return D_SpatialMomentum; }

    //! Get the multibody centroidal momentum
    /*! \param none
        \return a SpatialVector
             */
    SpatialVector getCentroidalMomentum(){ return CentroidalMomentum; }

    //! Get the multibody centroidal momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getCentroidalMomentumDifferentiation(){ return D_CentroidalMomentum; }

    VectorXr getGeneralizedCoordinates(){return q;}

    VectorXr getGeneralizedVelocity(){return dq;}


    protected:


};
} // end of namespace variational
} // end of namespace geo

#endif // HR_VARIATIONAL_DYNAMIC_H
