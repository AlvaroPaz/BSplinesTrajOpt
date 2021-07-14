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
 *	\file include/geombd/dynamics/Lie_operators.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Class to implement the inverse dynamics.
 */

#ifndef HR_DYNAMICS_INVERSE_DYNAMICS_H
#define HR_DYNAMICS_INVERSE_DYNAMICS_H

#include <memory>
#include "geombd/types.h"
#include "geombd/core.h"
#include "geombd/dynamics/Lie_operators.h"
#include "geombd/dynamics/kinematics.h"

namespace hr{
namespace core{


class InverseDynamics : public LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor.
    InverseDynamics(std::shared_ptr< MultiBody > robot) : robotKinematics{ robot } {

        this->robot = robot;
        this->n = robot->getDoF();
        this->m = robot->getDifferentiationSize();

        this->Diff_Tau = MatrixXr::Zero(n,m);
        this->DDiff_Tau = MatrixXr::Zero(n,m*m);

        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->ddq = VectorXr::Zero(n,1);
        this->Tau = VectorXr::Zero(n,1);

        this->Identity3x3.setIdentity();
        this->e_sq.setIdentity();
        this->bodies = robot->getPubBodies();

        this->FK_updated = false;
        this->CoM_updated = false;
        this->KCM_built = false;

        //! Variables for differentiation //!<$ Default <- wrt_state

        this->D_X = MatrixXr::Identity(m,m);
        this->D_q = D_X.topRows(n);
        this->D_dq = D_X.bottomRows(n);
        this->D_ddq = MatrixXr::Zero(n,m);


        this->DD_q = MatrixXr::Zero(n,m*m);


        Diff_Twist.conservativeResize(Eigen::NoChange, m);
        Diff_Acceleration.conservativeResize(Eigen::NoChange, m);
        Diff_Wrench.conservativeResize(Eigen::NoChange, m);

        inputI.conservativeResize(Eigen::NoChange, m);
        inputII.conservativeResize(Eigen::NoChange, m);
        inputIII.conservativeResize(Eigen::NoChange, m*m);
        inputIV.conservativeResize(Eigen::NoChange, m*m);
        inputV.conservativeResize(Eigen::NoChange, m*m);

        DDiff_Twist.conservativeResize(Eigen::NoChange, m*m);
        DDiff_Acceleration.conservativeResize(Eigen::NoChange, m*m);
        DDiff_Wrench.conservativeResize(Eigen::NoChange, m*m);
        ChildDDiff_Wrench.conservativeResize(Eigen::NoChange, m*m);


        //! Set root Twist and Accelerations
        bodyIterator = bodies.begin();
        tempBody = (DynamicBody *)(*bodyIterator);
        tempBody->setAcceleration(robot->getGravity());
        tempBody->setDiff_Twist(MatrixXr::Zero(6,m));
        tempBody->setDiff_Acceleration(MatrixXr::Zero(6,m));


        DD_Twist.resize(n);      DD_Acceleration.resize(n);      DD_Wrench.resize(n);

        SK.clear();             SK_2.clear();
        //! Forward recursion
        int ID = 0;
        for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
        {
            tempBody = (DynamicBody *)(*bodyIterator);

            sk = skew(tempBody->getScrewAxes().segment(3,3));

            SK.push_back(sk);
            SK_2.push_back(sk*sk);

            adjointScrew = ad(tempBody->getScrewAxes());
            tempBody->setadjointScrew(adjointScrew);

            adjointDualScrew = adDual(-tempBody->getScrewAxes());
            tempBody->setadjointDualScrew(adjointDualScrew);

            //! Compute the constant term kron(D_q,D_q)
            DD_q.row(ID) = Eigen::kroneckerProduct(D_q.row(ID),D_q.row(ID));

            //! Pre-allocation to spped up computation
            tempBody->setDiff_Twist(Diff_Twist);
            tempBody->setDiff_Acceleration(Diff_Acceleration);
            tempBody->setDiff_Wrench(Diff_Wrench);

            DD_Twist.at(ID).conservativeResize(Eigen::NoChange, m*m);
            DD_Acceleration.at(ID).conservativeResize(Eigen::NoChange, m*m);
            DD_Wrench.at(ID).conservativeResize(Eigen::NoChange, m*m);

            ID++;
        }


        //! Build Kronecker Commutation Matrix
        TripletXrList TripletKCM;
        int i_ = 0;
        for(int k_ = 0 ; k_ < m ; k_++){
            for(int j_ = k_ ; j_ < m*m ; j_ += m){
                TripletKCM.push_back( TripletXr(i_,j_,1) );
                i_++;
            }
        }
        SparseMatrixXr KCMaux(m*m,m*m);
        KCMaux.setFromTriplets(TripletKCM.begin(), TripletKCM.end());
        this->KCM = KCMaux;
        KCMaux.setIdentity();
        this->KCM += KCMaux;


        //! Variables for Center of Mass
        this->multibodyMass = robot->getTotalMass();
        this->G_CoMi.setIdentity();        

    }


    ~InverseDynamics(){}


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

    //! Robot kinematics object
    Kinematics robotKinematics;


    ///! Initialize variables
    short int j, i;

    VectorXr q, dq, ddq, Tau;
    Matrix3r sk, sk_2, e_wq, T_wq, Identity3x3;

    Matrix4r e_sq;

    SpatialMatrix Adjoint, adjoint, AdjointDual, adjointDual, SpatialInertia, SpatialInertiaTmp;
    SpatialVector Screwi;

    typename std::vector< Body* >::iterator bodyIterator;
    DynamicBody* tempBody;
    DynamicBody* parentBody;
    std::vector<Body*> bodies;
    std::vector<short int> bodyChildren;
    DynamicBody* children;
    typename std::vector< short int >::iterator childrenIterator;

    //! Variables for first differentiation
    MatrixXr D_X, D_q, D_dq, D_ddq;
    MatrixXr DD_q;

    SpatialVector Twist, Acceleration, Wrench;
    MatrixXr Aux_I, Aux_II, Diff_Tau;
    D_SpatialVector Diff_Twist, Diff_Acceleration, Diff_Wrench;
    SpatialMatrix adjointScrew, adjointDualScrew;
    D_SpatialVector inputI, inputII;
    D_SpatialVector inputIII,inputIV,inputV;

    std::vector< Matrix3r > SK, SK_2;

    //! Variables for second differentiation
    D_SpatialVector DDiff_Twist, DDiff_Acceleration, DDiff_Wrench, ChildDDiff_Wrench;
    std::vector< D_SpatialVector > DD_Twist, DD_Acceleration, DD_Wrench;
    MatrixXr DDiff_Tau;
    SparseMatrixXr KCM;

    //! Variables for Center of Mass
    Vector3r multibodyCoM;          // The multibody center of mass
    real_t multibodyMass;           // The multibody total mass
    Matrix4r G_CoMi, G_aux;         // Local and global i-th center of mass
    short int tempBodyId;           // Temporal body ID
    std::vector< MatrixXr > D_G;    // Differentiation of G
    MatrixXr D_multibodyCoM;        // Differentiation of the multibody center of mass
    short int sizeDecisionVector, parentBodyId;

    //! Flags for updating
    bool FK_updated, CoM_updated, KCM_built;

    //! Variables for centroidal momentum
    SpatialVector SpatialMomentum, CentroidalMomentum;
    D_SpatialVector D_SpatialMomentum, D_CentroidalMomentum;

    public:



    // --------------------------------------------
    // Methods
    // --------------------------------------------

    private:

    public:

    //! Compute my own sparse Kronecker product
         /*! \param ID for sparsity pattern, right and left operands
         * \return void
         */
    void myOwnKroneckerProduct(const short int ID__, const D_real_t &leftElement, const D_SpatialVector &rightElement);

    //! Compute the Inverse Dynamics for the multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeInverseDynamics(const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute the multibody spatial momentum at inertial frame
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeSpatialMomentum(const bool &computePartialDerivatives);

    //! Compute the multibody centroidal momentum
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCentroidalMomentum(const bool &computePartialDerivatives);

    //! Compute the center of mass of the whole multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCenterOfMass(const bool &computePartialDerivatives);

    //! Set generalized coordinates
         /*! \param configuration, generalized velocity, generalized acceleration
         * \return void
         */
    void setGeneralizedCoordinates( const VectorXr &q, const VectorXr &dq, const VectorXr &ddq );

    //! Set generalized coordinates differentiation
         /*! \param D configuration, D generalized velocity, D generalized acceleration
         * \return void
         */
    void setGeneralizedCoordinatesDifferentiation( const MatrixXr &D_q, const MatrixXr &D_dq, const MatrixXr &D_ddq );

    //! Set generalized coordinates differentiation
         /*! \param D configuration, D generalized velocity, D generalized acceleration
         * \return void
         */
    void setGeneralizedCoordinatesSecondDifferentiation( const MatrixXr &DD_q );

    //! Get generalized torques
    /*! \param none
        \return a real_t eigen vector
             */
    VectorXr getGeneralizedTorques(){ return Tau; }

    //! Get generalized torques first differentiation
    /*! \param none
        \return a real_t eigen matrix
             */
    MatrixXr getGeneralizedTorquesFirstDifferentiation(){ return Diff_Tau; }

    //! Get generalized torques second differentiation
    /*! \param none
        \return a real_t eigen matrix
             */
    MatrixXr getGeneralizedTorquesSecondDifferentiation(){ return DDiff_Tau; }

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


protected:


};
} // end of namespace core
} // end of namespace hr

#endif // HR_DYNAMICS_INVERSE_DYNAMICS_H
