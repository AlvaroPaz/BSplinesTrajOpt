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

#ifndef GEOMBD_DYNAMICS_INVERSE_DYNAMICS_H
#define GEOMBD_DYNAMICS_INVERSE_DYNAMICS_H

#include <memory>
#include "geombd/core.h"
#include "geombd/VanishProduct"

//int Vanisher_i = 0, Vanisher_n = 1, Vanisher_m = 2;
//int Vanisher_i, Vanisher_n, Vanisher_m;

namespace geo{

class InverseDynamics : public Kinematics
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor.
    InverseDynamics(std::shared_ptr< MultiBody > robot) : Kinematics(robot) {

        this->robot = robot;
        this->n = robot->getDoF();
        this->m = robot->getDifferentiationSize();

        ::Vanisher_n = this->n;   ::Vanisher_m = this->m;

        this->Diff_Tau = MatrixXr::Zero(n,m);
        this->DDiff_Tau = MatrixXr::Zero(n,m*m);
        this->DDiff_Tau_Contracted = MatrixXr::Zero(n,m*(m+1)/2);

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
        this->DD_q_Van = MatrixXr::Zero(n,m*(m+1)/2);


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
        int ID = 0;  hasChild.resize(n);
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
            DD_q_Van.row(ID) = Eigen::vanishProduct(D_q.row(ID),D_q.row(ID));

            //! Pre-allocation to spped up computation
            tempBody->setDiff_Twist(Diff_Twist);
            tempBody->setDiff_Acceleration(Diff_Acceleration);
            tempBody->setDiff_Wrench(Diff_Wrench);

            DD_Twist.at(ID).conservativeResize(Eigen::NoChange, m*(m+1)/2);
            DD_Acceleration.at(ID).conservativeResize(Eigen::NoChange, m*(m+1)/2);
            DD_Wrench.at(ID).conservativeResize(Eigen::NoChange, m*(m+1)/2);

            hasChild(ID) = robot->getBodyChildren(tempBody->getId()).size();

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

        //! Build Vanished Kronecker Commutation Matrix
        KCM_van.resize(m*m, m*(m+1)/2);
        int j_ = 0, k_ = 0, w_ = 0;
        for(i_ = 0 ; i_ < m*m ; i_ += m){
            KCM_van.middleCols(w_, m-k_) = KCM.middleCols(j_, m-k_);
            w_ += m - k_;
            k_++;
            j_ += m + 1;
        }

        //! Variables for Center of Mass
        this->multibodyMass = robot->getTotalMass();
        this->G_CoMi(3) = 1;

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
    MatrixXr DD_q, DD_q_Van;

    SpatialVector Twist, Acceleration, Wrench;
    MatrixXr Aux_I, Aux_II, Diff_Tau;
    D_SpatialVector Diff_Twist, Diff_Acceleration, Diff_Wrench;
    SpatialMatrix adjointScrew, adjointDualScrew;
    D_SpatialVector inputI, inputII;
    D_SpatialVector inputIII,inputIV,inputV;
    VectorXr hasChild;

    std::vector< Matrix3r > SK, SK_2;

    //! Variables for second differentiation
    D_SpatialVector DDiff_Twist, DDiff_Acceleration, DDiff_Wrench, ChildDDiff_Wrench;
    std::vector< D_SpatialVector > DD_Twist, DD_Acceleration, DD_Wrench;
    MatrixXr DDiff_Tau, DDiff_Tau_Contracted;
    SparseMatrixXr KCM;
    SparseMatrixXr KCM_van;

    //! Variables for Center of Mass
    Vector3r multibodyCoM;          // The multibody center of mass
    real_t multibodyMass;           // The multibody total mass
    Vector4r G_CoMi, G_aux;         // Local and global i-th center of mass
    Matrix4r D_G_x;
    short int tempBodyId;           // Temporal body ID
    std::vector< MatrixXr > D_G, D_Gv;    // Differentiation of G
    std::vector< MatrixXr > DD_G;   // Second differentiation of G
    MatrixXr D_multibodyCoM;        // Differentiation of the multibody center of mass
    MatrixXr DD_multibodyCoM;       // Second differentiation of the multibody center of mass
    MatrixXr DD_ContractedCoM;       // Second differentiation of the multibody center of mass
    short int sizeDecisionVector, parentBodyId;

    //! Flags for updating
    bool FK_updated, CoM_updated, KCM_built;

    //! Variables for centroidal momentum
    SpatialVector SpatialMomentum, CentroidalMomentum;
    D_SpatialVector D_SpatialMomentum, D_CentroidalMomentum;
    D_SpatialVector DD_SpatialMomentum, DD_ContractedSpatialMomentum, DD_CentroidalMomentum, DD_ContractedMu;

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

    //! Compute the Inverse Dynamics for the multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeContractedInverseDynamics(const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute the multibody spatial momentum at inertial frame
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeContractedSpatialMomentum(const bool &firstDerivative, const bool &secondDerivative);

    //! Compute the multibody spatial momentum at inertial frame
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeSpatialMomentum(const bool &firstDerivative, const bool &secondDerivative);

    //! Compute the multibody centroidal momentum
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCentroidalMomentumII(const bool &firstDerivative, const bool &secondDerivative);

    //! Compute the multibody centroidal momentum
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCentroidalMomentum(const bool &firstDerivative, const bool &secondDerivative);

    //! Compute the multibody centroidal momentum
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeContractedCentroidalMomentum(const bool &firstDerivative, const bool &secondDerivative);

    //! Compute the center of mass of the whole multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCenterOfMass(const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute the center of mass of the whole multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeCenterOfMassII(const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute the center of mass of the whole multibody system
         /*! \param boolean flag for computing partials
         * \return void
         */
    void computeContractedCenterOfMass(const bool &computeFirstDerivative, const bool &computeSecondDerivative);

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

    //! Set generalized coordinates differentiation
         /*! \param D configuration, D generalized velocity, D generalized acceleration
         * \return void
         */
    void setContractedDD_q( const MatrixXr &DD_q_Van );

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

    //! Get generalized torques second differentiation
    /*! \param none
        \return a real_t eigen matrix
             */
    MatrixXr getContractedDDTorque(){ return DDiff_Tau_Contracted; }

    //! Get the multibody center of mass
    /*! \param none
        \return a Vector3r
             */
    Vector3r getRobotCoM(){ return multibodyCoM; }

    //! Get the multibody center of mass differentiation
    /*! \param none
        \return a MatrixXr
             */
    MatrixXr getRobotD_CoM(){ return D_multibodyCoM; }

    //! Get the multibody center of mass differentiation
    /*! \param none
        \return a MatrixXr
             */
    MatrixXr getRobotDD_CoM(){ return DD_multibodyCoM; }

    //! Get the multibody center of mass differentiation
    /*! \param none
        \return a MatrixXr
             */
    MatrixXr getRobotDD_ContractedCoM(){ return DD_ContractedCoM; }

    //! Get the multibody spatial momentum
    /*! \param none
        \return a SpatialVector
             */
    SpatialVector getSpatialMomentum(){ return SpatialMomentum; }

    //! Get the multibody spatial momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getD_SpatialMomentum(){ return D_SpatialMomentum; }

    //! Get the multibody spatial momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getDD_SpatialMomentum(){ return DD_SpatialMomentum; }

    //! Get the multibody spatial momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getDD_ContractedSpatialMomentum(){ return DD_ContractedSpatialMomentum; }

    //! Get the multibody centroidal momentum
    /*! \param none
        \return a SpatialVector
             */
    SpatialVector getCentroidalMomentum(){ return CentroidalMomentum; }

    //! Get the multibody centroidal momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getD_CentroidalMomentum(){ return D_CentroidalMomentum; }

    //! Get the multibody centroidal momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getDD_CentroidalMomentum(){ return DD_CentroidalMomentum; }

    //! Get the multibody centroidal momentum differentiation
    /*! \param none
        \return a D_SpatialVector
             */
    D_SpatialVector getDD_ContractedCentroidalMomentum(){ return DD_ContractedMu; }


protected:


};
} // end of namespace geo

#endif // HR_DYNAMICS_INVERSE_DYNAMICS_H
