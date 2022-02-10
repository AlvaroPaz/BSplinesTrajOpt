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
 *	Class to implement the forward dynamics.
 */

#ifndef GEOMBD_DYNAMICS_FORWARD_DYNAMICS_H
#define GEOMBD_DYNAMICS_FORWARD_DYNAMICS_H

#include <memory>
#include "geombd/core.h"

namespace geo{

class ForwardDynamics : public LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor.
    /*! \param id The id of the body.
        \param name The name of the body.*/
    ForwardDynamics(std::shared_ptr< MultiBody > robot){

        this->robot = robot;
        this->n = robot->getDoF();
        this->ddq = VectorXr::Zero(n,1);
        this->D_ddq = MatrixXr::Zero(n,2*n);
        this->inv_H = MatrixXr::Zero(n,n);

        //!---------------------------------

        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->Tau = VectorXr::Zero(n,1);

        this->Identity3x3.setIdentity();
        this->Zeros6x1.setZero();
        this->e_sq.setIdentity();
        this->bodies = robot->getPubBodies();

        //! Variables for differentiation
        this->D_X = MatrixXr::Identity(2*n,2*n);
        this->D_q = D_X.block(0,0,n,2*n);
        this->D_dq = D_X.block(n,0,n,2*n); // Change when considering B-Splines

        this->Aux_I = MatrixXr::Zero(6*2*n,6);
        this->Aux_II = MatrixXr::Zero(6,2*n);
        this->Diff_Tau = MatrixXr::Zero(n,2*n);

        this->adjZ.setZero();    this->adjZ.block(0,3,3,3).setIdentity();    this->adjZ.block(3,0,3,3).setIdentity();
        this->screwIndexes << 0,1,2,3,4,5;

        this->Diff_Ufd = MatrixXr::Zero(6,2*n);
        this->Diff_Ia = MatrixXr::Zero(6,6*n);
        this->Diff_Tau = MatrixXr::Zero(n,2*n);

        this->Diff_Wrench = MatrixXr::Zero(6,2*n);
        this->Diff_Acceleration = MatrixXr::Zero(6,2*n);

        this->IdentityNxN = MatrixXr::Identity(n,n);

        this->Diff_invD = RowVectorXr::Zero(1,2*n);
        this->Diff_ufd = RowVectorXr::Zero(1,2*n);

        this->Diff_In.assign( n, MatrixXr::Zero(6,6*n) );

        this->s_tree = robot->getSubTree();

        //! Set root Twist and Accelerations
        bodyIterator = bodies.begin();
        tempBody = (DynamicBody *)(*bodyIterator);
        tempBody->setAcceleration(robot->getGravity());
        tempBody->setDiff_Twist(MatrixXr::Zero(6,2*n));
        tempBody->setDiff_Acceleration(MatrixXr::Zero(6,2*n));

        //! Pre-allocation for enhancement
        PRE.resize(n);      SUC.resize(n);      PRESUC.resize(n);
        this->NP = VectorXr::Zero(n,1);  this->NS = VectorXr::Zero(n,1);

        S_i = Eigen::MatrixXi::Zero(n,2);

        SK.clear();             SK_2.clear();
        i = 0;
        for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
        {
            tempBody = (DynamicBody *)(*bodyIterator);

            j = robot->getParentId(tempBody->getId());
            parentBody = (DynamicBody *)(bodies.at(j));

            nP = robot->getPredecessorJoints(tempBody->getId()).size(); //! Both consider the i-th joint
            nS = robot->getSubTreeBody(tempBody->getId()).size();

            Pre = robot->getPredecessorJoints(tempBody->getId());
            for(int pk=0; pk<Pre.size(); pk++){ Pre.at(pk)--; }

            Suc = robot->getSubTreeBody(tempBody->getId());
            Suc.erase (Suc.begin());
            for(int pk=0; pk<Suc.size(); pk++){ Suc.at(pk)--; }

            PRE.at(i) = Pre;    SUC.at(i) = Suc;    PreSuc = Pre;
            PreSuc.insert(std::end(PreSuc), std::begin(Suc), std::end(Suc));
            PRESUC.at(i) = PreSuc;


            tempBody->D_q_V.conservativeResize(Eigen::NoChange, nP);
            tempBody->D_dq_V.conservativeResize(Eigen::NoChange, nP);
            tempBody->D_q_c.conservativeResize(Eigen::NoChange, nP);
            tempBody->D_dq_c.conservativeResize(Eigen::NoChange, nP);
            tempBody->D_q_p.conservativeResize(Eigen::NoChange, nP);
            tempBody->D_dq_p.conservativeResize(Eigen::NoChange, nP);

            tempBody->D_U_h.conservativeResize(Eigen::NoChange, nS-1);
            tempBody->D_U_v.conservativeResize(6*(nS-1), Eigen::NoChange);
            tempBody->D_invD.conservativeResize(Eigen::NoChange, nS-1);
            tempBody->D_q_u.conservativeResize(Eigen::NoChange, nP+nS-1);
            tempBody->D_dq_u.conservativeResize(Eigen::NoChange, nP+nS-1);

            tempBody->D_Ia.conservativeResize(6*(nS-1), Eigen::NoChange);
            tempBody->D_q_Pa.conservativeResize(Eigen::NoChange, nP+nS-1);
            tempBody->D_dq_Pa.conservativeResize(Eigen::NoChange, nP+nS-1);

            tempBody->D_In.conservativeResize(6*(nS-1), Eigen::NoChange);
            tempBody->D_q_Pn.conservativeResize(Eigen::NoChange, nP+nS-1);
            tempBody->D_dq_Pn.conservativeResize(Eigen::NoChange, nP+nS-1);

            tempBody->D_q_dVa.conservativeResize(Eigen::NoChange, n);
            tempBody->D_dq_dVa.conservativeResize(Eigen::NoChange, n);
            tempBody->D_q_ddq.conservativeResize(Eigen::NoChange, n);
            tempBody->D_dq_ddq.conservativeResize(Eigen::NoChange, n);
            tempBody->D_q_dV.conservativeResize(Eigen::NoChange, n);
            tempBody->D_dq_dV.conservativeResize(Eigen::NoChange, n);

            tempBody->D_IaC_h.conservativeResize(Eigen::NoChange, nS-1);
            tempBody->D_IaC_v.conservativeResize(6*(nS-1), Eigen::NoChange);

            if(i==0){
                tempBody->D_q_dVa = Zeros6x1;
                tempBody->D_dq_dVa = Zeros6x1;
            }


            NP(i) = nP;   NS(i) = nS;

            Screwi = tempBody->getScrewAxes();

            real_t tempScrew = Screwi.transpose()*screwIndexes;
            if(Screwi.sum() == 1){S_i.row(i)<<1, (int) tempScrew;}


            sk = skew(Screwi.segment(3,3));

            SK.push_back(sk);
            SK_2.push_back(sk*sk);

            adjointScrew = ad(Screwi);
            tempBody->setadjointScrew(adjointScrew);

            adjointDualScrew = adDual(-Screwi);
            tempBody->setadjointDualScrew(adjointDualScrew);

            i++;
        }

    }


    ~ForwardDynamics(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    protected:

    //! Dof
    short int n;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;


    public:

    //! Joint acceleration vector
    MatrixXr ddq;

    //! First derivative of joint acceleration vector
    MatrixXr D_ddq;

    //! Inverse inertia matrix
    MatrixXr inv_H;

    ////////////////////////////////////////////////////////

    protected:

    short int j, i;

    VectorXr q, dq, Tau;
    Matrix3r sk, sk_2, e_wq, T_wq, Identity3x3;
    Matrix4r e_sq;
    SpatialMatrix Adjoint, adjoint, AdjointDual, adjointDual, AdjointDualadDual;
    SpatialVector quasiAcceleration, Zeros6x1;

    real_t D, invD, u;
    int idx_01, idx_02, idx_03;
    SpatialVector Screwi, ScrewiInvD, U, UinvD, Pa, Twist, C_bias;
    SpatialMatrix Ia, In, UtransU;

    typename std::vector< Body* >::iterator bodyIterator;
    DynamicBody* tempBody;
    DynamicBody* parentBody;
    std::vector<Body*> bodies;
    std::vector<short int> bodyChildren;
    DynamicBody* children;
    typename std::vector< short int >::iterator childrenIterator;

    //! Variables for differentiation
    MatrixXr D_X, D_q, D_dq;
    MatrixXr Aux_I, Aux_II, Diff_Tau;
    D_SpatialVector Diff_Twist, Diff_Acceleration, Diff_Wrench, Diff_Pa, Diff_Pn, Diff_VectorAux;
    SpatialMatrix adjointScrew, adjointDualScrew, adjZ, SpatialMatrixAux;
    SpatialVector screwIndexes;
    short int in_screw;
    MatrixXr Diff_Ufd, Diff_Ia, IdentityNxN;
    RowVectorXr Diff_invD, Diff_ufd;
    std::vector< MatrixXr > Diff_In;

    std::vector<short int> s_tree;

    //! ------------------------------------------
    //! Set variables for Enhanced differentiation

    std::vector< Matrix3r > SK, SK_2;

    D_SpatialVector D_Vb;

    short int nP, nS, PS, nPj, nSj;
    VectorXr NP, NS;
    SpatialMatrix D_p_aux;

    std::vector<short int> Pre, Suc, PreSuc;
    std::vector<std::vector<short int>> PRE, SUC, PRESUC;

    Eigen::MatrixXi S_i;

    //! Variables for spatial inertias transformation
    RowMatrix3r Rpro;
    Vector3r Ppro;
    Vector9r Rv;
    Matrix9r Routter;
    SpatialMatrix Rstack;
    Matrix3r Ablock, Bblock, Cblock;
    SpatialVector ABCblocks, Avec, Bvec, Cvec;
    real_t p0, p1, p2;

    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //! Compute the Forward Dynamics for the multibody system
         /*! \param none
         * \return void
         */
    void computeForwardDynamics(const bool &computePartialDerivatives);

    //! Compute the Enhanced Forward Dynamics for the multibody system
         /*! \param none
         * \return void
         */
    void computeEnhancedForwardDynamics(const bool &computePartialDerivatives);

    //! Compute the Inverse of the Inertia Matrix
         /*! \param none
         * \return void
         */
    void computeInverseInertiaMatrix();

    //! Get the joint acceleration vector
    /*! \param none
        \return a Xr vector
             */
    VectorXr getJointAcceleration(){ return ddq; }

    //! Get the first derivative of joint acceleration vector
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getDiffJointAcceleration(){ return D_ddq; }

    //! Get the inverse of inertia matrix
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getInverseInertiaMatrix(){ return inv_H; }

    //! Spatial inertia matrix transformation
         /*! \param Spatial inertia and a SE(3) element
         * \return void
         */
    void transformSpatialInertiaMatrix(SpatialMatrix& SpatialInertia, Matrix4r& G);

    //! Spatial inertia matrix transformation
         /*! \param Spatial inertia
         * \return void
         */
    void transformSpatialInertiaMatrix(SpatialMatrix& SpatialInertia);


    protected:


};
} // end of namespace geo

#endif // GEOMBD_DYNAMICS_FORWARD_DYNAMICS_H
