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

#ifndef HR_DYNAMICS_FORWARD_DYNAMICS_H
#define HR_DYNAMICS_FORWARD_DYNAMICS_H

#include <memory>
#include "geombd/types.h"
#include "geombd/core.h"
#include "geombd/dynamics/Lie_operators.h"

namespace hr{
namespace core{


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
        D_V_a.resize(n);    D_V_b.resize(n);    D_c_a.resize(n);    D_c_b.resize(n);
        D_p_a.resize(n);    D_p_b.resize(n);    D_Pn_a.resize(n);   D_Pn_b.resize(n);
        D_dVa_a.resize(n);  D_dVa_b.resize(n);  D_dV_a.resize(n);   D_dV_b.resize(n);
        D_U_v.resize(n);    D_U_h.resize(n);    D_Ia.resize(n);     D_In.resize(n);     D_invD.resize(n);
        D_u_a.resize(n);    D_u_b.resize(n);    D_ddq_a.resize(n);  D_ddq_b.resize(n);
        this->NP = VectorXr::Zero(n,1);  this->NS = VectorXr::Zero(n,1);

        S_i = Eigen::MatrixXi::Zero(n,2);

        i = 0;
        for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
        {
            tempBody = (DynamicBody *)(*bodyIterator);

            j = robot->getParentId(tempBody->getId());
            parentBody = (DynamicBody *)(bodies.at(j));

            nP = robot->getPredecessorJoints(tempBody->getId()).size(); //! Both consider the i-th joint
            nS = robot->getSubTreeBody(tempBody->getId()).size();

            Pre = robot->getPredecessorJoints(tempBody->getId());
            for(short int pk=0; pk<Pre.size(); pk++){ Pre.at(pk)--; }

            Suc = robot->getSubTreeBody(tempBody->getId());
            Suc.erase (Suc.begin());
            for(short int pk=0; pk<Suc.size(); pk++){ Suc.at(pk)--; }

            PRE.at(i) = Pre;    SUC.at(i) = Suc;    PreSuc = Pre;
            PreSuc.insert(std::end(PreSuc), std::begin(Suc), std::end(Suc));
            PRESUC.at(i) = PreSuc;

            D_V_a.at(i).conservativeResize(Eigen::NoChange, nP);    D_V_b.at(i).conservativeResize(Eigen::NoChange, nP);
            D_c_a.at(i).conservativeResize(Eigen::NoChange, nP);    D_c_b.at(i).conservativeResize(Eigen::NoChange, nP);
            D_p_a.at(i).conservativeResize(Eigen::NoChange, nP);    D_p_b.at(i).conservativeResize(Eigen::NoChange, nP);
            D_Pn_a.at(i).conservativeResize(Eigen::NoChange, nP+nS-1);   D_Pn_b.at(i).conservativeResize(Eigen::NoChange, nP+nS-1);
            D_In.at(i).resize(6*(nS-1),6);     D_Ia.at(i).resize(6*(nS-1),6);    D_invD.at(i).resize(1,nS-1);
            D_U_v.at(i).resize(6*(nS-1),1);    D_U_h.at(i).resize(6*(nS-1),1);
            D_u_a.at(i).resize(1,n);   D_u_b.at(i).resize(1,n);
            D_ddq_a.at(i).resize(1,n);   D_ddq_b.at(i).resize(1,n);
            D_dVa_a.at(i).conservativeResize(Eigen::NoChange, n);   D_dVa_b.at(i).conservativeResize(Eigen::NoChange, n);
            D_dV_a.at(i).conservativeResize(Eigen::NoChange, n);    D_dV_b.at(i).conservativeResize(Eigen::NoChange, n);

            NP(i) = nP;   NS(i) = nS;

            real_t tempScrew = tempBody->getScrewAxes().transpose()*screwIndexes;
            if(tempBody->getScrewAxes().sum() == 1){S_i.row(i)<<1, (int) tempScrew;}

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
    SpatialMatrix Adjoint, adjoint, AdjointDual, adjointDual;
    SpatialVector quasiAcceleration;

    real_t D, invD, u;
    int idx_01, idx_02, idx_03;
    SpatialVector Screwi, U, Pa, Twist, C_bias;
    SpatialMatrix Ia, In;

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
    MatrixXrColMajor Aux_III;
    D_SpatialVector Diff_Twist, Diff_Acceleration, Diff_Wrench, Diff_Pa, Diff_Pn;
    SpatialMatrix adjointScrew, adjZ;
    SpatialVector screwIndexes;
    short int in_screw;
    MatrixXr Diff_Ufd, Diff_Ia, IdentityNxN;
    RowVectorXr Diff_invD, Diff_ufd;
    std::vector< MatrixXr > Diff_In;

    std::vector<short int> s_tree;

    //! ------------------------------------------
    //! Set variables for Enhanced differentiation
    typedef std::vector<D_SpatialVector> D_Vector;
    typedef std::vector<MatrixXrColMajor> D_MatrixColMajor;
    typedef std::vector<MatrixXr> D_Matrix;

    D_Vector D_V_a, D_V_b, D_c_a, D_c_b, D_p_a, D_p_b, D_Pn_a, D_Pn_b, D_dVa_a, D_dVa_b, D_dV_a, D_dV_b;
    D_MatrixColMajor D_U_v, D_U_h;
    D_Matrix D_Ia, D_In, D_invD, D_u_a, D_u_b, D_ddq_a, D_ddq_b;
    D_SpatialVector D_Vb;

    short int nP, nS, PS, nPj;
    VectorXr NP, NS;
    SpatialMatrix D_p_aux;
    MatrixXr D_Pa_a, D_Pa_b, D_aux;

    std::vector<short int> Pre, Suc, PreSuc;
    std::vector<std::vector<short int>> PRE, SUC, PRESUC;

    Eigen::MatrixXi S_i;

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


    protected:


};
} // end of namespace core
} // end of namespace hr

#endif // HR_DYNAMICS_FORWARD_DYNAMICS_H
