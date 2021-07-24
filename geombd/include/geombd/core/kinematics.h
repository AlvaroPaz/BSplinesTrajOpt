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
 *	\file include/geombd/dynamics/kinematics.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Class to implement the multibody kinematics.
 */

#ifndef GEOMBD_CORE_KINEMATICS_H
#define GEOMBD_CORE_KINEMATICS_H

#include <memory>
#include "geombd/core.h"

namespace geo{

class Kinematics : public LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Default constructor.
    /*! \param .
        \param */
    Kinematics( ){}

    //! Custom constructor.
    Kinematics(std::shared_ptr< MultiBody > robot){

        this->robot = robot;
        this->n = robot->getDoF();

        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->ddq = VectorXr::Zero(n,1);

        this->Identity3x3.setIdentity();
        this->e_sq.setIdentity();
        this->bodies = robot->getPubBodies();

        //! Variables for differentiation //!<$ Default <- wrt_state

        this->D_X = MatrixXr::Identity(2*n,2*n);
        this->D_q = D_X.block(0,0,n,2*n);
        this->D_dq = D_X.block(n,0,n,2*n);
        this->D_ddq = MatrixXr::Zero(n,2*n);

        SK.clear();    SK_2.clear();
        //! Forward recursion
        for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
        {
            tempBody = (DynamicBody *)(*bodyIterator);

            sk = skew(tempBody->getScrewAxes().segment(3,3));

            SK.push_back(sk);
            SK_2.push_back(sk*sk);

        }

    }


    ~Kinematics(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    protected:

    //! Dof
    short int n;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;


    ///! Initialize variables
    VectorXr q, dq, ddq;
    Matrix3r sk, sk_2, e_wq, T_wq, Identity3x3;

    Matrix4r e_sq;

    typename std::vector< Body* >::iterator bodyIterator;
    DynamicBody* tempBody;
    DynamicBody* parentBody;
    std::vector<Body*> bodies;
    std::vector<short int> bodyChildren;
    DynamicBody* children;
    typename std::vector< short int >::iterator childrenIterator;

    //! Variables for differentiation

    MatrixXr D_X, D_q, D_dq, D_ddq;

    std::vector< Eigen::Matrix<real_t, 3, 3> > SK, SK_2;

    public:



    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //! Compute the whole multibody forward kinematics
         /*! \param null
         * \return void
         */
    void computeForwardKinematics( );    

    //! Compute the predeccessors forward kinematics
         /*! \param DynamicBody*
         * \return void
         */
    void computeForwardKinematics( DynamicBody* currentBody );

    //! Compute the whole multibody forward kinematics and twist simultaneously
         /*! \param null
         * \return void
         */
    void computeKinematics( );

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


protected:


};
} // end of namespace geo

#endif // GEOMBD_CORE_KINEMATICS_H
