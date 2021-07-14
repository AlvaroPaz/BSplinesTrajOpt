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
 *	Definition of Lie operators.
 */

#ifndef HR_DYNAMICS_LIE_OPERATORS_H
#define HR_DYNAMICS_LIE_OPERATORS_H

#include "geombd/types.h"
#include "geombd/core.h"

#define PI 3.141592653589793

using std::cout;
using std::cin;
using std::endl;
using std::string;

namespace hr{
namespace core{

class LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------


    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Default constructor.
    /*! \param .
        \param */
    LieOperators(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------


    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:
        /*! Description: Computes the 3x3 skew symmetric matrix from a 3D vector
             * \parameters: 3D vector
             * \return: Skew symmetric matrix
             */
        Matrix3r skew(const Vector3r &vec);

        /*! Description: Upgrades the screw vector into its matrix form
             * \parameters: Spatial screw vector
             * \return: Screw matrix representation
             */
        Matrix4r upScrew(const SpatialVector &screw);

        /*! Description: Computes the Adjoint operator Ad: se(3) -> se(3)
             * \parameters: An element of the group SE(3)
             * \return: Adjoint operator
             */
        SpatialMatrix Ad(const Matrix4r &G);

        /*! Description: Computes the dual Adjoint operator Ad*: se*(3) -> se*(3)
             * \parameters: An element of the group SE(3)
             * \return: Dual Adjoint operator
             */
        SpatialMatrix AdDual(const Matrix4r &G);

        /*! Description: Computes the inverse dual Adjoint operator Ad*: se*(3) -> se*(3)
             * \parameters: An element of the group SE(3)
             * \return: Inverse dual Adjoint operator
             */
        SpatialMatrix AdDualInv(const Matrix4r &G);

        /*! Description: Computes the adjoint operator ad: se(3) -> se(3)
             * \parameters: An element of the algebra se(3)
             * \return: adjoint operator
             */
        SpatialMatrix ad(const SpatialVector &V);

        /*! Description: Computes the dual adjoint operator ad*: se*(3) -> se*(3)
             * \parameters: An element of the algebra se(3)
             * \return: adjoint operator
             */
        SpatialMatrix adDual(const SpatialVector &V);

        /*! Description: Computes the bar adjoint operator ad^{_}
             * \parameters: An element of the algebra se*(3)
             * \return: adjoint operator
             */
        SpatialMatrix adBar(const SpatialVector &F);

};
} // end of namespace core
} // end of namespace hr

#endif // HR_DYNAMICS_LIE_OPERATORS_H
