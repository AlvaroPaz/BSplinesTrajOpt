/*
 * MIT License
 *
 * Copyright (c) 2020 Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Alvaro Paz <alvaro.paz@cinvestav.edu.mx>
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
 *	\file include/geombd/types.h
 *	\author Gustavo Arechavaleta, Alvaro Paz
 *	\version 1.0
 *	\date 2020
 *
 *	Definition of Types
 */

#ifndef GEOMBD_TYPES_H
#define GEOMBD_TYPES_H

//#define EIGEN_NO_MALLOC
//#define EIGEN_RUNTIME_NO_MALLOC
//#define EIGEN_NO_DEBUG

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/../unsupported/Eigen/KroneckerProduct"
//#include "geombd/VanishProduct"
//#include "Eigen/src/Core/Map.h"
#include <chrono>


//! Namespace geombd. Main namespace of GeoMBD
/*!
    More information about geombd.
*/
namespace geo{

//! Defines TYPE real_t for facilitating the switching between double and float.
//#ifdef  __USE_SINGLE_PRECISION__
//        using real_t = float;  //c++11
//! Defines real_t TYPE as float.
//typedef float real_t;
//#else
//! Defines real_t TYPE as double.
typedef double real_t;
//#endif

//! Spatial algebra Vector, 6-d vector for spatial transformations.
typedef Eigen::Matrix< real_t, 6, 1> SpatialVector;

//! Spatial algebra Matrix, 6x6 matrix for spatial transformations.
typedef Eigen::Matrix< real_t, 6, 6> SpatialMatrix;

//! 9-dimensional vector for spatial transformations.
typedef Eigen::Matrix< real_t, 9, 1> Vector9r;

//! 9x9 matrix for spatial transformations.
typedef Eigen::Matrix< real_t, 9, 9> Matrix9r;

//! Derivative of Spatial Vector.
typedef Eigen::Matrix< real_t, 6, Eigen::Dynamic> D_SpatialVector;

//! Derivative of Spatial Matrix.
typedef Eigen::Matrix< real_t, Eigen::Dynamic, 6> D_SpatialMatrix;

//! Derivative of scalar real_t
typedef Eigen::Matrix< real_t, 1, Eigen::Dynamic> D_real_t;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t row major matrix. RowMajor defines the way the matrix is stored
 */
typedef  Eigen::Matrix<real_t,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXr;

typedef  Eigen::Matrix<real_t,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXrColMajor;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector.
 */
typedef  Eigen::Matrix<real_t,Eigen::Dynamic, 1> VectorXr;

typedef  Eigen::Matrix<real_t, 1, Eigen::Dynamic> RowVectorXr;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector of size 3x1.
 */
typedef Eigen::Matrix<real_t , 3, 1> Vector3r;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector of size 2x1.
 */
typedef Eigen::Matrix<real_t , 2, 1> Vector2r;

typedef Eigen::Matrix<real_t , 4, 1> Vector4r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 4x4.
 */
typedef Eigen::Matrix<real_t , 4 , 4> Matrix4r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 3x3.
 */
typedef Eigen::Matrix<real_t , 3 , 3> Matrix3r;


typedef Eigen::Matrix<real_t , 3 , 3 , Eigen::RowMajor> RowMatrix3r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 2x2.
 */
typedef Eigen::Matrix<real_t , 2 , 2> Matrix2r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 2x2.
 */
typedef Eigen::SparseMatrix< real_t > SparseMatrixXr; // Default ColMajor


typedef Eigen::Triplet< real_t > TripletXr;


typedef std::vector< TripletXr > TripletXrList;


typedef  Eigen::SparseVector< real_t > SparseVectorXr;


//! Type definition of Eigen::AngleAxis< Scalar > as Eigen::AngleAxis< real_t >
typedef Eigen::AngleAxis<real_t> AngleAxisr;

}

#endif // GEOMBD_TYPES_H
