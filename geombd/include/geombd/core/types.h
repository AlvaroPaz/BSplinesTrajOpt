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
 *	\file include/openhrc/core/types.h
 *	\author Gustavo Arechavaleta, Alvaro Paz
 *	\version 1.0
 *	\date 2020
 *
 *	Definition of Types
 */

#ifndef HR_CORE_TYPES_H
#define HR_CORE_TYPES_H

#include "geombd/types.h"


namespace hr{

//! Namespace core. Namespace of the core Module
/*!
    More information about core.
*/
typedef Eigen::Matrix<real_t, 3 , Eigen::Dynamic , Eigen::RowMajor> Matrix3Xr;

namespace core{

    //! Type definition of Eigen::AngleAxis< Scalar > as Eigen::AngleAxis< real_t >
    typedef Eigen::AngleAxis<real_t> AngleAxisr;


//! Namespace nao. Namespace for nao settings
/*!
    More information about nao.
*/
namespace nao{

enum body{

            kLeftFoot = 24,
            kRightFoot = 30,
            kLeftHand = 13,
            kRightHand = 18,
            kTorso = 6,
            kHead = 8
};

enum Pose{
    kStandZero,
    kStandInit
};


} // end of namespace nao

} // end of namespace core

} // end of namespace hr


#endif // HR_CORE_TYPES_H
