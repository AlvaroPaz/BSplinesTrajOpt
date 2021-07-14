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
 *	\file include/geombd/core/operational_handle.h
 *	\author Gustavo Arechavaleta, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 * OperationalHandle class.
 */

#ifndef HR_CORE_OPERATIONAL_HANDLE_H
#define HR_CORE_OPERATIONAL_HANDLE_H

#include <string>
#include "types.h"

namespace hr {
namespace core {

//! OperationalHandle class definition
/*! \ingroup core_module
    This class represents a coordinate frame that can be attached to
    a given Body or and object of the environment. \sa Body, Joint
    */
class OperationalHandle {

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param name The name of the operational handle.
      */
    OperationalHandle(const unsigned int id, const std::string &name = "End-effector");

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param pos The position of the operational handle, relative to the associated body.
        \param name The name of the operational handle.

      */
    OperationalHandle(const unsigned int id, const Vector3r &pos,
                      const std::string &name = "End-effector");

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param name The name of the operational handle.
        \param pos The position of the operational handle, relative to the associated body
        \param rot The corresponding rotation (Rx,Ry,Rz) of the operational handle.
      */
    OperationalHandle(const unsigned int id, const std::string &name,
                      const Vector3r &pos, const Matrix3r &rot);

    //! Default destructor.
    ~OperationalHandle();

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:
    //! OperationalHandle identifier
    unsigned int id;

    //! OperationalHandle name
    std::string name;

    //! OperationalHandle position (x,y,z) relative to the associated body
    Vector3r position;

    //! OperationalHandle rotation matrix (Rx,Ry,Rz)
    /*! This matrix represents the orientation relative to the associated body
     */
    Matrix3r orientation;

    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! Get OperationalHandle ID
     /*! \param none
       \returns short int identifying the opHandle.
     */
    short int getId(){ return id; }

    //! Set OperationalHandle ID
     /*! \param short int identifying the opHandle
         \returns void
     */
    void setId(const int &id){ this->id = id; }

    //! Get OperationalHandle name
     /*! \param name
      \returns string with the opHandle name.
     */
    std::string getName(){ return this->name; }

    //!  Set OperationalHandle name
     /*! \param string with the opHandle name.
         \returns void
     */
    void setName(const std::string &name){ this->name = name; }

    //! Get the position of the OperationalHandle.
     /*! \param none
      \returns 3-dimensional real_t precision array (3D position)
     */
    Vector3r getPosition(){ return position; }

    //! Set the position of the OperationalHandle
    /*! \param a 3D point of real_t precision numbers.
      \returns void
     */
    void setPosition(const Vector3r &position){ this->position = position; }

    //! Get the orientation of the OperationalHandle.
     /*! \param none
       \returns a 3x3 real_t precision matrix R in SO(3)
     */
    Matrix3r getOrientation(){ return orientation; }

    //! \brief Set the orientation of the OperationalHandle
     /*! \param a 3x3 real_t precision matrix R in SO(3).
       \returns void
     */
    void setOrientation(const Matrix3r &orientation){ this->orientation = orientation; }
};

} // end of namespace core
} // end of namespace hr

#endif // OPERATIONALHANDLE_H_
