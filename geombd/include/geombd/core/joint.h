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
 *	\file include/geombd/core/joint.h
 *	\author Gustavo Arechavaleta, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 * Joint class.
 */

#ifndef GEOMBD_CORE_JOINT_H
#define GEOMBD_CORE_JOINT_H

#include "geometry.h"
#include "types.h"
namespace geo{
//! Enum used to specify the type of the joint.
/*!
  \ingroup core_module
  A joint has an associated type. It can be revolute, prismatic or a more specialized type such as glue,
  if there is no specified type then it is default to Null. */
enum JointType
{
    kRevoluteJoint = 0,  /*!< The type of the joint is revolute. */
    kPrismaticJoint = 1, /*!< The type of the joint is prismatic. */
    kGlueJoint = 2,      /*!< The type of the joint is glue. */
    kNullJoint = -1      /*!< There is no information about the type. */  // TODO 03/2015: review the name.
};

//! Joint class
/*! \ingroup core_module
    A joint is one of the basic components of a robotic tree structure. It has upper and lower limits that correspond to the mechanical limits of the robot.
    The state of the joint is described by its position, velocity and acceleration.
    \sa Body, OperationalHandle
*/
class Joint
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Custom constructor. Initializes the jointId with the given id.
    /*! \param id the id of the joint.*/
    Joint(short int id);
    //! Custom constructor. Initializes the jointType with the given type.
    /*! \param jType the type of the joint.*/
    Joint(JointType jType);

    //! Default destructor
    ~Joint(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:
    //! Joint identifier number.
    short int jointId;

    //! Joint type (Revolute, Prismatic, etc.):
    JointType jointType;

    //! Joint position.
    real_t q;

    //! Joint velocity.
    real_t dq;

    //! Joint acceleration.
    real_t ddq;

    //! Joint torque.
    real_t tau;

    //! The axis in which the joint operates.
    Vector3r actionAxis;

    //! Joint upper limit in position.
    real_t upLimit;

    //! Joint lower limit in position.
    real_t lowLimit;

    //! Joint upper speed limit.
    real_t upSpeedLimit;

    //! Joint lower speed limit.
    real_t lowSpeedLimit;


    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! set upLimit to the specified value.
    void setUpLimit(real_t value){upLimit = value;}

    //! return the upLimit specified value.
    real_t getUpLimit(){return upLimit;}

    //! set lowLimit to the specified value.
    void setLowLimit(real_t value){lowLimit = value;}

    //! return the lowLimit specified value.
    real_t getLowLimit(){return lowLimit;}


    //! set upSpeedLimit to the specified value.
    void setUpSpeedLimit(real_t value){upSpeedLimit = value;}

    //! return the upSpeedLimit value.
    real_t getUpSpeedLimit(){return upSpeedLimit;}

    //! set lowSpeedLimit to the specified value.
    void setLowSpeedLimit(real_t value){lowSpeedLimit = value;}

    //! return the lowSpeedLimit value.
    real_t getLowSpeedLimit(){return lowSpeedLimit;}




    /*! Get the joint identifier.
         \param none
         \return the identification number.
         */
    short int getId(){ return jointId; }

    /*! Update joint configuration
         \param a real number.
         \return void
         */
    void setConfiguration(real_t value){ q = value; }

    /*! Get the joint configuration.
         \param none
         \return the configuration of the joint as geo::real_t
         */
    real_t getConfiguration(){ return q; }

    /*! Update joint velocity
         \param a real number.
         \return void
         */
    void setVelocity(real_t value){ dq = value; }

    /*! Get the joint velocity.
         \param none
         \return the velocity of the joint as geo::real_t
         */
    real_t getVelocity(){ return dq; }

    /*! Update joint acceleration
         \param a real number.
         \return void
         */
    void setAcceleration(real_t value){ ddq = value; }

    /*! Get the joint acceleration.
          \param none
          \return the acceleration of the joint as geo::real_t
         */
    real_t getAcceleration(){ return ddq; }

    /*! Update joint torque
         \param a real number.
         \return void
         */
    void setTorque(real_t value){ tau = value; }

    /*! Get the joint torque
          \param none
          \return the torque of the joint as geo::real_t
         */
    real_t getTorque(){ return tau; }

    /*! Set the type of the joint
         \param the joint Type: 0 for a revolute joint and 1 for prismatic one.
         \return void
         */
    void setJointType(JointType jType){ jointType = jType; }

    /*! Get the type of the joint.
         \param none
         \return the joint type as a jointType variable
         */
    JointType getJointType(){ return jointType; }


    /*! Set the action axis of the joint
         \param a vector with the action axis
         \return void
         */
    void setActionAxis(const Vector3r &axis){ actionAxis = axis; }



    /*! Get the action axis of the joint.
         \param none
         \return a coordinate instance with the action axis
         */
    const Vector3r& getActionAxis(){ return actionAxis; }





};
}   // end of namespace geo
#endif // GEOMBD_CORE_JOINT_H
