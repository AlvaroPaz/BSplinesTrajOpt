/*
 * MIT License
 *
 * Copyright (c) 2020 Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Carla Villanueva <carla.villanueva@cinvestav.edu.mx>,
 * Alvaro Paz <alvaro.paz@cinvestav.edu.mx>
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
 *	\file include/geombd/core/multibody.h
 *	\author Gustavo Arechavaleta, Carla Villanueva, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 * Multibody class.
 */

#ifndef HR_CORE_MULTIBODY_H
#define HR_CORE_MULTIBODY_H


#include "joint.h"
#include "dynamic_body.h"
#include "types.h"


namespace hr{
namespace core{

//! The type of the Jacobian.
/*! \ingroup core_module  */
enum JacobianType
{
    kGeometricJacobian = 0, /*!< The complete geometric jacobian. */
    kGeometricJacobian_Jv = 1, /*!< The linear velocity part of the jacobian. */
    kGeometricJacobian_Jw = 2 /*!< The angular velocity part of the jacobian. */
};



//! MultiBody class
/*! \ingroup core_module  \sa Body, DynamicBody, Joint, OperationalHandle */
class MultiBody
{
    //friend class World;

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Default constructor.
    MultiBody(){}

    //! Custom constructor.
    /*!
        \param id The identifier number of the multiBody.
        \param name The name of the multiBody.
      */
    MultiBody(short int id, const std::string &name)
    {
        multiBodyId = id;
        DoF = 0;
        this->name = name;
        this->gravity << 0, 0, 9.81, 0, 0, 0;
    }

    //! Default destructor.
    ~MultiBody(){}



    // --------------------------------------------
    // Members
    // --------------------------------------------
public:

    real_t totalMass;

protected:
    //! Multibody identifier.
    short int multiBodyId;

    //! Multibody name.
    std::string name;

    //! Multibody degrees of freedom.
    short int DoF;

    //! Multibody differentiation size
    int DiffSize;

    // lambda(i): The parent of body i
    //! Parents vector.
    std::vector<short int> parent;

    // kappa(i): The set of joints in the path between body i and the root
    //! Predecessor joints vector.
    std::vector< std::vector<short int> > predecessorJoints;

    // mu(i): The set of children of body i
    //! Children vector.
    std::vector< std::vector<short int>* > childrenBodies;

    // nu(i) the set of bodies in the subtree starting at body i.
    //! Subtree bodies vector.
    std::vector< std::vector<short int>* > subtreeBodies;

    // leafs of bodies
    //! Subtree bodies vector.
    std::vector<short int> subTree;

    // p(i) the predecessor body of joint i
    //! joint predecessor vector.
    std::vector<short int> jointPredecessor;

    // s(i) the succecessor body of joint i, i.e. the body i
    //! Joint successor vector.
    std::vector<short int> jointSuccessor;

    // X_T(i) the position of the frame of body i relative to the frame of body lambda(i)
    //! Relative body position vector.
    std::vector< Vector3r > relativeBodyPosition;

    //! The vector of pointers to bodies.
    std::vector<Body*> bodies;

    //! the vector of pointers to joints.
    std::vector< Joint* > joints;

    //! Spatial vector of gravity.
    SpatialVector gravity;

public:


    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! Get multiBody ID
         /*! \param none
         * \return short int identifying the multiBody.
         */
    inline short int getId(){ return multiBodyId; }

    //! Get multiBody total mass
         /*! \param none
         * \return real_t
         */
    real_t getTotalMass(){ return totalMass; }

    //! Set the multiBody total mass
         /*! \param real_t
         * \return void
         */
    inline void setTotalMass( real_t totalMass ){ this->totalMass = totalMass; }

    //! Get the multiBody name
         /*! \param none
         * \return the multibody name as a string.
         */
    const std::string& getName(){ return this->name; }

    //! Get the number of degrees of freedom
         /*! \param none
         * \return the number of degrees of freedom
         */
    short int getDoF(){ return DoF; }

    //! Get the differentiation size
         /*! \param none
         * \return the number partials size
         */
    int getDifferentiationSize(){ return DiffSize; }

    //! Set the number of degrees of freedom
         /*! \param the number of degrees of freedom
         * \return void
         */
    inline void setDoF( short int dof ){ DoF = dof; }

    //! Set the differentiation size
         /*! \param the differentiation size
         * \return void
         */
    inline void setDifferentiationSize( int DiffSizeNew ){ DiffSize = DiffSizeNew; }

    //! Get the configuration of the multibody systemore
         /*! \param none
         * \return a vector of dimension n of real_ts
         */
    VectorXr getConfiguration();

    //! Get the gravity of multibody
         /*! \param none
         * \return spatial vector
         */
    SpatialVector getGravity(){ return gravity; }

    //! Set the configuration of the multibody system
         /*! \param a vector of dimension n of real_ts
         * \return void
         */
    void setConfiguration( const VectorXr &q);

    //! Get the multibody generalized velocity vector
         /*! \param none
         * \return a vector of dimension n of real_ts
         */
    VectorXr getGeneralizedVelocity();

    //! Set the generalized velocity of the multibody
         /*! \param a vector of dimension n of real_ts
         * \return void
         */
    void setGeneralizedVelocity( const VectorXr &dq);

    //! Get the multibody generalized acceleration vector
         /*! \param none
         * \return a vector of dimension n of real_ts
         */
    VectorXr getGeneralizedAcceleration();

    //! Set the generalized acceleration of the multibody
         /*! \param a vector of dimension n of real_ts
         * \return void
         */
    void setGeneralizedAcceleration( const VectorXr &ddq);

    //! Get the multibody generalized torque vector
         /*! \param none
         * \return a vector of dimension n of real_ts
         */
    VectorXr getGeneralizedTorques();

    //! Set the generalized torques of the multibody
         /*! \param a vector of dimension n of real_ts
         * \return void
         */
    void setGeneralizedTorques( const VectorXr &torques);

    //! Load the multiBody tree with URDF File
         /*! \param a string with the name of the xml file
         * \param a string with the name of the URDFxml file
         * \return bool status
         */
    bool loadTreeURDF( const std::string &xmlFileName, const std::string &URDFxmlFileName );

    //! Create a new body in the multibody system
         /*! \param the name of the new body
         * \param the xmlFile name of the new body
         * \param the Id of the new body's parent
         * \return bool status
         */
    bool newBody(std::string newBodyName, const std::string &xmlFileName, short int parentBodyId);

    //! Create a new body in the multibody system
         /*! \param the name of the new body
         * \param the xmlFile name of the new body
         * \param the Id of the new body's parent
         * \param the body's relative position
         * \return bool status
         */
    bool newBodyURDF(std::string newBodyName, const std::string &xmlFileName, short int parentBodyId, const Vector3r &relativePosition);

    //! Get the children Bodies of body i.
         /*! \param the id of the body
         * \return the body children vector
         */
    std::vector<short int> getBodyChildren( short int bodyId){ return *childrenBodies.at( bodyId ); }

    //! Get the vector of subTree
         /*! \param none
         * \return the vector of pointers to joint successors
         */
    std::vector<short int> getSubTreeBody(short int bodyId){ return *subtreeBodies.at( bodyId ); }
    std::vector<short int> getSubTree(){ return subTree; }

    //! Get the id of parent body of body i.
         /*! \param the id of the body
         * \return the id of the parent body
         */
    int getParentId(short int bodyId){ return parent.at( bodyId ); }

    //! Get the position of the body i relative to its parent.
         /*! \param the id of the body
         * \return the relative position vector
         */
    Vector3r getRelativeBodyPosition(short int bodyId){ return relativeBodyPosition.at( bodyId ); }

    //! Get the id of the predecessor body of joint i.
         /*! \param joint id
         * \return the id of the predecessor body
         */
    int getPredecessorId( short int jointId ) { return jointPredecessor.at(jointId-1); }

    //! Get the id of the successor body of joint i.
         /*! \param joint id
         * \return the id of the successor body
         */
    int getSuccessorId( short int jointId ) { return jointSuccessor.at(jointId-1); }

    //! Get the predecessor joints vector of body i.
         /*! \param body id
         * \return the predecessor joints vector
         */
    std::vector<short int> getPredecessorJoints( short int bodyId ) { return predecessorJoints.at(bodyId); }

    //! Create a new joint in the multibody system
         /*! \param string with the joint type
         * \param the action axis
         * \param the predecessor body name
         * \param the successor body name
         * \param the xmlFile name
         * \return bool status
         */
    bool newJoint(std::string JointType, const Vector3r &actionAxis, std::string jointPredecessor,
                   std::string jointSuccessor, real_t ulim, real_t llim, const std::string &xmlFileName );

    //! Create a new joint in the multibody system
         /*! \param string with the joint type
         * \param the action axis
         * \param the predecessor body name
         * \param the successor body name
         * \param the xmlFile name
         * \return bool status
         */
    bool newJoint(std::string JointType, const Vector3r &actionAxis, std::string jointPredecessor,
                   std::string jointSuccessor, real_t ulim, real_t llim, real_t maxSpeed, const std::string &xmlFileName );

    //! Create a new joint in the multibody system
         /*! \param string with the joint type
         * \param the action axis
         * \param the predecessor body name
         * \param the successor body name
         * \param the xmlFile name
         * \param the body's relative position
         * \return bool status
         */
    bool newJointURDF(std::string JointType, const Vector3r &actionAxis, std::string jointPredecessor,
                   std::string jointSuccessor, real_t ulim, real_t llim, real_t maxSpeed, const std::string &xmlFileName, const Vector3r &relativePosition);

    //! Clear vectors of bodies, joints, joint predecessors and successors
         /*! \param none
         * \return void
         */
    void clear();

    //! Print the multibody structure
         /*! \param none
         * \return void
         */
    void printStructure();

    //! Get the skewMatix from the given vector
         /*! \param e 3d vector
         * \return the corresponding MatrixXr
         */
    Matrix3r getSkewMatrix(const Vector3r &vector);

    //! Get the size of vector of Bodies from the multibody system.
         /*! \param none
         * \return the number of bodies composing the multibody system.
         */
    int getNumberOfBodies(){ return bodies.size(); }

    //! Evaluate if there exist collisions between bodies
         /*! \param none
         * \return void
         */
    void autoCollisionQuery();

    //! Get the upper limit of a joint from the multibody system.
         /*! \param joint id
         * \return the joint upper limit
         */
    real_t getJointUpLimit( short int jointId ){ return this->getJoint(jointId)->getUpLimit(); }

    //! Get the lower limit of a joint from the multibody system.
         /*! \param joint id
         * \return the joint lower limit
         */
    real_t getJointLowLimit( short int jointId ){ return this->getJoint(jointId)->getLowLimit(); }

    //! Get the upper limits of the joints from the multibody system.
         /*!
         * \return the joint upper limits in a VectorXr
         */
    VectorXr getJointUpLimits();

    //! Get the lower limits of the joints from the multibody system.
         /*!
         * \return the joint lower limits in a VectorXr
         */
    VectorXr getJointLowLimits();


    //! Get the upper speed limits of the joints from the multibody system.
         /*!
         * \return the joint upper limits in a VectorXr
         */
    VectorXr getJointUpSpeedLimits();

    //! Get the lower speed limits of the joints from the multibody system.
         /*!
         * \return the joint lower limits in a VectorXr
         */
    VectorXr getJointLowSpeedLimits();

    //! Set the CoM handles to every DynamicBody
    virtual void setCoMHandles()= 0;

    //! Get the position of a given body from the multibody system.
         /*! \param body id
         * \return the current 3d position of the corresponding body
         */
    Vector3r getBodyPosition(const short int &bodyId);

    //! Get the orientation of a given body from the multibody system.
         /*! \param body id
         * \return the rotation matrix attached to the body
         */
    Matrix3r getBodyOrientation(const short int &bodyId);

    //! Get the geometry of a given body from the multibody system.
         /*! \param body id
         * \return the geometry object attached to the body
         */
    Geometry* getBodyGeometry(const short int &bodyId){
        return this->getBody(bodyId)->getGeometry();
    }

    //! Get the position of a given operational Handle attached to a body from the multibody system.
         /*! \param body id
         * \param handle id
         * \return the current 3d position of the corresponding handle
         */
    Vector3r getOperationalHandlePosition(const short int &bodyId, const short int &handleId);

    //! Get the orientation of a given operational Handle attached to a body from the multibody system.
         /*! \param body id
         * \param handle id
         * \return the current orientation of the corresponding handle
         */
    Matrix3r getOperationalHandleOrientation(const short int &bodyId, const short int &handleId);

    //! Get the configuration of a given joint from the multibody system.
         /*! \param joint id
         * \return the current value of the joint configuration
         */
    real_t getJointConfiguration(const short int &jointId){
        return this->getJoint(jointId)->getConfiguration();
    }

    //! Set the configuration of a given joint in the multibody system.
         /*! \param joint id
         * \param configuration value of jointMULTIBODY_H
         * \return void
         */
    void setJointConfiguration(const short int &jointId, real_t value){
        this->getJoint(jointId)->setConfiguration(value);
    }

    //! Get the joint type between two bodies in the multibody system.
         /*! \param joint id
         * \return the type of joint
         */
    JointType getJointType(const short int &jointId){
        return this->getJoint(jointId)->getJointType();
    }

    //! Add an operational handle to a body in the multibody system.
         /*! \param body id
         * \param the position of the handle relative to the body frame
         * \return the id of the new operational handle
         */
    int addOperationalHandle(const short int &bodyId, const Vector3r &handlePosition){
        return (this->getBody(bodyId))->addOperationalHandle(handlePosition);
    }

    //! Get the pointer to the vector of Bodies from the multibody system.
         /*! \param none
         * \return the vector of pointers to bodies
         */
    std::vector<Body*> getPubBodies(){ return bodies; }

protected:

    //! Find and get a Body from the multibody system.
         /*! \param body id
         * \return a pointer to the corresponding body
         */
    Body* getBody( short int id );    //TODO 03/2015: Change to DynamicBody

    //! Find and get a Body from the multibody system.
         /*! \param body name
         * \return a pointer to the corresponding body
         */
    Body* getBody( std::string bodyName );   //TODO 03/2015: Change to DynamicBody

    //! Get the pointer to the vector of Bodies from the multibody system.
         /*! \param none
         * \return the vector of pointers to bodies
         */
    std::vector<Body*>* getBodies(){ return &bodies; }   //TODO 03/2015: Change to DynamicBody

    //! Find and get a Joint from the multibody system.
         /*! \param joint id
         * \return a pointer to the corresponding joint
         */
    Joint* getJoint( short int jointId );

    //! Get the vector of joints from the multibody system.
         /*! \param none
         * \return the vector of pointers to joints
         */
    std::vector<Joint*>* getJoints(){ return &joints; }

    //! Get the vector of joint predecessors from the multibody system.
         /*! \param none
         * \return the vector of pointers to joint predecessors
         */
    std::vector<short int>* getJointPredecessors(){ return &jointPredecessor; }

    //! Get the vector of joint successors from the multibody system.
         /*! \param none
         * \return the vector of pointers to joint successors
         */
    std::vector<short int>* getJointSuccessors(){ return &jointSuccessor; }

    //! Fill the vector of subtree bodies
         /*! \param none
         * \return void
         */
    void fillSubtreeBodies();

    //! Compute the pose of body i from the parent pose.
         /*! \param a pointer to the body
         * \param the position vector of the body parent
         * \param the orientation matrix of the body parent
         * \return void
         */
    void computeBodyPoseFromParent(Body* body, const Vector3r &parentPosition,
                                   const Matrix3r &parentOrientation);

    bool collisionQueryAllowed(Body* body1, Body* body2);

    //! Check if two bodies are in collision.
    bool inCollision(Body* body1, Body* body2);

};

//! MultiBody for NAO Humanoids derived from MultiBody class.
class MultiBodyNAO : public MultiBody
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Custom constructor
    MultiBodyNAO(short int id, const std::string &name) : MultiBody(id,name) {}


    void setCoMHandles();


};

//! MultiBody for HRP-2 Humanoids derived from MultiBody class.
class MultiBodyHRP2 : public MultiBody
{

public:
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
    MultiBodyHRP2(short int id, const std::string &name) : MultiBody(id,name) {}

    void setCoMHandles();

};



} //end of namespace core
} //end of namespace hr
#endif // HR_CORE_MULTIBODY_H
