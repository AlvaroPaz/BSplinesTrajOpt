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
 *	\file include/openhrc/core/world.h
 *	\author Gustavo Arechavaleta, Carla Villanueva, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 *	World class.
 */
#ifndef HR_CORE_WORLD_H
#define HR_CORE_WORLD_H

#include <memory>
#include "multibody.h"

namespace hr{
namespace core{

//! Enum used to specify the robot type. A robot can be NAO, HRP2, YOUBOT.
/*! \ingroup core_module  */
enum RobotType
{
    kNAO = 0, /*!< The multibody is a humanoid robot NAO */
    kHRP2 = 1, /*!< The multibody is a humanoid robot HRP2 */
    kYOUBOT = 2, /*!< The multibody is a YOUBOT */
    kPA10 = 3 /*!< The multibody is a PA10 */
};

//! World class. A world can hold many bodies and multibodies, and test the collision between the bodies.
/*! \ingroup core_module  */
class World
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Default constructor.
    World(){}   

    // --------------------------------------------
    // Members
    // --------------------------------------------
private:

    //! Vector of pointers to Multibody.
    /*! This multibodies are the robots in the world. */
    std::vector< std::shared_ptr< MultiBody > > robots;

    //! Vector of pointers to Body.
    /*! This bodies can be considered obstacles in the world. */
    std::vector< Body* > obstacles;

    //! Iterator for the vector of pointers to robots.
    std::vector< std::shared_ptr< MultiBody > >::iterator robots_iter;

    //! Iterator for the vector of pointers to bodies.
    std::vector< Body* >::iterator robotBodies_iter;


    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:


//    //! Vector of pointers to Body.
//    /*! These bodies can be considered obstacles in the world. */
//    /*! \param sFileName The XML file name of the robot geometry specficiation.
//        \param robotID The robot identifier.
//        \param robotType The robot type (NAO, HRP2, etc.).
//        \return void. */
//    void loadMultiBody(std::string sFileName, int robotID, const RobotType robotType);

    //! Return a robot.
    /*! \param id The robot identifier.
        \return A pointer to the Multibody with the corresponding name.*/
    std::shared_ptr< MultiBody > getRobot(short int id);

    //! Return a robot.
    /*! \param robotName The name of the robot.
        \return A pointer to the Multibody with the corresponding name.*/
    std::shared_ptr< MultiBody > getRobot(std::string robotName);

    //! Return the vector of robots.
    std::vector< std::shared_ptr< MultiBody > >* getRobotsVector(){return &robots;}

    //! The the collision of the robots and bodies in the world.
    void collisionQuery();

    //! Vector of pointers to Body.
    /*! These bodies can be considered obstacles in the world. */
    /*! \param sFileName The URDF_XML file name of the robot specficiation.
        \param robotID The robot identifier.
        \param robotType The robot type (NAO, HRP2, etc.).
        \return void. */
    void loadMultiBodyURDF(std::string sFileName, int robotID, const RobotType robotType);


    //bool inCollision(Body* body1, Body* body2);
    //bool collisionQueryAllowed(Body* body1, int robot1_ID, Body* body2, int robot2_ID);
    //void partialCollisionQuery(Body* body, int robotID);

};
}   // end of namespace core
}   // end of namesapce hr

#endif // HR_CORE_WORLD_H

