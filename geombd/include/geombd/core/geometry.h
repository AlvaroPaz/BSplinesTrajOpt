/*
 *
 * Copyright (C) 2015
 * Julio Jarquin <jjarquin.ct@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx> , Gerardo Jarquin
 * CINVESTAV - Saltillo Campus
 *
 * This file is part of OpenHRC
 * OpenHRC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * OpenHRC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/**
 *	\file include/geombd/core/geometry.h
 *	\author Gerardo Jarquin, Gustavo Arechavaleta, Julio Jarquin.
 *	\version 1.0
 *	\date 2015
 *
 * Geometry classes (Triangle, Geometry Piece, Geometry).
 * Geometry classes are implemented using Vector3f for points.
 */

#ifndef HR_CORE_GEOMETRY_H
#define HR_CORE_GEOMETRY_H

#include <vector>
//#include <cstdlib>
#include <string>
//#include "PQP/PQP.h"
#include "Eigen/Dense"
#include "string_util.h"
#include <iostream>

namespace hr
{
namespace core{

//! Triangle class.
/*! \ingroup core_module
    A simple class to store the information of a simple 3D triangle.
    A triangle is composed of three vertex (3D points), vertexA, vertexB, vertexC.
    Also, it is important to calculate its normal vector.
    \sa GeometryPiece Geometry
  */
class Triangle
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Default constructor.
    Triangle(){}

    //! Default destructor.
    ~Triangle(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:
    //! Vertex A of the triangle.
    Eigen::Vector3f vertexA;
    //! Vertex B of the triangle.
    Eigen::Vector3f vertexB;
    //! Vertex C of the triangle.
    Eigen::Vector3f vertexC;
    //! Normal vector of the triangle.
    Eigen::Vector3f normal;

    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! Set the three vertex of the triangle and calculate its normal vector.
    void setVertices(const Eigen::Vector3f &vertexA,
                     const Eigen::Vector3f &vertexB,
                     const Eigen::Vector3f &vertexC);

    //! Returns the vertexA of the triangle.
    Eigen::Vector3f* getVertexA(){ return &vertexA; }

    //! Returns the vertexB of the triangle.
    Eigen::Vector3f* getVertexB(){ return &vertexB; }

    //! Returns the vertexC of the triangle.
    Eigen::Vector3f* getVertexC(){ return &vertexC; }

    //! Returns the normal vector of the triangle.
    Eigen::Vector3f* getNormal(){ return &normal; }

    //! Print method for debugging, consider deleting it.
    void printA(){std::cout << vertexA;} // TODO 03/2015: Debug method?
};

//! \ingroup core_module GeometryPiece class.
/*! A simple class to store the information of a 3D geometric body piece.
  \sa Geometry Triangle
  */
class GeometryPiece
{
public:
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
    //! Default constructor.
    GeometryPiece(){}

    //! Default destructor.
    ~GeometryPiece(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------
    //! RGB Color of the body piece
    Eigen::Vector3f colorRGB;

    //! Vector of vertices of the body piece
    std::vector<Eigen::Vector3f> vertices;

    //! Vector of triangles of the body piece
    std::vector<Triangle> triangles;

    // --------------------------------------------
    // Methods
    // --------------------------------------------

    //! Set the vertices of the body piece.
    /*!
      \param verticesString the string containting the body vertices.
      */
    void setVertices(std::string verticesString,
                     bool rotFlag, const Eigen::AngleAxisf &rotation,
                     bool scaleFlag, const Eigen::Vector3f &scale,
                     bool centerFlag, const Eigen::Vector3f &center,
                     bool transFlag, const Eigen::Vector3f &translation);

    //! Set the vertices of the body pieces.
    /*!
      \param vertexIndicesString the string containting the vertices to form the triangles.
      */
    void setTriangles(std::string vertexIndicesString);

    //! Get a triangle from the triangle vector.
    /*!
      \param triangleId the id of the triangle.
      \return A pointer to the triangle.
      */
    Triangle* getTriangle( int triangleId ){ return &triangles.at( triangleId ); }

    void printTriangle(int tId){ std::cout << "Tri " << tId << std::endl; triangles.at(tId).printA(); }  // TODO 03/2015: Debug method?
};

//! \ingroup core_module Geometry class.
/*! A  class to store a vector of 3D geometric body pieces.
    \sa GeometryPiece Triangle
  */
class Geometry
{
public:
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
    //! Default constructor.
    Geometry(){ collisionFlag = false; }
    //! Default destructor.
    ~Geometry(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------
    //! Vector of geometric pieces
    std::vector<GeometryPiece*> geomPieces;  //TODO: use smart pointer

    // --------------------------------------------
    // Methods
    // --------------------------------------------
    //! Set the collisionStatus of the geometry
    /*! \param colFlag whether the geometry is in collision or not  */
    void setCollisionStatus( bool colFlag){collisionFlag = colFlag; }
    //! Get the collisionStatus of the geometry
    /*! \return whether the geometry is in collision or not*/
    bool getCollisionStatus(){ return collisionFlag; }

private:
    bool collisionFlag;
};
} // end of namespace core
} // end of namespace hr
#endif // HR_CORE_GEOMETRY_H
