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
 *	\file include/openhrc/core/string_util.h
 *	\author Gustavo Arechavaleta, Alvaro Paz.
 *	\version 1.0
 *	\date 2020
 *
 * Inline string manipulation methods.
 */

#ifndef HR_CORE_STRING_UTIL_H
#define HR_CORE_STRING_UTIL_H

#include <string>
#include <vector>
#include "Eigen/Dense"
#include "types.h"

//using namespace std;

namespace hr
{
namespace core
{
inline int Split(const std::string& input, const std::string& delimiter,
                 std::vector<std::string>& results, bool includeEmpties)
{
    std::vector<int> positions;
    int offset;
    int numFound = 0;
    int iPos = 0;
    int newPos = -1;
    int sizeS2 = (int)delimiter.size();
    int isize = (int)input.size();

    if( isize == 0 || sizeS2 == 0 )
        return 0;
    newPos = input.find (delimiter, 0);
    if( newPos < 0 )
        return 0;
    while( newPos >= iPos ) {
        positions.push_back(newPos);
        iPos = newPos;
        newPos = input.find (delimiter, iPos+sizeS2);
    }
    for ( int i=0; i <= (int)positions.size(); ++i ) {
        std::string s("");
        offset = 0;
        if( i == 0 ) {
            s = input.substr( i, positions[i] );
            offset = positions[0] + sizeS2;
        }
        else
            offset = positions[i-1] + sizeS2;
        if ( offset < isize ) {
            if( i == (int)positions.size() )
                s = input.substr(offset);
            else if( i > 0 )
                s = input.substr( positions[i-1] + sizeS2, positions[i] - positions[i-1] - sizeS2 );
        }
        if( includeEmpties || ( s.size() > 0 ) ) {
            numFound++;
            results.push_back(s);
        }
    }
    return numFound;
}

inline void TrimLeft(std::string& str, const char* chars2remove)
{
    if (!str.empty()) {
        std::string::size_type pos = str.find_first_not_of(chars2remove);
        if (pos != std::string::npos)
            str.erase(0,pos);
        else
            str.erase( str.begin() , str.end() ); // make empty
    }
}

inline void TrimRight(std::string& str, const char* chars2remove)
{
    if (!str.empty()){
        std::string::size_type pos = str.find_last_not_of(chars2remove);
        if (pos != std::string::npos)
            str.erase(pos+1);
        else
            str.erase( str.begin() , str.end() ); // make empty
    }
}

inline void Trim(std::string& str, const char* chars2remove)
{
    if (!str.empty()) {
        TrimRight( str, chars2remove );
        TrimLeft( str, chars2remove );
    }
}

inline int stringToVector(std::vector<std::string> &resultVector, std::string inputString,
                          const std::string& delimiter)
{
    Trim(inputString, " ,\t\r\n");
    if (inputString != "") {
        int vectorSize = Split(inputString, delimiter, resultVector, false);
        return vectorSize;
    }
    return 0;
}

inline Eigen::VectorXf stringToVectorXf(bool &status, std::string inputString,
                                       const std::string& delimiter)
{
    Eigen::VectorXf outputVector;
    std::vector<std::string> resultVector;

    short int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize == 0) {
        status = false;
        return outputVector;
    }
    outputVector.resize( vectorSize );
    for (short int i=0; i<vectorSize; i++)
        outputVector(i) = atof( resultVector.at(i).c_str() );
    status = true;
    return outputVector;
}

inline Eigen::Vector3f stringToVector3f(bool &status, std::string inputString, const std::string& delimiter)
{
    Eigen::Vector3f outputVector;
    std::vector<std::string> resultVector;

    int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize != 3) {
        status = false;
        return outputVector.setZero();
    }
    for (short int i=0; i<3; i++)
        outputVector(i) = atof( resultVector.at(i).c_str() );
    status = true;
    return outputVector;
}

inline Vector3r stringToVector3r(bool &status, std::string inputString, const std::string& delimiter)
{
    Vector3r outputVector;
    std::vector<std::string> resultVector;

    int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize != 3) {
        status = false;
        return outputVector.setZero();
    }
    for (short int i=0; i<3; i++)
        outputVector(i) = atof( resultVector.at(i).c_str() );
    status = true;
    return outputVector;
}

inline VectorXr stringToVectorXr(bool &status, std::string inputString,
                                       const std::string& delimiter)
{
    VectorXr outputVector;
    std::vector<std::string> resultVector;

    short int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize == 0) {
        status = false;
        return outputVector;
    }
    outputVector.resize( vectorSize );
    for (short int i=0; i<vectorSize; i++)
        outputVector(i) = atof( resultVector.at(i).c_str() );
    status = true;
    return outputVector;
}


inline Eigen::AngleAxisf stringToAngleAxisf(bool &status, std::string inputString, const std::string& delimiter)
{
    Eigen::AngleAxisf outputAngleAxis;
    std::vector<std::string> resultVector;
    Eigen::Vector3f axis;

    int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize != 4) {
        status = false;
        return outputAngleAxis.Identity();
    }
    for (short int i=0; i<3; i++)
        axis(i) = atof( resultVector.at(i).c_str() );
    outputAngleAxis.axis() = axis;
    outputAngleAxis.angle() = atof( resultVector.at(3).c_str() );
    status = true;
    return outputAngleAxis;
}

inline AngleAxisr stringToAngleAxisr(bool &status, std::string inputString, const std::string& delimiter)
{
    AngleAxisr outputAngleAxis;
    std::vector<std::string> resultVector;
    Vector3r axis;

    int vectorSize = stringToVector(resultVector, inputString, delimiter);
    if (vectorSize != 4) {
        status = false;
        return outputAngleAxis.Identity();
    }
    for (short int i=0; i<3; i++)
        axis(i) = atof( resultVector.at(i).c_str() );
    outputAngleAxis.axis() = axis;
    outputAngleAxis.angle() = atof( resultVector.at(3).c_str() );
    status = true;
    return outputAngleAxis;
}


} // end of namespace core
} // end of namespace hr
#endif //HR_CORE_STRING_UTIL_H
