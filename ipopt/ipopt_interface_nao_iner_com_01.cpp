// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file ipopt_interface_01.cpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Source class to include the Trajectory Optimization Problem into the Ipopt Solver
 */

#include "geombd/core.h"
#include "geombd/dynamics.h"
#include "geombd/trajectoryOptimization.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <memory>
#include <fstream>
#include <sstream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iomanip>

//! Coppelia headers files
extern "C" {
    #include "remoteApi/extApi.h"
    #include "remoteApi/extApiPlatform.h"
    #include "remoteApi/simConst.h"
//    #include "remoteApi/extApiInternal.h"
//    #include "remoteApi/shared_memory.h"
}

using std::cout;
using std::endl;
using namespace Eigen;


#include <cassert>
#include <iostream>

#include "ipopt_interface_nao_iner_com_01.hpp"

typedef Eigen::Map<const geo::VectorXr> MapVec;

std::ofstream myfile_1, myfile_2, myfile_3;

// using namespace Ipopt;

//! constructor
PracticeNLP::PracticeNLP( std::shared_ptr< geo::MultiBody > robot, std::shared_ptr< geo::robotSettingsTrajectoryOptimization > robotSettings ) : robotNonlinearProblem{nullptr} {

    this->robot = robot;
    this->robotSettings = robotSettings;
    this->StackConstraints = robotSettings->StackConstraints;
    this->nDoF = robotSettings->n;
    this->numberControlPoints = robotSettings->numberControlPoints;
    this->numberPartitions = robotSettings->numberPartitions;
    this->si = robotSettings->si;
    this->sf = robotSettings->sf;

    // Retrieve time vector
    this->S = robotSettings->S;

    // Motion object instance
    this->robotNonlinearProblem = new geo::DirectCollocation( robot, robotSettings );

    this->numberConstraints = robotSettings->numberConstraints;

    // Build basis function for B-Splines
    robotNonlinearProblem->buildBasisFunctions(  );

    // Allocating control points
    this->controlPoints = geo::VectorXr::Zero( numberControlPoints*nDoF, 1 );

    // Set initial configuration
    this->initialConfiguration = robotSettings->initialConfiguration;

    // Set final configuration
    this->finalConfiguration = robotSettings->finalConfiguration;

    // Set initial generalized velocity
    this->initialGeneralizedVelocity = robotSettings->initialGeneralizedVelocity;

    // Set final generalized velocity
    this->finalGeneralizedVelocity = robotSettings->finalGeneralizedVelocity;

    myfile_1.open ("qData.txt");
    myfile_2.open ("dqData.txt");
    myfile_3.open ("ddqData.txt");

}

//! destructor
PracticeNLP::~PracticeNLP()
{
    if(robotNonlinearProblem){delete robotNonlinearProblem;}
}

//! returns the size of the problem
bool PracticeNLP::get_nlp_info(Ipopt::Index& n, 
                               Ipopt::Index& m,
                               Ipopt::Index& nnz_jac_g,
                               Ipopt::Index& nnz_h_lag,
                               Ipopt::TNLP::IndexStyleEnum& index_style)
{

    //! Update settings
    robotNonlinearProblem->updateSettings();

    //! size of decision variable
    n = controlPoints.size();

    //! total numer of equality and inequality constraint
    m = numberConstraints;

    //! nonzeros in the constraints Jacobian
    nnz_jac_g = m*n;

    //! for the Hessian, we only need the lower left corner (since it is symmetric)
    nnz_h_lag = (n+1)*n/2;  // kid Gauss formula

    //! use the C style indexing (0-based)
    index_style = Ipopt::TNLP::C_STYLE;

    return true;
}

//! returns the variable bounds
bool PracticeNLP::get_bounds_info(Ipopt::Index n,
                                  Ipopt::Number* x_l,
                                  Ipopt::Number* x_u,
                                  Ipopt::Index m,
                                  Ipopt::Number* g_l,
                                  Ipopt::Number* g_u)
{

//    assert(n == controlPoints.size());
//    assert(m == 4*nDoF);

    //! upper joint limits
    Map<geo::VectorXr>( x_u, n, 1 ) = robot->getJointUpLimits( ).replicate( numberControlPoints, 1 );

    //! lower joint limits
    Map<geo::VectorXr>( x_l, n, 1 ) = robot->getJointLowLimits( ).replicate( numberControlPoints, 1 );

    //! all constraints are assumed bilateral
    geo::VectorXr constraintsUP(m,1);
    geo::VectorXr constraintsLOW(m,1);

    constraintsUP.setOnes();
    constraintsLOW.setOnes();

    //! This is customizable too
    geo::real_t qiBound_u = 1e-3;             //! Position 1e-3
    geo::real_t qfBound_u = 1e-3;
    geo::real_t dqiBound_u = 1e-3;            //! Velocity 1e-3
    geo::real_t dqfBound_u = 1e-3;
    geo::real_t pelvisBound_u = 1e-3;         //! Pelvis symmetry 1e-10
    geo::real_t comBound_u = 0.026;           //! CoM    0.026
    geo::real_t muBound_u = 1e-2;             //! Centroidal momentum 1e-2

    geo::real_t qiBound_l = -1e-3;            //! Position -1e-3
    geo::real_t qfBound_l = -1e-3;
    geo::real_t dqiBound_l = -1e-3;           //! Velocity -1e-3
    geo::real_t dqfBound_l = -1e-3;
    geo::real_t pelvisBound_l = -1e-3;        //! Pelvis symmetry -1e-10
    geo::real_t comBound_l = -0.026;          //! CoM    -0.026
    geo::real_t muBound_l = -1e-2;            //! Centroidal momentum -1e-2


    //! Fill up constraints
    int innerIndex = 0;
    for(geo::ConstraintsStack::iterator it = StackConstraints.begin(); it != StackConstraints.end(); ++it) {
        switch ( *it ) {
        case geo::constraint_initialConfiguration: {
            constraintsUP.segment(innerIndex,nDoF) *= qiBound_u;
            constraintsLOW.segment(innerIndex,nDoF) *= qiBound_l;
            innerIndex += nDoF;

            break;
        }
        case geo::constraint_finalConfiguration: {
            constraintsUP.segment(innerIndex,nDoF) *= qfBound_u;
            constraintsLOW.segment(innerIndex,nDoF) *= qfBound_l;
            innerIndex += nDoF;

            break;
        }
        case geo::constraint_initialGeneralizedVelocity: {
            constraintsUP.segment(innerIndex,nDoF) *= dqiBound_u;
            constraintsLOW.segment(innerIndex,nDoF) *= dqiBound_l;
            innerIndex += nDoF;

            break;
        }
        case geo::constraint_finalGeneralizedVelocity: {
            constraintsUP.segment(innerIndex,nDoF) *= dqfBound_u;
            constraintsLOW.segment(innerIndex,nDoF) *= dqfBound_l;
            innerIndex += nDoF;

            break;
        }
        case geo::constraint_pelvisSymmetry: {
            constraintsUP.segment(innerIndex,S.size()) *= pelvisBound_u;
            constraintsLOW.segment(innerIndex,S.size()) *= pelvisBound_l;
            innerIndex += S.size();

            break;
        }
        case geo::constraint_centerOfMass: {
            constraintsUP.segment(innerIndex,2*S.size()) *= comBound_u;
            constraintsLOW.segment(innerIndex,2*S.size()) *= comBound_l;
            innerIndex += 2*S.size();

            break;
        }
        case geo::constraint_centroidalMomentum: {
            constraintsUP.segment(innerIndex,6*S.size()) *= muBound_u;
            constraintsLOW.segment(innerIndex,6*S.size()) *= muBound_l;
            innerIndex += 6*S.size();

            break;
        }
        default:
            cout<<"Unrecognized constraint"<<endl;

            break;
        }
    }

    Map<geo::VectorXr>( g_u, m, 1 ) = constraintsUP;
    Map<geo::VectorXr>( g_l, m, 1 ) = constraintsLOW;

    return true;
}

//! returns the initial point for the problem
bool PracticeNLP::get_starting_point(Ipopt::Index n, 
                                     bool init_x,
                                     Ipopt::Number* x,
                                     bool init_z,
                                     Ipopt::Number* z_L,
                                     Ipopt::Number* z_U,
                                     Ipopt::Index m,
                                     bool init_lambda,
                                     Ipopt::Number* lambda)
{
    //! You can provide starting values for the dual variables if you wish

    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    //! set primal solution as linear interpolation
    geo::MatrixXr primalSolution = initialConfiguration.replicate(numberControlPoints,1) + Eigen::kroneckerProduct(geo::VectorXr::LinSpaced(numberControlPoints, 0, 1),finalConfiguration-initialConfiguration);
    Map<geo::MatrixXr>( x, n, 1 ) = primalSolution;

    return true;
}

//! returns the value of the objective function
bool PracticeNLP::eval_f(Ipopt::Index n, 
                         const Ipopt::Number* x,
                         bool new_x, Ipopt::Number& obj_value)
{

    assert(n == controlPoints.size());

    //! casting decision variable
    MapVec controlPoints( x, n );

    //! objective function evaluation, partial derivatives disabled
    robotNonlinearProblem->computeObjectiveFunction(controlPoints, false, false);

    //! retrieve cost
    obj_value = robotNonlinearProblem->getCost();

    return true;

}

//! return the gradient of the objective function grad_{x} f(x)
bool PracticeNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{

    assert(n == controlPoints.size());

    //! casting decision variable
    MapVec controlPoints( x, n );

    //! objective function evaluation, first partial derivative enabled
    robotNonlinearProblem->computeObjectiveFunction(controlPoints, true, false);

    //! casting the retrieved cost gradient
    Map<geo::MatrixXr>( grad_f, 1, n ) = robotNonlinearProblem->getCostGradient();

    return true;
}

//! return the value of the constraints: g(x)
bool PracticeNLP::eval_g(Ipopt::Index n, 
                         const Ipopt::Number* x,
                         bool new_x,
                         Ipopt::Index m,
                         Ipopt::Number* g)
{

//    assert(n = controlPoints.size());
//    assert(m = 4*nDoF);

    //! Casting decision variable
    MapVec controlPoints( x, n );

    //! constraints evaluation, partial derivatives disabled
    robotNonlinearProblem->computeConstraints(controlPoints, false, false);

    //! casting the retrieved constraints
    Map<geo::MatrixXr>( g, m, 1 ) = robotNonlinearProblem->getConstraints();

    return true;
}

//! return the structure or values of the jacobian
bool PracticeNLP::eval_jac_g(Ipopt::Index n, 
                             const Ipopt::Number* x, bool new_x,
                             Ipopt::Index m,
                             Ipopt::Index nele_jac,
                             Ipopt::Index* iRow,
                             Ipopt::Index *jCol,
                             Ipopt::Number* values)
{
    if (values == NULL) {

        Ipopt::Index idx=0;
        for (Ipopt::Index row = 0; row < m; row++) {
            for (Ipopt::Index col = 0; col < n; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

    }
    else {

        //! Casting decision variable
        MapVec controlPoints( x, n );

        //! constraints evaluation, partial derivatives enabled
        robotNonlinearProblem->computeConstraints(controlPoints, true, false);

        //! casting the retrieved constraints Jacobian
        Map<geo::MatrixXr>( values, m, n ) = robotNonlinearProblem->getConstraintsJacobian();

    }

    return true;
}

//! Hessian of the Lagrangian
bool PracticeNLP::eval_h(Ipopt::Index n,
                         const Ipopt::Number* x,
                         bool new_x,
                         Ipopt::Number obj_factor,
                         Ipopt::Index m,
                         const Ipopt::Number* lambda,
                         bool new_lambda,
                         Ipopt::Index nele_hess,
                         Ipopt::Index* iRow,
                         Ipopt::Index* jCol,
                         Ipopt::Number* values)
{
    if (values == NULL) {

        Ipopt::Index idx=0;
        Ipopt::Index j__ = 0;
        for (Ipopt::Index i_ = 0; i_ < n; i_++) {
            for (Ipopt::Index j_ = 0; j_ < (n-j__); j_++) {
                iRow[idx] = j_;
                jCol[idx] = j_ + j__;
                idx++;
            }
            j__++;
        }

        assert(idx == nele_hess);
    }
    else {
        //! Casting decision variable
        MapVec controlPoints( x, n );

        //! Casting lambda
        MapVec eigenLambda( lambda, m );

        //! Cost-function Hessian evaluation
        if (new_x) {
            if (obj_factor!=0) robotNonlinearProblem->computeObjectiveFunction(controlPoints, true, true);
            robotNonlinearProblem->computeConstraints(controlPoints, true, true);
        }

        //! Triangulating and casting the retrieved symmetric cost Hessian
        geo::MatrixXr costHessian = obj_factor*robotNonlinearProblem->getCostHessian();

        //! Triangulating and casting the retrieved symmetric constraints Hessian
        geo::MatrixXr gHessian = eigenLambda.transpose()*robotNonlinearProblem->getConstraintsHessian();
        gHessian.resize(n, n);

        //! Hessian of the Lagrangian
        geo::MatrixXr HessLagrange = costHessian + gHessian;

        geo::VectorXr stackHessLagrange(nele_hess);

        Ipopt::Index k_ = 0;
        for (int i__=0;i__<n;i__++){
            stackHessLagrange.segment(k_,n-i__) = HessLagrange.diagonal(i__);
            k_ += n-i__;
        }

        Map<geo::MatrixXr>( values, nele_hess, 1 ) = stackHessLagrange;

    }

    return true;
}

void PracticeNLP::finalize_solution(Ipopt::SolverReturn status,
                                    Ipopt::Index n,
                                    const Ipopt::Number* x,
                                    const Ipopt::Number* z_L,
                                    const Ipopt::Number* z_U,
                                    Ipopt::Index m,
                                    const Ipopt::Number* g,
                                    const Ipopt::Number* lambda,
                                    Ipopt::Number obj_value,
                                    const Ipopt::IpoptData* ip_data,
                                    Ipopt::IpoptCalculatedQuantities* ip_cq)
{

    //! Casting decision variable
    MapVec controlPoints( x, n );

    bool printSolution, saveData, BSplineInterpolation, remoteApiCoppelia;
    printSolution           = false;
    saveData                = false;
    BSplineInterpolation    = false;
    remoteApiCoppelia       = true;

    //! Print optimal result
    if(printSolution){
        cout << endl << endl << "Solution of the primal variables, x" << endl;
        cout << controlPoints.transpose() << endl;

        robotNonlinearProblem->computeObjectiveFunction(controlPoints, true, false);
        cout << endl << endl << "f(x*) = " << endl;
        cout << robotNonlinearProblem->getCost() << endl;
        //    cout << endl << endl << "Gradient of cost function" << endl;
        //    cout << robotNonlinearProblem->getCostGradient().transpose() << endl;

        robotNonlinearProblem->computeConstraints(controlPoints, true, true);
        cout << endl << endl << "CONSTRAINTS = ----------------------------------------------" << endl;
        cout << endl << "Initial config: " << endl;         short int ID = 0;
        cout << robotNonlinearProblem->getConstraints().segment(ID,nDoF).transpose() << endl << endl;
        cout << endl << "Final config: " << endl;           ID += nDoF;
        cout << robotNonlinearProblem->getConstraints().segment(ID,nDoF).transpose() << endl << endl;
        cout << endl << "Initial velocity: " << endl;       ID += nDoF;
        cout << robotNonlinearProblem->getConstraints().segment(ID,nDoF).transpose() << endl << endl;
        cout << endl << "Final velocity: " << endl;         ID += nDoF;
        cout << robotNonlinearProblem->getConstraints().segment(ID,nDoF).transpose() << endl << endl;
        cout << endl << "Pelvis symmetry: " << endl;        ID += nDoF;
        cout << robotNonlinearProblem->getConstraints().segment(ID,numberPartitions+1).transpose() << endl << endl;
        cout << endl << "Center of mass: " << endl;         ID += numberPartitions+1;
        cout << robotNonlinearProblem->getConstraints().segment(ID,2*(numberPartitions+1)).transpose() << endl << endl;
        cout << endl << "Centroidal momentum: " << endl;    ID += 2*(numberPartitions+1);
        cout << robotNonlinearProblem->getConstraints().segment(ID,6*(numberPartitions+1)).transpose() << endl << endl;

        //    cout << endl << endl << "Jacobian of constraints" << endl;
        //    cout << robotNonlinearProblem->getConstraintsJacobian().transpose() << endl;

        //    geo::MatrixXrColMajor aux = robotNonlinearProblem->getBasis()*controlPoints;
        //    aux.resize(nDoF,S.size());
        //    cout << endl << endl << "Optimal trajectory = " << endl;
        //    cout << std::scientific << std::setprecision(12) << aux.transpose() << endl;
    }

    //! Save q, dq and ddq data in files
    if(saveData){
        geo::MatrixXrColMajor qTrajectory, dqTrajectory, ddqTrajectory;
        qTrajectory = robotNonlinearProblem->getBasis()*controlPoints;
        dqTrajectory = robotNonlinearProblem->getDBasis()*controlPoints;
        ddqTrajectory = robotNonlinearProblem->getDDBasis()*controlPoints;

        qTrajectory.resize(nDoF,S.size());
        dqTrajectory.resize(nDoF,S.size());
        ddqTrajectory.resize(nDoF,S.size());

        qTrajectory.topRows(6) *= -1;
        dqTrajectory.topRows(6) *= -1;
        ddqTrajectory.topRows(6) *= -1;

//        myfile_1.open ("qData.txt");
//        myfile_2.open ("dqData.txt");
//        myfile_3.open ("ddqData.txt");

        myfile_1 << qTrajectory << endl;
        myfile_2 << dqTrajectory << endl;
        myfile_3 << ddqTrajectory << endl;
    }

    //! B-Spline extension of vector time
    if(BSplineInterpolation){
        numberPartitions = 99;
        S = geo::VectorXr::LinSpaced(numberPartitions+1, si, sf);
        robotNonlinearProblem->buildBasisFunctions(numberControlPoints, S);

        geo::MatrixXrColMajor qTrajectory, dqTrajectory, ddqTrajectory;
        qTrajectory = robotNonlinearProblem->getBasis()*controlPoints;
        dqTrajectory = robotNonlinearProblem->getDBasis()*controlPoints;
        ddqTrajectory = robotNonlinearProblem->getDDBasis()*controlPoints;

        qTrajectory.resize(nDoF,S.size());
        dqTrajectory.resize(nDoF,S.size());
        ddqTrajectory.resize(nDoF,S.size());

        qTrajectory.topRows(6) *= -1;
        dqTrajectory.topRows(6) *= -1;
        ddqTrajectory.topRows(6) *= -1;

//        myfile_1.open ("qData.txt");
//        myfile_2.open ("dqData.txt");
//        myfile_3.open ("ddqData.txt");

        myfile_1 << qTrajectory << endl;
        myfile_2 << dqTrajectory << endl;
        myfile_3 << ddqTrajectory << endl;
    }

    //! Coppelia remote streaming
    if(remoteApiCoppelia){
        //! B-Spline extension of vector time
        int prevNumberPartitions = numberPartitions;
        numberPartitions = 120;  //! 120
        S = geo::VectorXr::LinSpaced(numberPartitions+1, si, sf);
        robotNonlinearProblem->buildBasisFunctions(numberControlPoints, S);

        geo::MatrixXrColMajor qTrajectory, dqTrajectory, ddqTrajectory;
        qTrajectory = robotNonlinearProblem->getBasis()*controlPoints;
        dqTrajectory = robotNonlinearProblem->getDBasis()*controlPoints;
        ddqTrajectory = robotNonlinearProblem->getDDBasis()*controlPoints;

        qTrajectory.resize(nDoF,S.size());
        dqTrajectory.resize(nDoF,S.size());
        ddqTrajectory.resize(nDoF,S.size());

        qTrajectory.topRows(6) *= -1;
        dqTrajectory.topRows(6) *= -1;
        ddqTrajectory.topRows(6) *= -1;

        int portNb = 19997;
        int clientID = simxStart("127.0.0.1", portNb, true, true, 2500, 5);

        std::vector< const simxChar* > CharHandle{"LAnkleRoll3","LAnklePitch3","LKneePitch3","LHipPitch3","LHipRoll3","LHipYawPitch3",
                                                  "LShoulderPitch3","LShoulderRoll3","LElbowYaw3","LElbowRoll3","LWristYaw3",
                                                  "HeadYaw","HeadPitch",
                                                  "RShoulderPitch3","RShoulderRoll3","RElbowYaw3","RElbowRoll3","RWristYaw3",
                                                  "RHipYawPitch3","RHipRoll3","RHipPitch3","RKneePitch3","RAnklePitch3","RAnkleRoll3"};

        std::vector< int > IntHandle(nDoF);
        std::vector< int > FeasibleHandle(nDoF);
        int FeasibilitySum = 0;

        //! Enable handles
        for(short int ID = 0 ; ID < nDoF ; ID++){
            FeasibleHandle.at(ID) = simxGetObjectHandle(clientID, CharHandle.at(ID), &IntHandle.at(ID), simx_opmode_blocking);
            FeasibilitySum += FeasibleHandle.at(ID);
        }

        //! Send data init
        for(short int ID = 0 ; ID < nDoF ; ID++){
//            simxSetJointTargetPosition(clientID, IntHandle.at(ID), qTrajectory(ID,0),  simx_opmode_blocking);
        }

        //! Verify connection
//        cout << "client: " << clientID << endl;
//        cout << "feasibility: " << FeasibilitySum << endl;
        if (clientID == -1 || FeasibilitySum){
            cout << endl << "Coppelia connection failed" << endl;
            simxFinish(clientID);
        } else {
            cout << endl << "Connected to Coppelia remote API server" << endl;
        }

//        simxStartSimulation(clientID, simx_opmode_blocking); // start the simulation
        simxSynchronous(clientID, true); // Enable the synchronous mode

        geo::real_t inc_s = S(1) - S(0);
        int inc_s_int;  inc_s_int = (int)(inc_s*1000);

        //! Trajectory streaming
        for(int s_iter = 0; s_iter < S.size() ; s_iter++){

            //! Send data
            simxPauseCommunication(clientID, true);
            for(short int ID = 0 ; ID < nDoF ; ID++){
                simxSetJointTargetPosition(clientID, IntHandle.at(ID), qTrajectory(ID,s_iter), simx_opmode_oneshot);
            }
            simxPauseCommunication(clientID, false);

            simxSynchronousTrigger(clientID);

//            std::cout<<"time = "<<inc_s_int<<std::endl;
//            extApi_sleepMs(inc_s_int);  //! robotized movement

        }
        simxFinish(clientID);
        cout << "Streaming finished" << endl;

        //! Restore basis
        numberPartitions = prevNumberPartitions;
        S = geo::VectorXr::LinSpaced(numberPartitions+1, si, sf);
        robotNonlinearProblem->buildBasisFunctions(numberControlPoints, S);
    }

}

