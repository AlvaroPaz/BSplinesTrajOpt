// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file ipopt/ipopt_interface_nao_iner_com.cpp
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

#include "ipopt_interface_atlas_iner_com.hpp"

typedef Eigen::Map<const geo::VectorXr> MapVec;

std::ofstream myfile_1, myfile_2, myfile_3, myfile_4, myfile_5;

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

    // Retrieve the weights vector
    this->weights = robotSettings->weights;

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

    itHess = 0;

    myfile_1.open ("../../../../../Dropbox/Atlas_project/matlab_data/qData.txt");
    myfile_2.open ("../../../../../Dropbox/Atlas_project/matlab_data/dqData.txt");
    myfile_3.open ("../../../../../Dropbox/Atlas_project/matlab_data/ddqData.txt");
    myfile_4.open ("../../../../../Dropbox/Atlas_project/matlab_data/com.txt");
    myfile_5.open ("../../../../../Dropbox/Atlas_project/matlab_data/mu.txt");

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
    geo::real_t dqiBound_u = 1e-3;            //! Velocity 1e-3
    geo::real_t dqfBound_u = 1e-3;
    geo::real_t pelvisBound_u = 1e-3;         //! Pelvis symmetry 1e-10

    geo::real_t dqiBound_l = -1e-3;           //! Velocity -1e-3
    geo::real_t dqfBound_l = -1e-3;
    geo::real_t pelvisBound_l = -1e-3;        //! Pelvis symmetry -1e-10


    //! Fill up constraints
    int innerIndex = 0;
    for(geo::ConstraintsStack::iterator it = StackConstraints.begin(); it != StackConstraints.end(); ++it) {
        switch ( *it ) {
        case geo::constraint_initialConfiguration: { // 1 and 1.5 well
            constraintsUP.segment(innerIndex,nDoF) = robotSettings->qiBound_u;
            constraintsLOW.segment(innerIndex,nDoF) = robotSettings->qiBound_l;
            innerIndex += nDoF;

            break;
        }
        case geo::constraint_finalConfiguration: {
            constraintsUP.segment(innerIndex,nDoF) = robotSettings->qfBound_u;
            constraintsLOW.segment(innerIndex,nDoF) = robotSettings->qfBound_l;
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
            constraintsUP.segment(innerIndex,2*S.size()) = robotSettings->comBound_u.replicate( S.size(), 1 );
            constraintsLOW.segment(innerIndex,2*S.size()) = robotSettings->comBound_l.replicate( S.size(), 1 );
            innerIndex += 2*S.size();

            std::cout<<"===> Center of Mass Criteria Enabled"<<std::endl;

            break;
        }
        case geo::constraint_centroidalMomentum: {
            constraintsUP.segment(innerIndex,6*S.size()) = robotSettings->muBound_u.replicate( S.size(), 1 );
            constraintsLOW.segment(innerIndex,6*S.size()) = robotSettings->muBound_l.replicate( S.size(), 1 );
            innerIndex += 6*S.size();

            std::cout<<"===> Centroidal Momentum Criteria Enabled"<<std::endl;

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
    robotNonlinearProblem->computeObjectiveFunction(controlPoints, weights, false, false);

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

  if(robotSettings->deriveRoutine == geo::_Analytic_) {

      //! objective function evaluation, first partial derivative enabled
      robotNonlinearProblem->computeObjectiveFunction(controlPoints, weights, true, false);

      //! casting the retrieved cost gradient
      Map<geo::MatrixXr>( grad_f, 1, n ) = robotNonlinearProblem->getCostGradient();

    } else {

      //! first objective function evaluation
      robotNonlinearProblem->computeObjectiveFunction(controlPoints, weights, false, false);

      auto cost_exa = robotNonlinearProblem->getCost();
      numericD_cost = geo::MatrixXr::Zero(1, n);

      double inc_c = 1e-8;
      for(short int iter = 0 ; iter < n ; iter++){
          //! wrt c
          c_aux = controlPoints;
          c_aux(iter) += inc_c;

          //! objective function evaluation
          robotNonlinearProblem->computeObjectiveFunction(c_aux, weights, false, false);

          auto temporalCost = robotNonlinearProblem->getCost();

          numericD_cost(iter) = (temporalCost-cost_exa)/inc_c;
        }

      //! casting the retrieved cost gradient
      Map<geo::MatrixXr>( grad_f, 1, n ) = numericD_cost;
    }

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
    robotNonlinearProblem->computeConstraintsII(controlPoints, false, false);

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

      if(robotSettings->deriveRoutine == geo::_Analytic_) {

          //! constraints evaluation, first partial derivative enabled
          robotNonlinearProblem->computeConstraintsII(controlPoints, true, false);

          //! casting the retrieved constraints Jacobian
          Map<geo::MatrixXr>( values, m, n ) = robotNonlinearProblem->getConstraintsJacobian();

        } else {
          //! first constraints evaluation
          robotNonlinearProblem->computeConstraints(controlPoints, false, false);

          auto Const_exa = robotNonlinearProblem->getConstraints();
          numericD_rest = geo::MatrixXr::Zero(m, n);

          double inc_c = 1e-8;
          for(short int iter = 0 ; iter < n ; iter++){
              //! wrt c
              c_aux = controlPoints;
              c_aux(iter) += inc_c;

              //! constraints evaluation
              robotNonlinearProblem->computeConstraints(c_aux, false, false);

              auto temporalConst = robotNonlinearProblem->getConstraints();

              numericD_rest.col(iter) = (temporalConst-Const_exa)/inc_c;
            }

          //! casting the retrieved constraints Jacobian
          Map<geo::MatrixXr>( values, m, n ) = numericD_rest;
        }

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

      if(robotSettings->deriveRoutine == geo::_Analytic_) {

          Ipopt::Index idx=0;
          for (Ipopt::Index i_ = 0; i_ < n; i_++) {
              for (Ipopt::Index j_ = i_; j_ < n; j_++) {
                  iRow[idx] = i_;
                  jCol[idx] = j_;
                  idx++;
                }
            }
          assert(idx == nele_hess);

        } else {

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
    }
  else {

      if (itHess == 0) {

          stackHessLagrange = geo::VectorXr::Zero(nele_hess,1);

          //! Casting decision variable
          MapVec controlPoints( x, n );

          //! Casting lambda
          MapVec eigenLambda( lambda, m );

          if(robotSettings->deriveRoutine == geo::_Analytic_) {


              //! Cost-function Hessian evaluation
              if (obj_factor!=0) robotNonlinearProblem->computeContractedCostFunction(controlPoints, weights, true, true);
              robotNonlinearProblem->computeContractedConstraints(controlPoints, true, true);

              //! Triangulating and casting the retrieved symmetric cost Hessian
              stackHessLagrange = obj_factor*robotNonlinearProblem->getContractedCostHessian();

              //! Triangulating and casting the retrieved symmetric constraints Hessian
              stackHessLagrange += robotNonlinearProblem->getContractedConstraintsHessian().transpose() * eigenLambda;


            } else {
              //! first functions evaluation
              robotNonlinearProblem->computeObjectiveFunction(controlPoints, weights, false, false);
              robotNonlinearProblem->computeConstraints(controlPoints, false, false);

              auto cost_exa = robotNonlinearProblem->getCost();
              auto rest_exa = robotNonlinearProblem->getConstraints();

              //! f(x+h) and f(x+k) evaluation
              double h__ = 1e-4;
              double k__ = 1e-4;
              double hk__ = h__*k__;

              geo::VectorXr f_h(n), f_k(n);
              geo::MatrixXr g_h(m,n), g_k(m,n);

              for(short int iter = 0 ; iter < n ; iter++){
                  //! h variation
                  c_aux = controlPoints;
                  c_aux(iter) += h__;

                  //! functions evaluation
                  robotNonlinearProblem->computeObjectiveFunction(c_aux, weights, false, false);
                  f_h(iter) = robotNonlinearProblem->getCost();

                  robotNonlinearProblem->computeConstraintsII(c_aux, false, false);
                  g_h.col(iter) = robotNonlinearProblem->getConstraints();

                  //! k variation
                  c_aux = controlPoints;
                  c_aux(iter) += k__;

                  //! functions evaluation
                  robotNonlinearProblem->computeObjectiveFunction(c_aux, weights, false, false);
                  f_k(iter) = robotNonlinearProblem->getCost();

                  robotNonlinearProblem->computeConstraints(c_aux, false, false);
                  g_k.col(iter) = robotNonlinearProblem->getConstraints();
                }

              geo::MatrixXr temporalD_cost(1, n);
              geo::MatrixXr temporalD_rest(m, n);

              numericDD_cost = geo::MatrixXr::Zero(1, n*n);
              numericDD_rest = geo::MatrixXr::Zero(m, n*n);

              for(int i__ = 0 ; i__ < n ; i__++){
                  //! h variation
                  c_aux = controlPoints;
                  c_aux(i__) += h__;

                  for(int j__ = 0 ; j__ < n ; j__++){
                      //! k variation
                      c_aux_in = c_aux;
                      c_aux_in(j__) += k__;

                      //! f(x+h+k) evaluation
                      robotNonlinearProblem->computeObjectiveFunction(c_aux_in, weights, false, false);

                      double f_h_k = robotNonlinearProblem->getCost();

                      temporalD_cost(j__) = (f_h_k - f_h(i__) - f_k(j__) + cost_exa)/hk__;

                      //! g(x+h+k) evaluation
                      robotNonlinearProblem->computeConstraints(c_aux_in, false, false);

                      geo::VectorXr g_h_k = robotNonlinearProblem->getConstraints();

                      temporalD_rest.col(j__) = (g_h_k - g_h.col(i__) - g_k.col(j__) + rest_exa)/hk__;
                    }
                  numericDD_cost.middleCols(i__*n, n) = temporalD_cost;
                  numericDD_rest.middleCols(i__*n, n) = temporalD_rest;
                }

              //! Triangulating and casting the retrieved symmetric cost Hessian
              costHessian = obj_factor*numericDD_cost;
              costHessian.resize(n, n);

              //! Triangulating and casting the retrieved symmetric constraints Hessian
              gHessian = eigenLambda.transpose()*numericDD_rest;
              gHessian.resize(n, n);

              //! Hessian of the Lagrangian
              HessLagrange = costHessian + gHessian;

              Ipopt::Index k_ = 0;
              for (int i__=0;i__<n;i__++){
                  stackHessLagrange.segment(k_,n-i__) = HessLagrange.diagonal(i__);
                  k_ += n-i__;
                }
            }
        }

      Map<geo::MatrixXr>( values, nele_hess, 1 ) = stackHessLagrange;

      itHess++;
      if(itHess == 18) itHess = 0;  //! default 18 for com
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
    BSplineInterpolation    = true;
    remoteApiCoppelia       = false;

    //! Print optimal result
    if(printSolution){
        cout << endl << endl << "Solution of the primal variables, x" << endl;
        cout << controlPoints.transpose() << endl;

        robotNonlinearProblem->computeObjectiveFunction(controlPoints, weights, true, false);
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
        numberPartitions = 100;  // 20 per second, assuming inc_t = 0.05
        S = geo::VectorXr::LinSpaced(numberPartitions+1, si, sf);
        robotNonlinearProblem->buildBasisFunctions(numberControlPoints, S);

        geo::MatrixXrColMajor qTrajectory, dqTrajectory, ddqTrajectory;
        qTrajectory = robotNonlinearProblem->getBasis()*controlPoints;
        dqTrajectory = robotNonlinearProblem->getDBasis()*controlPoints;
        ddqTrajectory = robotNonlinearProblem->getDDBasis()*controlPoints;

        qTrajectory.resize(nDoF,S.size());
        dqTrajectory.resize(nDoF,S.size());
        ddqTrajectory.resize(nDoF,S.size());

        //! Save extra data
        //! --------------------------------------------------------------------------------
        auto robotDynamics = std::make_shared< geo::InverseDynamics >( robot );
        geo::VectorXr q__(nDoF,1), dq__(nDoF,1), ddq__(nDoF,1);

        for(int k__ = 0 ; k__ < S.size() ; k__++){
            q__ = qTrajectory.col(k__);  dq__ = dqTrajectory.col(k__);  ddq__ = ddqTrajectory.col(k__);

            robotDynamics->setGeneralizedCoordinates(q__, dq__, ddq__);
            robotDynamics->computeContractedCenterOfMass(false, false);
            robotDynamics->computeCentroidalMomentum(false, false);

            myfile_4 << robotDynamics->getRobotCoM().transpose() << endl;
            myfile_5 << robotDynamics->getCentroidalMomentum().transpose() << endl;
          }

        //! --------------------------------------------------------------------------------

        qTrajectory.topRows(6) *= -1;
        dqTrajectory.topRows(6) *= -1;
        ddqTrajectory.topRows(6) *= -1;

        myfile_1 << qTrajectory << endl;
        myfile_2 << dqTrajectory << endl;
        myfile_3 << ddqTrajectory << endl;

    }

    //! Coppelia remote streaming
    if(remoteApiCoppelia){
        //! B-Spline extension of vector time
        int prevNumberPartitions = numberPartitions;
        numberPartitions = 220;  //! 220
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

