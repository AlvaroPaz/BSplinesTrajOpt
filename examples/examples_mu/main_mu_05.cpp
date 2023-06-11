// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file examples_com/main_com_06.cpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Optimal motion generation for inertial Nao (sequential movements) -> Sky
  *      -> Before executing this example you must run the com_example_01 first
 */

#include "geombd/core.h"
#include "geombd/dynamics.h"
#include "geombd/trajectoryOptimization.h"

#include <Eigen/Dense>

#include <memory>
#include <fstream>
#include <sstream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//std::string naoFile = "../../../data/TROmodels/lwr.urdf";
std::string naoFile = "../../../data/TROmodels/atlas_cpp_reduced.urdf";

#include <IpIpoptApplication.hpp>
#include <iostream>
#include <iomanip>

#include "../ipopt/ipopt_interface_atlas_iner_com.hpp"


// using namespace Ipopt;


using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000


int main(int argv, char* argc[])
{

    //! Build multibody
    geo::World world = geo::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, geo::kNAO);


    std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);

    geo::VectorXr q = robot->getConfiguration();

    q.setZero();

    robot->setConfiguration(q);


    //! Optimization parameters
    std::shared_ptr< geo::robotSettingsTrajectoryOptimization > optSettings(new geo::robotSettingsTrajectoryOptimization);

    optSettings->n = robot->getDoF();
    optSettings->numberControlPoints = 4;//4
    optSettings->numberPartitions    = 7;//7  15
    optSettings->si = 0.0;//0.0
    optSettings->sf = 3.5;//0.1 or 0.2 or 25.0  5.0
    optSettings->S = geo::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = geo::wrt_controlPoints;
    geo::VectorXr weights;  weights.setOnes(optSettings->n);
    weights.segment(0,9) *= 1e-5;  // si minimiza el torque en soporte aumenta el equilibrio
    weights.tail(6) *= 1e-4;
    optSettings->weights = weights;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );

    //! Set derivation routine
    optSettings->deriveRoutine = geo::_Analytic_;
//    optSettings->deriveRoutine = geo::_Numeric_;
//    optSettings->deriveRoutine = geo::_BFGS_;


    geo::VectorXr mu_aux;
    mu_aux = geo::VectorXr::Ones(6,1);
    mu_aux << 5, 5, 28, 8, 8, 8;

    //! Boundaries
    geo::VectorXr bound_aux;
    bound_aux = geo::VectorXr::Ones(optSettings->n,1) * 1e-3;

    bound_aux.segment(9,7) = geo::VectorXr::Ones(7,1) * 0.6;  // arms 0.4
    bound_aux.segment(17,7) = geo::VectorXr::Ones(7,1) * 0.6;

    bound_aux.segment(6,3) = geo::VectorXr::Ones(3,1) * 1.2;

//    bound_aux.segment(6,18) = geo::VectorXr::Ones(18,1) * 2;  // 3

    optSettings->qiBound_u = geo::VectorXr::Ones(optSettings->n,1) * 1e-3;
    optSettings->qfBound_u = bound_aux;
    optSettings->comBound_u = geo::VectorXr::Ones(2,1) * 0.035;  // 0.045
    optSettings->muBound_u = mu_aux;  // 1e-2

    optSettings->qiBound_l = - geo::VectorXr::Ones(optSettings->n,1) * 1e-3;
    optSettings->qfBound_l = -bound_aux;
    optSettings->comBound_l = -geo::VectorXr::Ones(2,1) * 0.035;
    optSettings->muBound_l = -mu_aux;


    //! Generalized Configurations
    geo::VectorXr q_1(optSettings->n,1), q_2(optSettings->n,1), q_3(optSettings->n,1), q_4(optSettings->n,1), q_5(optSettings->n,1);



    q_4 << -0.00100001, 0.901, -2.099, 1.501, -0.00100001, -0.00100001, -0.663225, -0.219388, -0.101, 0.6, -0.00638687, 0.6, 0.00166667, 0.6, -0.462452, -0.599996, 0.000999941, -0.6, -0.6, 0, -0.6, -0.6, -0.6, -0.6, -0.00100001, 0.00100001, 0.65764, 1.499, 0.000999824, 0.00100001;


//    q_5 << 0, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.3, 0.6, -0.3, 0;

//    q_5 << 0, 0.9, -2.1, 1.5, 0, 0, 0, 0, -0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.657, 1.5, 0, 0;

    q_5 << 0, 0.1, -0.2, 0.1, -0.52, 0, 0, 0, -0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.7, -0.52, -1.61, 0, -1, 0;



    std::vector< geo::VectorXr > Q_input;
    Q_input.clear();
    Q_input.push_back(q_4);
    Q_input.push_back(q_5);

    //! Fill up stack of constraints
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
    optSettings->initialConfiguration = q_4;

    optSettings->StackConstraints.push_back(geo::constraint_finalConfiguration);
    optSettings->finalConfiguration = q_5;

    optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_centerOfMass);

    optSettings->StackConstraints.push_back(geo::constraint_centroidalMomentum);


    //! Create a new instance of your nlp (use a SmartPtr, not raw)
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new PracticeNLP(robot, optSettings);


    //! Create a new instance of IpoptApplication (use a SmartPtr, not raw)
    //! We are using the factory, since this allows us to compile this example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);

    //! Change some options
    app->Options()->SetNumericValue("tol", 1e-4);  // 1e-4
    app->Options()->SetIntegerValue("max_iter", 5000);
    app->Options()->SetStringValue("mu_strategy", "adaptive"); ///monotone adaptive

//    app->Options()->SetIntegerValue("file_print_level", 12);

//    app->Options()->SetStringValue("derivative_test", "second-order");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");

    //! Setting derivation routine
    if(optSettings->deriveRoutine == geo::_BFGS_) {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
      } else {
        app->Options()->SetStringValue("hessian_approximation", "exact");
      }


    //! Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }


    //! Ask Ipopt to solve the problem
    for(short int innerIter; innerIter < Q_input.size()-1; innerIter++ ){
        optSettings->initialConfiguration = Q_input.at(innerIter);
        optSettings->finalConfiguration = Q_input.at(innerIter+1);
        status = app->OptimizeTNLP(mynlp);

        if (status == Ipopt::Solve_Succeeded) {
            std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
        }
        else {
            std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
        }

    }


    //! As the SmartPtrs go out of scope, the reference count will be decremented and the objects will automatically be deleted.
    return (int) status;

//    return 0;

}
