// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file examples/main_01.cpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Optimal motion generation for inertial Nao -> Airplane-like movement
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

std::string naoFile = "../../BSplinesTrajOpt/data/nao_inertial_python.urdf";

#include <IpIpoptApplication.hpp>
#include <iostream>

#include "../ipopt/ipopt_interface.hpp"


// using namespace Ipopt;


using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000


int main(int argv, char* argc[])
{
    //! Build multibody
    //! ---------------------------------------------------------------------------------------------------------- !//
    geo::World world = geo::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, geo::kNAO);

    std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);

    geo::VectorXr q = robot->getConfiguration();
    q.setZero();
    robot->setConfiguration(q);


    //! Optimization parameters
    //! ---------------------------------------------------------------------------------------------------------- !//
    std::shared_ptr< geo::robotSettingsTrajectoryOptimization > optSettings(new geo::robotSettingsTrajectoryOptimization);

    optSettings->n = robot->getDoF();
    optSettings->numberControlPoints = 4;//4
    optSettings->numberPartitions    = 7;//7  30
    optSettings->si = 0.0;//0.0
    optSettings->sf = 25.0;//0.1 or 0.2 or 25.0
    optSettings->S = geo::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = geo::wrt_controlPoints;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );
    geo::VectorXr weights;  weights.setOnes(optSettings->n);
    optSettings->weights = weights;


    geo::VectorXr q_1(optSettings->n,1), q_2(optSettings->n,1), q_3(optSettings->n,1), q_4(optSettings->n,1), q_5(optSettings->n,1);
    q_1 << 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, 0, 0, 0, 0, 0;

    q_2 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, -0.379, 0, 0, 0, 0.379;

    q_3 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, -0.79, 0, 0, 0, 0.379;

    q_4 << -0.3000, 0, 0, 0, 0, 0, 1.4112, 0.2730, -1.3730, -0.9863, -0.0062, 0.0015, 0.0214, 1.3945, -0.2731, 1.3698, 0.9879, -0.0077, 0, 0.0016, -0.4510, 0.6995, -0.3528, 0;

    q_5 << -0.2000, 0.3513, -1.1000, 1.5300, 0, 0, 1.4112, 1.2000, -1.3730, -0.9863, -0.0062, 0.0015, -0.6000, 1.3945, -1.2000, 1.3698, 0.9879, -0.0077, 0, 0.0016, 0.4800, 0.6995, -0.3528, 0;

    std::vector< geo::VectorXr > Q_input;
    Q_input.clear();
    Q_input.push_back(q_1);
    Q_input.push_back(q_2);
    Q_input.push_back(q_3);
    Q_input.push_back(q_4);
    Q_input.push_back(q_5);


    //! Fill up stack of constraints
    //! ---------------------------------------------------------------------------------------------------------- !//
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
    optSettings->initialConfiguration = q_1;

    optSettings->StackConstraints.push_back(geo::constraint_finalConfiguration);
    optSettings->finalConfiguration = q_2;

    optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_pelvisSymmetry);

    optSettings->StackConstraints.push_back(geo::constraint_centerOfMass);

//    optSettings->StackConstraints.push_back(geo::constraint_centroidalMomentum);


    //! Create a new instance of your nlp (use a SmartPtr, not raw)
    //! ---------------------------------------------------------------------------------------------------------- !//
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new PracticeNLP(robot, optSettings);

    //! Create a new instance of IpoptApplication (use a SmartPtr, not raw)
    //! We are using the factory, since this allows us to compile this example with an Ipopt Windows DLL
    //! ---------------------------------------------------------------------------------------------------------- !//
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);

    //! Set optimization options
    //! ---------------------------------------------------------------------------------------------------------- !//
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetIntegerValue("max_iter", 5000);

    app->Options()->SetStringValue("mu_strategy", "adaptive"); ///monotone adaptive

    //! limited-memory for BFGS and exact for our analytic
    app->Options()->SetStringValue("hessian_approximation", "exact");
//    app->Options()->SetStringValue("derivative_test", "second-order");
//    app->Options()->SetStringValue("jac_c_constant", "yes");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");


    //! Initialize the IpoptApplication and process the options
    //! ---------------------------------------------------------------------------------------------------------- !//
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

}
