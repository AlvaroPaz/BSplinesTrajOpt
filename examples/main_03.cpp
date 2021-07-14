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
 *	Main file to include the Trajectory Optimization Problem into the Ipopt Solver
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

std::string naoFile = "../../BSplinesTrajOpt/data/NAOURDF_inertial.xml";


#include <IpIpoptApplication.hpp>
#include <iostream>

#include "ipopt_interface.hpp"


// using namespace Ipopt;


using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000


int main(int argv, char* argc[])
{

    //! Build multibody
    hr::core::World world = hr::core::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, hr::core::kNAO);


    std::shared_ptr<hr::core::MultiBody> robot = world.getRobot(0);

    hr::VectorXr q = robot->getConfiguration();

    q.setZero();

    robot->setConfiguration(q);


    //! Optimization parameters
    std::shared_ptr< hr::core::robotSettingsTrajectoryOptimization > optSettings(new hr::core::robotSettingsTrajectoryOptimization);

    optSettings->n = robot->getDoF();
    optSettings->numberControlPoints = 4;//4
    optSettings->numberPartitions    = 25;//7
    optSettings->si = 0.0;//0.0
    optSettings->sf = 1.0;//0.1 or 0.2 or 25.0
    optSettings->S = hr::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = hr::core::wrt_controlPoints;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );


    hr::VectorXr q_1(optSettings->n,1), q_2(optSettings->n,1), q_3(optSettings->n,1), q_4(optSettings->n,1), q_5(optSettings->n,1), q_6(optSettings->n,1), q_7(optSettings->n,1), q_8(optSettings->n,1), q_9(optSettings->n,1);
    q_1 << 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, 0, 0, 0, 0, 0;

    q_2 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, -0.379, 0, 0, 0, 0.379;

    q_3 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 1.2000, 0.5148, 0, 0, 0, 0.035, 0, 0, -0.79, 0, 0, 0, 0.379;

    q_4 << 0.25, 0.9200, -2.0000, 1.0800, -0.3794, 0, 0, 1.3264, 0, -0.0350, 0, 1.2000, 0.5148, 0, 0.3141, 0, 1.5445, 0, 0, -0.7904, 0, 0, 0, 0.397;

    q_5 << 0.25, 0.9200, -2.0000, 1.0800, -0.3794, 0, 0, -0.3141, 0, -1.5445, 0, -1.2000, 0.5148, 0, -1.3264, 0, 0.0350, 0, 0, -0.7904, 0, 0, 0, 0.397;

    q_6 << 0.25, 0, 0, 0, -0.3794, 0, 0, -0.3141, 0, -1.5445, 0, -1.2000, 0.5148, 0, -1.3264, 0, 0.0350, 0, 0, -0.7904, -1.0800, 2.0000, -0.9200, 0.2;

    q_7 << 0.25, -1, 0, 0, -0.3794, 0, 0, -0.3141, 0, -1.5445, 0, -1.2000, 0.5148, 0, -1.3264, 0, 0.0350, 0, 0, -0.7904, -1.0800, 2.0000, -0.9200, 0.0;

    q_8 << 0.379, 0, 0, 0, -0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, 0.379, 0, 0, 0, -0.379;

    q_9 << 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, 0, 0, 0, 0, 0;

    std::vector< hr::VectorXr > Q_input;
    Q_input.clear();
    Q_input.push_back(q_1);
    Q_input.push_back(q_2);
    Q_input.push_back(q_3);
    Q_input.push_back(q_4);
    Q_input.push_back(q_5);
    Q_input.push_back(q_6);
    Q_input.push_back(q_7);
    Q_input.push_back(q_8);
    Q_input.push_back(q_9);

    //! Fill up stack of constraints
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(hr::core::constraint_initialConfiguration);
    optSettings->initialConfiguration = q_1;

    optSettings->StackConstraints.push_back(hr::core::constraint_finalConfiguration);
    optSettings->finalConfiguration = q_2;

    optSettings->StackConstraints.push_back(hr::core::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_pelvisSymmetry);


    //! Create a new instance of your nlp (use a SmartPtr, not raw)
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new PracticeNLP(robot, optSettings);


    //! Create a new instance of IpoptApplication (use a SmartPtr, not raw)
    //! We are using the factory, since this allows us to compile this example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);

    //! Change some options
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetIntegerValue("max_iter", 5000);
    //        app->Options()->SetNumericValue("max_cpu_time", 200);
    //        app->Options()->SetIntegerValue("acceptable_iter", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive"); ///monotone adaptive
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory"); /// exact
//    app->Options()->SetStringValue("derivative_test", "second-order");
    //        app->Options()->SetStringValue("jac_c_constant", "yes");
    //        app->Options()->SetStringValue("linear_solver", "mumps");
    //        app->Options()->SetStringValue("accept_every_trial_step", "yes");
    //        app->Options()->SetNumericValue("bound_relax_factor", 0.2);
    //        app->Options()->SetStringValue("corrector_type", "primal-dual");
    //        app->Options()->SetNumericValue("obj_scaling_factor", 0.5);
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    //        app->Options()->SetStringValue("alpha_for_y", "full");



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

}
