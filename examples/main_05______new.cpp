// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file examples/main_02.cpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Optimal motion generation for Nao (sequential movements) - 24DoF
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

std::string naoFile = "../../BSplinesTrajOpt/data/nao_python.urdf";


#include <IpIpoptApplication.hpp>
#include <iostream>

#include "../ipopt/ipopt_interface_05.hpp"


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
    optSettings->numberPartitions    = 25;//7 12
    optSettings->si = 0.0;//0.0
    optSettings->sf = 4.0;//0.1 or 0.2 or 25.0
    optSettings->S = geo::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = geo::wrt_controlPoints;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );


    geo::VectorXr q_1(optSettings->n,1), q_2(optSettings->n,1), q_3(optSettings->n,1), q_4(optSettings->n,1), q_5(optSettings->n,1), q_6(optSettings->n,1), q_7(optSettings->n,1), q_8(optSettings->n,1), q_9(optSettings->n,1);



    q_1 << 0, 0, 0, 0, 0, 0,
           0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0.035, 0, -0.035, 0,
           0, -0.035, 0, 0.035, 0;

    q_2 << 0, 0, 0, 0, 0, 0,
           0, 0,
           0.2, 0, -1.535, 2.112, -1.189, 0,
           0.2, 0, -1.535, 2.12, -1.186, 0,
           0, 0.035, 0, -0.035, 0,
           0, -0.035, 0, 0.035, 0;

    q_3 << 0, 0, 0, 0, 0, 0,
           0, 0,
           -0.3, 0, -1.535, 2.112, -1.189, 0,
           -0.3, 0, -1.535, 2.12, -1.186, 0,
           0, 0.035, 0, -0.035, 0,
           0, -0.035, 0, 0.035, 0;

    q_4 << 0, 0, 0, 0, 0, 0,
           0, 0,
           -0.3, 0, -1.535, 2.112, -1.189, 0,
           -0.3, 0, -1.535, 2.12, -1.186, 0,
           -0.11, 0.035, -1.56, -1.56, 0,
           -0.11, -0.035, 1.56, 1.56, 0;

    q_5 << 0, 0, 0, 0, 0, 0,
           0, -0.671,
           0, 0, -1.535, -0.092, 0.922, 0,
           0, 0, 0.484, 0, 0, 0,
           -0.8, 0.035, -1.56, -1.56, 0,
           -0.8, -0.035, 1.56, 1.56, 0;

    q_6 << 0, 0, 0, 0, 0, 0,
           0, -0.671,
           0, 0, -1.535, -0.092, 0.922, 0,
           0, 0, 0.484, 2.1201, 0.9320, 0,
           -0.8, 0.035, -1.56, -1.56, 0,
           -0.8, -0.035, 1.56, 1.56, 0;

    q_7 << 0, 0, 0, 0, 0, 0,
           0, -0.671,
           -1.145/2, 0.79, -1.535, -0.092, 0.922, 0,
           0.740/2, -0.79, 0.484, 2.1201, 0.9320, 0,
           -0.8, 0.035, -1.56, -1.56, 0,
           -0.8, -0.035, 1.56, 1.56, 0;

    q_8 << 0, 0, 0, 0, 0, 0,
           0, -0.671,
           -1.145/2, 0.79, -1.535, -0.092, 0.922, 0,
           0.740/2, -0.79, 0.484, 2.1201, 0.9320, 0,
           -0.8+0.4*0, 0.035+1.0, -1.56, -1.56, 0,
           -0.8+0.4*0, -0.035-1.4, 1.56, 1.56, 0;





    std::vector< geo::VectorXr > Q_input;
    Q_input.clear();
    Q_input.push_back(q_1);
    Q_input.push_back(q_2);
    Q_input.push_back(q_3);
    Q_input.push_back(q_4);
    Q_input.push_back(q_5);
    Q_input.push_back(q_6);
    Q_input.push_back(q_7);
    Q_input.push_back(q_8);
//    Q_input.push_back(q_9);

    //! Fill up stack of constraints
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
    optSettings->initialConfiguration = q_1;

    optSettings->StackConstraints.push_back(geo::constraint_finalConfiguration);
    optSettings->finalConfiguration = q_2;

    optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

//    optSettings->StackConstraints.push_back(geo::constraint_pelvisSymmetry);


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
    //    app->Options()->SetStringValue("output_file", "ipopt.out");

        //! limited-memory for BFGS and exact for our analytic
        app->Options()->SetStringValue("hessian_approximation", "exact");

        //! finite-difference-values for numeric aprox. and exact for our analytic
        app->Options()->SetStringValue("jacobian_approximation", "exact");

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




//q_1 << 0, 0, 0, 0, 0, 0,
//       0, 0,
//       0, 0, 0, 0, 0, 0,
//       0, 0, 0, 0, 0, 0,
//       0, 0.035, 0, -0.035, 0,
//       0, -0.035, 0, 0.035, 0;

//q_2 << 0, 0, 0, 0, 0, 0,
//       0, 0,
//       0.2, 0, -1.535, 2.112, -1.189, 0,
//       0.2, 0, -1.535, 2.12, -1.186, 0,
//       0, 0.035, 0, -0.035, 0,
//       0, -0.035, 0, 0.035, 0;

//q_3 << 0, 0, 0, 0, 0, 0,
//       0, 0,
//       -0.3, 0, -1.535, 2.112, -1.189, 0,
//       -0.3, 0, -1.535, 2.12, -1.186, 0,
//       0, 0.035, 0, -0.035, 0,
//       0, -0.035, 0, 0.035, 0;

//q_4 << 0, 0, 0, 0, 0, 0,
//       0, 0,
//       -0.3, 0, -1.535, 2.112, -1.189, 0,
//       -0.3, 0, -1.535, 2.12, -1.186, 0,
//       -0.11, 0.035, -1.56, -1.56, 0,
//       -0.11, -0.035, 1.56, 1.56, 0;

//q_5 << 0, 0, 0, 0, 0, 0,
//       0, -0.671,
//       0, 0, -1.535, -0.092, 0.922, 0,
//       0, 0, 0.484, 0, 0, 0,
//       -0.8, 0.035, -1.56, -1.56, 0,
//       -0.8, -0.035, 1.56, 1.56, 0;

//q_6 << 0, 0, 0, 0, 0, 0,
//       0, -0.671,
//       0, 0, -1.535, -0.092, 0.922, 0,
//       0, 0, 0.484, 2.1201, 0.9320, 0,
//       -0.8, 0.035, -1.56, -1.56, 0,
//       -0.8, -0.035, 1.56, 1.56, 0;

//q_7 << 0, 0, 0, 0, 0, 0,
//       0, -0.671,
//       -1.145/2, 0.79, -1.535, -0.092, 0.922, 0,
//       0.740/2, -0.79, 0.484, 2.1201, 0.9320, 0,
//       -0.8, 0.035, -1.56, -1.56, 0,
//       -0.8, -0.035, 1.56, 1.56, 0;

//q_8 << 0, 0, 0, 0, 0, 0,
//       0, -0.671,
//       -1.145/2, 0.79, -1.535, -0.092, 0.922, 0,
//       0.740/2, -0.79, 0.484, 2.1201, 0.9320, 0,
//       -0.8+0.4*0, 0.035+1.0, -1.56, -1.56, 0,
//       -0.8+0.4*0, -0.035-1.4, 1.56, 1.56, 0;
