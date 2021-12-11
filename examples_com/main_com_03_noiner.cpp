// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file examples/main_03.cpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Optimal motion generation for Nao robot (30DoF)
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

#include "../ipopt/ipopt_interface_nao_com_test.hpp"


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
    optSettings->sf = 10.0;//0.1 or 0.2 or 25.0
    optSettings->S = geo::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = geo::wrt_controlPoints;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );


    geo::VectorXr q_1(optSettings->n,1), q_2(optSettings->n,1), q_3(optSettings->n,1), q_4(optSettings->n,1), q_5(optSettings->n,1);
    q_1 << 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, 0, 0, 0, 0, 0;
    q_2 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, -0.379, 0, 0, 0, 0.379;
    q_3 << -0.379, 0, 0, 0, 0.379, 0, 0, 0, 0, -0.035, 0, 0, 0, 0, 0, 0, 0.035, 0, 0, -0.79, 0, 0, 0, 0.379;


    q_1 << 0, 0, 0, 0, 0, 0,
           0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0,
           0, 0.035, 0, -0.035, 0,
           0, -0.035, 0, 0.035, 0;

    q_2 << 0, 0, 0, 0, 0, 0,
           0, 0,
           0, 0, -1.535, 2.112, -1.189, 0,
           0, 0, -1.535, 2.12, -1.186, 0,
           0, 0.035, 0, -0.035, 0,
           0, -0.035, 0, 0.035, 0;


    //! Airplane like pose

    q_4 << 0,0,0,0,0,0,
           0.0015, 0.0214,
           0, 0, 0, 0, 0, 0.3,
           0, 0.0016, -0.4510, 1.5, -0.3528, 0,
           1.4112, 0.2730, -1.3730, -0.9863, -0.0062,
           1.3945, -0.2731, 1.3698, 0.9879, -0.0077;


    q_5 << 0,0,0,0,0,0,
           0.0015, -0.6000,
           0, 0, -1.5300, 1.500, -0.3513, -0.390,
           0, 0.0016, 0.4800, 1.0, -0.3528, 0,
           1.4112, 1.2000, -1.3730, -0.9863, -0.0062,
           1.3945, -1.2000, 1.3698, 0.9879, -0.0077;


//    std::cout<<robot->getJointUpLimits().transpose()<<std::endl;
//    std::cout<<robot->getJointLowLimits().transpose()<<std::endl;


//    geo::InverseDynamics* InvDyn = new geo::InverseDynamics( robot );
//    InvDyn->setGeneralizedCoordinates(q_4,q_4,q_4);
//    InvDyn->computeCenterOfMass(false, false);
//    std::cout<<"CoM = "<<InvDyn->getRobotCoM().transpose()<<std::endl;

//    InvDyn->computeForwardKinematics();

//    geo::DynamicBody * BB = (geo::DynamicBody *)( robot->getPubBodies().at(6) );
//    std::cout<<"LL pos = " <<std::endl<< BB->getName() << std::endl;
//    std::cout<<"LL pos = " <<std::endl<< BB->getGlobalConfiguration() << std::endl;


    //! Fill up stack of constraints
    //! ---------------------------------------------------------------------------------------------------------- !//
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
    optSettings->initialConfiguration = q_4;

    optSettings->StackConstraints.push_back(geo::constraint_finalConfiguration);
    optSettings->finalConfiguration = q_5;

    optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

//    optSettings->StackConstraints.push_back(geo::constraint_pelvisSymmetry);

//    optSettings->StackConstraints.push_back(geo::constraint_centerOfMass);

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
//    app->Options()->SetNumericValue("max_cpu_time", 200);
//    app->Options()->SetIntegerValue("acceptable_iter", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive"); ///monotone adaptive
//    app->Options()->SetStringValue("output_file", "ipopt.out");

    //! limited-memory for BFGS and exact for our analytic
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

//    //! finite-difference-values for numeric aprox. and exact for our analytic
//    app->Options()->SetStringValue("jacobian_approximation", "exact");

//    app->Options()->SetStringValue("derivative_test", "second-order");
//    app->Options()->SetStringValue("jac_c_constant", "yes");
//    app->Options()->SetStringValue("linear_solver", "mumps");
//    app->Options()->SetStringValue("accept_every_trial_step", "yes");
//    app->Options()->SetNumericValue("bound_relax_factor", 0.2);
//    app->Options()->SetStringValue("corrector_type", "primal-dual");
//    app->Options()->SetNumericValue("obj_scaling_factor", 0.5);
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
//    app->Options()->SetStringValue("alpha_for_y", "full");


    //! Initialize the IpoptApplication and process the options
    //! ---------------------------------------------------------------------------------------------------------- !//
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }


    //! Ask Ipopt to solve the problem
    //! ---------------------------------------------------------------------------------------------------------- !//
    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    //! As the SmartPtrs go out of scope, the reference count will be decremented and the objects will automatically be deleted.
    return (int) status;

}
