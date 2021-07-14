/**
 *	\file examples/example_02.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the cost function, its gradient and its Hessian wrt control points
 */

#include "geombd/core.h"
#include "geombd/dynamics.h"
#include "geombd/trajectoryOptimization.h"

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

std::string naoFile = "../../data/NAOURDF_inertial.xml";

using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000


int main(){

    cout << "Example to test the Inverse Dynamics and its two first partial derivatives wrt control points" << endl;

    //! Build multibody
    hr::core::World world = hr::core::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, hr::core::kNAO);


    std::shared_ptr< hr::core::MultiBody > robot = world.getRobot(0);

    hr::VectorXr q = robot->getConfiguration();

    q.setZero();

    robot->setConfiguration(q);


    //! Optimization parameters
    std::shared_ptr< hr::core::robotSettingsTrajectoryOptimization > optSettings(new hr::core::robotSettingsTrajectoryOptimization);

    optSettings->n =  robot->getDoF();
    optSettings->numberControlPoints = 4;//4
    optSettings->numberPartitions    = 7;//7
    optSettings->si = 0.0;
    optSettings->sf = 0.1;
    optSettings->S = hr::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
    optSettings->DifferentiationWRT = hr::core::wrt_controlPoints;
    robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );

    //! Fill up stack of constraints
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(hr::core::constraint_initialConfiguration);
    optSettings->initialConfiguration = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_finalConfiguration);
    optSettings->finalConfiguration = hr::VectorXr::Ones(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_pelvisSymmetry);
    optSettings->StackConstraints.push_back(hr::core::constraint_centerOfMass);
    optSettings->StackConstraints.push_back(hr::core::constraint_centroidalMomentum);



    int t_total = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();


    //! Trajectory optimization instantiation
    hr::core::DirectCollocation robotNonlinearProblem(robot, optSettings);


    //! Build basis functions
    robotNonlinearProblem.buildBasisFunctions(  );


    hr::VectorXr c(optSettings->n*optSettings->numberControlPoints); // control points

    c.setRandom();

    //! Objective Function Gradient wrt Control Points Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotNonlinearProblem.computeObjectiveFunction(c, true, false);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Objective Function Gradient time = " << t_total << endl;



    auto cost = robotNonlinearProblem.getCost();
    auto D_cost = robotNonlinearProblem.getCostGradient();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericD_Cost(D_cost.rows(),D_cost.cols());
    hr::real_t inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < c.size() ; iter++){
        //! wrt c
        auto c_aux1 = c;
        c_aux1(iter) += inc_s;

        robotNonlinearProblem.computeObjectiveFunction(c_aux1, false, false);

        auto temporalCost = robotNonlinearProblem.getCost();

        numericD_Cost(iter) = (temporalCost-cost)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Objective Function Gradient time = " << t_total << endl;

    cout << "-------------------------- COST GRADIENT" << endl;
    cout << "Analytic: " << D_cost.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericD_Cost.cwiseAbs().sum() << endl;
    auto error_in_D_Cost = D_cost - numericD_Cost;
    cout << "Error: " << error_in_D_Cost.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Objective Function Hessian wrt Control Points Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotNonlinearProblem.computeObjectiveFunction(c, true, true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Objective Function Hessian time = " << t_total << endl;


    auto DD_cost = robotNonlinearProblem.getCostHessian();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericDD_Cost(DD_cost.rows(),DD_cost.cols());
    inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < c.size() ; iter++){
        //! wrt c
        auto c_aux1 = c;
        c_aux1(iter) += inc_s;

        robotNonlinearProblem.computeObjectiveFunction(c_aux1, true, false);

        auto temporalD_Cost = robotNonlinearProblem.getCostGradient();

        numericDD_Cost.row(iter) = (temporalD_Cost-D_cost)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Objective Function Hessian time = " << t_total << endl;

    cout << "-------------------------- COST HESSIAN" << endl;
    cout << "Analytic: " << DD_cost.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericDD_Cost.cwiseAbs().sum() << endl;
    auto error_in_DD_Cost = DD_cost - numericDD_Cost;
    cout << "Error: " << error_in_DD_Cost.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Constraints Jacobian wrt Control Points Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotNonlinearProblem.computeConstraints(c, true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Constraints Jacobian time = " << t_total << endl;


    auto Constraints = robotNonlinearProblem.getConstraints();
    auto D_constraints = robotNonlinearProblem.getConstraintsJacobian();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericD_constraints(D_constraints.rows(),D_constraints.cols());
    inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < c.size() ; iter++){
        //! wrt c
        auto c_aux1 = c;
        c_aux1(iter) += inc_s;

        robotNonlinearProblem.computeConstraints(c_aux1, false);

        auto temporalConstraints = robotNonlinearProblem.getConstraints();

        numericD_constraints.col(iter) = (temporalConstraints-Constraints)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Constraints Jacobian time = " << t_total << endl;

    cout << "-------------------------- CONSTRAINTS JACOBIAN" << endl;
    cout << "Analytic: " << D_constraints.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericD_constraints.cwiseAbs().sum() << endl;
    auto error_in_D_constraints = D_constraints - numericD_constraints;
    cout << "Error: " << error_in_D_constraints.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////

    //!   MODIFY THE CONSTRAINTS   !//

    //! Fill up stack of constraints
    optSettings->StackConstraints.clear();

    optSettings->StackConstraints.push_back(hr::core::constraint_finalGeneralizedVelocity);
    optSettings->finalGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_initialConfiguration);
    optSettings->initialConfiguration = hr::VectorXr::Zero(optSettings->n,1);

    optSettings->StackConstraints.push_back(hr::core::constraint_initialGeneralizedVelocity);
    optSettings->initialGeneralizedVelocity = hr::VectorXr::Zero(optSettings->n,1);


    //! Update last changes
    robotNonlinearProblem.updateSettings();

return 0;
}
