/**
 *	\file examples/example_03.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the ABA Differentiation wrt state
 */

#include "geombd/core.h"
#include "geombd/dynamics.h"

#include <Eigen/Dense>

#include <iomanip>
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

//! Define time variables
hr::VectorXr S;
short int n;
int M;
hr::real_t si, sf, s_k;


int main(){

    cout << "Example to test the ABA Differentiation wrt state" << endl;

    //! Files to write up
    std::ofstream myfile;
    myfile.open ("inv_H.txt");

    //! Build multibody
    hr::core::World world = hr::core::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, hr::core::kNAO);

    std::shared_ptr<hr::core::MultiBody> robot = world.getRobot(0);


    int t_total = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();


    int n = robot->getDoF();
    hr::VectorXr q(n), dq(n), tau(n);

    q.setRandom();    dq.setRandom();    tau.setRandom();


    hr::core::ForwardDynamics robotMotion(robot);
    robot->setConfiguration(q);
    robot->setGeneralizedVelocity(dq);
    robot->setGeneralizedTorques(tau);


    //! Forward Dynamics Gradient wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeForwardDynamics(true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Forward Dynamics Gradient time = " << t_total << endl;


    auto ddq = robotMotion.getJointAcceleration();
    auto D_ddq = robotMotion.getDiffJointAcceleration();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericGradient(D_ddq.rows(),D_ddq.cols());
    hr::real_t inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robot->setConfiguration(q_aux1);
        robot->setGeneralizedVelocity(dq);
        robot->setGeneralizedTorques(tau);

        robotMotion.computeForwardDynamics(false);

        auto temporal_ddq = robotMotion.getJointAcceleration();

        numericGradient.col(iter) = (temporal_ddq-ddq)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robot->setConfiguration(q);
        robot->setGeneralizedVelocity(dq_aux1);
        robot->setGeneralizedTorques(tau);

        robotMotion.computeForwardDynamics(false);

        temporal_ddq = robotMotion.getJointAcceleration();

        numericGradient.col(iter+n) = (temporal_ddq-ddq)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Forward Dynamics Gradient time = " << t_total << endl;

    cout << "-------------------------- DDQ GRADIENT" << endl;
    cout << "Analytic: " << D_ddq.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericGradient.cwiseAbs().sum() << endl;
    auto error_in_D_ddq = D_ddq - numericGradient;
    cout << "Error: " << error_in_D_ddq.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////

    //! Inertia Matrix Inverse

    robot->setConfiguration(q);
    robot->setGeneralizedVelocity(dq);
    robot->setGeneralizedTorques(tau);

    //! Perform the inverse of the inertia matrix
    robotMotion.computeInverseInertiaMatrix();

    //! Retrieve the inverse of the inertia matrix
    auto inv_H = robotMotion.getInverseInertiaMatrix();

    //! Write up the txt files with high accuracy
    myfile << std::scientific << std::setprecision(20) << inv_H;

return 0;
}
