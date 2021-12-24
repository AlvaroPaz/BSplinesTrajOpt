/**
 *	\file examples/example_04.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the Enhanced ABA Differentiation wrt state
 */

#include "geombd/core.h"
#include "geombd/dynamics.h"

#include <Eigen/Dense>

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>

std::string naoFile = "../../data/NAOURDF_inertial.xml";
//std::string naoFile = "../../data/NAOURDF_inertial_ext.xml"; /// for simulating 28DoF

using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000

//! Define time variables
geo::VectorXr S;
short int n;
int M;
geo::real_t si, sf, s_k;


int main(){

    cout << "Example to test the ABA Differentiation wrt state" << endl;

    //! Build multibody
    geo::World world = geo::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();

    world.loadMultiBodyURDF(sFile,robotId, geo::kNAO);

    std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);

    geo::VectorXr q(robot->getDoF());
    geo::VectorXr dq(robot->getDoF());
    geo::VectorXr ddq(robot->getDoF());
    geo::VectorXr tau(robot->getDoF());

    q.setZero();        robot->setConfiguration(q);
    dq.setZero();       robot->setGeneralizedVelocity(dq);
    ddq.setZero();      robot->setGeneralizedAcceleration(ddq);
    tau.setZero();      robot->setGeneralizedTorques(tau);

    //! Get initial configuration and compute forward kinematics
    robot->setConfiguration(q);


    //! Set time variables
    n = robot->getDoF();     // degrees of freedom
    M = 99;    // number of partitions; number of collocation points = M+1
    si = 0;    // first collocation point
    sf = 1;    // last collocation point


    //! Compute time parameter vector
    S.setLinSpaced(M+1,si,sf);


    //! Vector of coeficients
    geo::VectorXr C;
    C.setLinSpaced(n,0.1,0.1*n);


    //! Forward dynamics object
    geo::ForwardDynamics ForwardDynamics(robot);

    int t_total = 0, t_total_00 = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();


    //! Time loop for forward dynamics S.size()
    for (int k = 0 ; k<S.size() ; k++) {
        s_k = S(k);

        //! Set generalized vectors
        q = C*sin(s_k);
        robot->setConfiguration(q);

        dq = C*cos(s_k);
        robot->setGeneralizedVelocity(dq);

        tau = C*sin(s_k);
        robot->setGeneralizedTorques(tau);

        //! Perform the Enhanced Forward Dynamics
        t1 = std::chrono::high_resolution_clock::now();
        ForwardDynamics.computeEnhancedForwardDynamics(false);
        t2 = std::chrono::high_resolution_clock::now();

        t_total_00 += (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

        //! Perform the Enhanced Forward Dynamics and its partial derivative
        t1 = std::chrono::high_resolution_clock::now();
        ForwardDynamics.computeEnhancedForwardDynamics(true);
        t2 = std::chrono::high_resolution_clock::now();

        t_total += (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    }

    cout<<"Forward dynamics = "<<t_total_00<<endl;
    cout<<"Forward dynamics with derivatives = "<<t_total<<endl;

return 0;
}

