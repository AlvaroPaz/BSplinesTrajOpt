/**
 *	\file examples/example_01.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the Inverse Dynamics and its two first partial derivatives wrt state
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

    cout << "Example to test the Inverse Dynamics and its two first partial derivatives wrt state" << endl;

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
    hr::VectorXr q(n), dq(n), ddq(n);

    q.setRandom();    dq.setRandom();    ddq.setRandom();


    hr::core::InverseDynamics robotMotion(robot);
    robotMotion.setGeneralizedCoordinates(q, dq, ddq);


    //! Inverse Dynamics Gradient wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeInverseDynamics(true, false);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Inverse Dynamics Gradient time = " << t_total << endl;



    auto Tau = robotMotion.getGeneralizedTorques();
    auto D_Tau = robotMotion.getGeneralizedTorquesFirstDifferentiation();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericGradient(D_Tau.rows(),D_Tau.cols());
    hr::real_t inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q_aux1, dq, ddq);
        robotMotion.computeInverseDynamics(false, false);

        auto temporalTau = robotMotion.getGeneralizedTorques();

        numericGradient.col(iter) = (temporalTau-Tau)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q, dq_aux1, ddq);
        robotMotion.computeInverseDynamics(false, false);

        temporalTau = robotMotion.getGeneralizedTorques();

        numericGradient.col(iter+n) = (temporalTau-Tau)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Inverse Dynamics Gradient time = " << t_total << endl;

    cout << "-------------------------- TAU GRADIENT" << endl;
    cout << "Analytic: " << D_Tau.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericGradient.cwiseAbs().sum() << endl;
    auto error_in_D = D_Tau - numericGradient;
    cout << "Error: " << error_in_D.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Inverse Dynamics Hessian wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeInverseDynamics(true, true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Inverse Dynamics Hessian time = " << t_total << endl;



    D_Tau = robotMotion.getGeneralizedTorquesFirstDifferentiation();
    auto DD_Tau = robotMotion.getGeneralizedTorquesSecondDifferentiation();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericHessian(DD_Tau.rows(),DD_Tau.cols());

    inc_s = pow(2,-14);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q_aux1, dq, ddq);
        robotMotion.computeInverseDynamics(true, false);

        auto temporalD_Tau = robotMotion.getGeneralizedTorquesFirstDifferentiation();

        numericHessian.middleCols(2*n*iter,2*n) = (temporalD_Tau-D_Tau)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q, dq_aux1, ddq);
        robotMotion.computeInverseDynamics(true, false);

        temporalD_Tau = robotMotion.getGeneralizedTorquesFirstDifferentiation();

        numericHessian.middleCols(2*n*(iter+n),2*n) = (temporalD_Tau-D_Tau)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Inverse Dynamics Hessian time = " << t_total << endl;

    cout << "-------------------------- TAU HESSIAN" << endl;
    cout << "Analytic: " << DD_Tau.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericHessian.cwiseAbs().sum() << endl;
    auto error_in_DD = DD_Tau - numericHessian;
    cout << "Error: " << error_in_DD.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Center of Mass wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeCenterOfMass(true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Center of Mass Gradient time = " << t_total << endl;


    auto CoM = robotMotion.getMultibodyCoM();
    auto D_CoM = robotMotion.getMultibodyCoMDifferentiation();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericD_CoM(D_CoM.rows(),D_CoM.cols());

    inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q_aux1, dq, ddq);
        robotMotion.computeCenterOfMass(false);

        auto temporalCoM = robotMotion.getMultibodyCoM();

        numericD_CoM.col(iter) = (temporalCoM-CoM)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q, dq_aux1, ddq);
        robotMotion.computeCenterOfMass(false);

        temporalCoM = robotMotion.getMultibodyCoM();

        numericD_CoM.col(iter+n) = (temporalCoM-CoM)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Center of Mass Gradient time = " << t_total << endl;

    cout << "-------------------------- COM GRADIENT" << endl;
    cout << "Analytic: " << D_CoM.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericD_CoM.cwiseAbs().sum() << endl;
    auto error_in_D_CoM = D_CoM - numericD_CoM;
    cout << "Error: " << error_in_D_CoM.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Spatial Momentum wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeSpatialMomentum(true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Spatial Momentum Gradient time = " << t_total << endl;


    auto SMu = robotMotion.getSpatialMomentum();
    auto D_SMu = robotMotion.getSpatialMomentumDifferentiation();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericD_SMu(D_SMu.rows(),D_SMu.cols());

    inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q_aux1, dq, ddq);
        robotMotion.computeSpatialMomentum(false);

        auto temporalSMu = robotMotion.getSpatialMomentum();

        numericD_SMu.col(iter) = (temporalSMu-SMu)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q, dq_aux1, ddq);
        robotMotion.computeSpatialMomentum(false);

        temporalSMu = robotMotion.getSpatialMomentum();

        numericD_SMu.col(iter+n) = (temporalSMu-SMu)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Spatial Momentum Gradient time = " << t_total << endl;

    cout << "-------------------------- SMu GRADIENT" << endl;
    cout << "Analytic: " << D_SMu.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericD_SMu.cwiseAbs().sum() << endl;
    auto error_in_D_SMu = D_SMu - numericD_SMu;
    cout << "Error: " << error_in_D_SMu.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;


    ////////////////////////////////////////////////////////////////////


    //! Centroidal Momentum wrt State Verification

    t1 = std::chrono::high_resolution_clock::now();

    robotMotion.computeCentroidalMomentum(true);

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Analytic Centroidal Momentum Gradient time = " << t_total << endl;


    auto Mu = robotMotion.getCentroidalMomentum();
    auto D_Mu = robotMotion.getCentroidalMomentumDifferentiation();


    //! Numeric verification through finite differences //!
    hr::MatrixXr numericD_Mu(D_Mu.rows(),D_Mu.cols());

    inc_s = pow(2,-24);

    t1 = std::chrono::high_resolution_clock::now();

    for(short int iter = 0 ; iter < n ; iter++){
        //! wrt q
        auto q_aux1 = q;
        q_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q_aux1, dq, ddq);
        robotMotion.computeCentroidalMomentum(false);

        auto temporalMu = robotMotion.getCentroidalMomentum();

        numericD_Mu.col(iter) = (temporalMu-Mu)/inc_s;

        //! wrt dq
        auto dq_aux1 = dq;
        dq_aux1(iter) += inc_s;

        robotMotion.setGeneralizedCoordinates(q, dq_aux1, ddq);
        robotMotion.computeCentroidalMomentum(false);

        temporalMu = robotMotion.getCentroidalMomentum();

        numericD_Mu.col(iter+n) = (temporalMu-Mu)/inc_s;
    }

    t2 = std::chrono::high_resolution_clock::now();

    t_total = (int) std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    cout << "Numeric Centroidal Momentum Gradient time = " << t_total << endl;

    cout << "-------------------------- MU GRADIENT" << endl;
    cout << "Analytic: " << D_Mu.cwiseAbs().sum() << endl;
    cout << "Numeric: " << numericD_Mu.cwiseAbs().sum() << endl;
    auto error_in_D_Mu = D_Mu - numericD_Mu;
    cout << "Error: " << error_in_D_Mu.cwiseAbs().sum() << endl;
    cout << "--------------------------" << endl;

return 0;
}

