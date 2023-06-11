/**
 *	\file examples/example_01.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the Inverse Dynamics and its two first partial derivatives wrt state
 */

#define EIGEN_NO_DEBUG

#include "geombd/core.h"
#include "geombd/dynamics.h"
#include "geombd/trajectoryOptimization.h"

//std::string naoFile = "../../data/nao_inertial_python.urdf";


#define __FU_PATH_PREFIX__ "../../../data/TROmodels/"
std::string urdf_dir = __FU_PATH_PREFIX__;

//std::string urdf_dir = __FU_PATH_PREFIX__ "lwr.urdf"; //7
//std::string urdf_dir = __FU_PATH_PREFIX__ "nao_inertial_python.urdf"; //24
//std::string urdf_dir = __FU_PATH_PREFIX__ "HRP2.urdf"; //28
//std::string urdf_dir = __FU_PATH_PREFIX__ "atlas.urdf"; //30

using std::cout;
using std::endl;
using namespace Eigen;

#define infty 100000

//! Set time variables
int n, k = 0;

bool firstDerivative, secondDerivative;
int diffSize;
double inc_s = 1e-8;    // sensitivity

geo::VectorXr c, q, dq, ddq;
geo::MatrixXr D_q, D_dq, D_ddq, DD_q;
geo::MatrixXr numericD_CoM, numericDD_CoM, numericD_Mu, numericDD_Mu, numericD_Tau, numericDD_Tau;


int main( int argc, char** argv ){

  if (argc > 1) {
      urdf_dir.append( argv[1] );
    } else {
      urdf_dir.append( "nao_inertial_python.urdf" );
    }

  cout << "//!--------------------------------------------------------------------------------------!//" << endl;
  cout << "Unit test for accuracy in CoM, momentum, inverse dynamics and their derivatives wrt random c" << endl;
  cout << "//!--------------------------------------------------------------------------------------!//" << endl;


  //! Build multibody
  //!------------------------------------------------------------------------------!//
  geo::World world = geo::World();
  int robotId = world.getRobotsVector()->size();
  world.loadMultiBodyURDF(urdf_dir,robotId, geo::kNAO);

  std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);
  n = robot->getDoF();

  geo::VectorXr q = robot->getConfiguration();

  q.setZero();

  robot->setConfiguration(q);

  //! Optimization settings pointer
  //!------------------------------------------------------------------------------!//
  int numberControlPoints = 4;
  int diffSize = n*numberControlPoints;
  robot->setDifferentiationSize( diffSize );

  //! Random c
  //!------------------------------------------------------------------------------!//
  c = geo::VectorXr::Random(diffSize);

  //! Dynamics object pointer
  //!------------------------------------------------------------------------------!//
  auto robotDynamics = std::make_shared< geo::InverseDynamics >( robot );

  //! Build random basis functions
  //!------------------------------------------------------------------------------!//
  geo::MatrixXr D_q, D_dq, D_ddq;
  D_q   = Eigen::kroneckerProduct(geo::RowVectorXr::Random(numberControlPoints),
                                  geo::MatrixXr::Identity(n,n));
  D_dq  = Eigen::kroneckerProduct(geo::RowVectorXr::Random(numberControlPoints),
                                  geo::MatrixXr::Identity(n,n));
  D_ddq = Eigen::kroneckerProduct(geo::RowVectorXr::Random(numberControlPoints),
                                  geo::MatrixXr::Identity(n,n));

  q = D_q * c;  dq = D_dq * c;  ddq = D_ddq * c;

  robotDynamics->setGeneralizedCoordinates(q, dq, ddq);

  robotDynamics->setGeneralizedCoordinatesDifferentiation(D_q, D_dq, D_ddq);

  geo::MatrixXr DD_q(n, diffSize*diffSize), DD_q_Van(n, diffSize*(diffSize+1)/2);
  for(int _x_ = 0; _x_ < n; _x_++) {
      DD_q.row(_x_) = Eigen::kroneckerProduct(D_q.row(_x_), D_q.row(_x_));
      DD_q_Van.row(_x_) = Eigen::vanishProduct(D_q.row(_x_), D_q.row(_x_));
    }

  robotDynamics->setGeneralizedCoordinatesSecondDifferentiation( DD_q );
  robotDynamics->setContractedDD_q( DD_q_Van );


  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Center of Mass First-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//
  robotDynamics->computeCenterOfMass(true, true);

  auto CoM = robotDynamics->getRobotCoM();
  auto D_CoM = robotDynamics->getRobotD_CoM();
  auto DD_CoM = robotDynamics->getRobotDD_CoM();


  //! Numeric verification through finite differences //!
  geo::MatrixXr numericD_CoM(D_CoM.rows(),D_CoM.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeCenterOfMass(false,false);

      auto temporalCoM = robotDynamics->getRobotCoM();

      numericD_CoM.col(iter) = (temporalCoM-CoM)/inc_s;
    }

  cout << endl << "//!----------------------------->> Center of Mass First-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_CoM.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_CoM.cwiseAbs().sum() << endl;
  auto error_in_D_CoM = D_CoM - numericD_CoM;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_CoM.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Center of Mass Second-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//

  //! Numeric verification through finite differences //!
  geo::MatrixXr numericDD_CoM(DD_CoM.rows(),DD_CoM.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeCenterOfMass(true,false);

      auto temporalD_CoM = robotDynamics->getRobotD_CoM();

      numericDD_CoM.middleCols(iter*diffSize,diffSize) = (temporalD_CoM-D_CoM)/inc_s;
    }

  cout << "//!----------------------------->> Center of Mass Second-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_CoM.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_CoM.cwiseAbs().sum() << endl;
  auto error_in_DD_CoM = DD_CoM - numericDD_CoM;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_CoM.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Spatial Momentum First-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//
  robotDynamics->setGeneralizedCoordinates(q, dq, ddq);
  robotDynamics->computeSpatialMomentum(true, true);

  auto SMu = robotDynamics->getSpatialMomentum();
  auto D_SMu = robotDynamics->getD_SpatialMomentum();
  auto DD_SMu = robotDynamics->getDD_SpatialMomentum();


  //! Numeric verification through finite differences //!
  geo::MatrixXr numericD_SMu(D_SMu.rows(),D_SMu.cols());

  inc_s = pow(2,-24);  // -8, -14, -24

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeSpatialMomentum(false,false);

      auto temporalSMu = robotDynamics->getSpatialMomentum();

      numericD_SMu.col(iter) = (temporalSMu-SMu)/inc_s;
    }

  cout << "//!----------------------------->> Spatial Momentum First-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_SMu.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_SMu.cwiseAbs().sum() << endl;
  auto error_in_D_SMu = D_SMu - numericD_SMu;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_SMu.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Spatial Momentum Second-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//

  //! Numeric verification through finite differences //!
  geo::MatrixXr numericDD_SMu(DD_SMu.rows(), DD_SMu.cols());

  inc_s = pow(2,-24);   // 8

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeSpatialMomentum(true,false);

      auto temporalD_SMu = robotDynamics->getD_SpatialMomentum();

      numericDD_SMu.middleCols(iter*diffSize,diffSize) = (temporalD_SMu-D_SMu)/inc_s;
    }

  cout << "//!----------------------------->> Spatial Momentum Second-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_SMu.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_SMu.cwiseAbs().sum() << endl;
  auto error_in_DD_SMu = DD_SMu - numericDD_SMu;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_SMu.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Centroidal Momentum First-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//
  robotDynamics->setGeneralizedCoordinates(q, dq, ddq);
  robotDynamics->computeCentroidalMomentum(true, true);

  auto Mu = robotDynamics->getCentroidalMomentum();
  auto D_Mu = robotDynamics->getD_CentroidalMomentum();
  auto DD_Mu = robotDynamics->getDD_CentroidalMomentum();


  //! Numeric verification through finite differences //!
  geo::MatrixXr numericD_Mu(D_Mu.rows(),D_Mu.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeCentroidalMomentum(false,false);

      auto temporalMu = robotDynamics->getCentroidalMomentum();

      numericD_Mu.col(iter) = (temporalMu-Mu)/inc_s;
    }

  cout << "//!----------------------------->> Centroidal Momentum First-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_Mu.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_Mu.cwiseAbs().sum() << endl;
  auto error_in_D_Mu = D_Mu - numericD_Mu;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_Mu.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Centroidal Momentum Second-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//

  //! Numeric verification through finite differences //!
  geo::MatrixXr numericDD_Mu(DD_Mu.rows(), DD_Mu.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeCentroidalMomentum(true,false);

      auto temporalD_Mu = robotDynamics->getD_CentroidalMomentum();

      numericDD_Mu.middleCols(iter*diffSize,diffSize) = (temporalD_Mu-D_Mu)/inc_s;
    }

  cout << "//!----------------------------->> Centroidal Momentum Second-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_Mu.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_Mu.cwiseAbs().sum() << endl;
  auto error_in_DD_Mu = DD_Mu - numericDD_Mu;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_Mu.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Inverse Dynamics First-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//
  robotDynamics->setGeneralizedCoordinates(q, dq, ddq);
  robotDynamics->computeInverseDynamics(true, true);

  auto Tau = robotDynamics->getGeneralizedTorques();
  auto D_Tau = robotDynamics->getGeneralizedTorquesFirstDifferentiation();
  auto DD_Tau = robotDynamics->getGeneralizedTorquesSecondDifferentiation();


  //! Numeric verification through finite differences //!
  geo::MatrixXr numericD_Tau(D_Tau.rows(),D_Tau.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeInverseDynamics(false,false);

      auto temporalTau = robotDynamics->getGeneralizedTorques();

      numericD_Tau.col(iter) = (temporalTau-Tau)/inc_s;
    }

  cout << "//!----------------------------->> Inverse Dynamics First-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_Tau.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_Tau.cwiseAbs().sum() << endl;
  auto error_in_D_Tau = D_Tau - numericD_Tau;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_Tau.cwiseAbs().sum() << endl << endl;



  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Inverse Dynamics Second-Order Derivative wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//

  //! Numeric verification through finite differences //!
  geo::MatrixXr numericDD_Tau(DD_Tau.rows(), DD_Tau.cols());

  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      auto q_aux1 = q;   auto dq_aux1 = dq;   auto ddq_aux1 = ddq;

      q_aux1 = D_q * c_aux1;   dq_aux1 = D_dq * c_aux1;   ddq_aux1 = D_ddq * c_aux1;

      robotDynamics->setGeneralizedCoordinates(q_aux1, dq_aux1, ddq_aux1);
      robotDynamics->computeInverseDynamics(true,false);

      auto temporalD_Tau = robotDynamics->getGeneralizedTorquesFirstDifferentiation();

      numericDD_Tau.middleCols(iter*diffSize,diffSize) = (temporalD_Tau-D_Tau)/inc_s;
    }

  cout << "//!----------------------------->> Inverse Dynamics Second-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_Tau.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_Tau.cwiseAbs().sum() << endl;
  auto error_in_DD_Tau = DD_Tau - numericDD_Tau;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_Tau.cwiseAbs().sum() << endl << endl << endl;



  return 0;
}

