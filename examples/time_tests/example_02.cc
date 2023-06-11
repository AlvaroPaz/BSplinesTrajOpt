/**
 *	\file examples/example_02.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2021
 *
 *	Example to test the dynamic objects and their two first partial (numeric) derivatives wrt to control points with random c
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
const int M = 1e2;      // sample size
int n, k = 0;

bool firstDerivative, secondDerivative;
int diffSize;
double inc_c = 1e-8;    // sensitivity

geo::VectorXr c, q, dq, ddq, c_aux, q_aux, dq_aux, ddq_aux, c_aux_in;
geo::MatrixXr D_q, D_dq, D_ddq, DD_q;
geo::MatrixXr numericD_CoM, numericDD_CoM, numericD_Mu, numericDD_Mu, numericD_Tau, numericDD_Tau;


//! Time loop for center of mass
//!------------------------------------------------------------------------------!//
void loop_CoM ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeCenterOfMass(false, false);

      if (firstDerivative) {
          auto CoM_exa = robotDynamics->getRobotCoM();
          numericD_CoM = geo::MatrixXr::Zero(3, diffSize);

          for(short int iter = 0 ; iter < diffSize ; iter++){
              //! wrt c
              c_aux = c;
              c_aux(iter) += inc_c;

              q_aux   = D_q*c_aux;
              dq_aux  = D_dq*c_aux;
              ddq_aux = D_ddq*c_aux;
              robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
              robotDynamics->computeCenterOfMass(false, false);

              auto temporalCoM = robotDynamics->getRobotCoM();

              numericD_CoM.col(iter) = (temporalCoM-CoM_exa)/inc_c;
            }

          if (secondDerivative) {

              //! f(x+h) and f(x+k) evaluation
              double h__ = 1e-4;
              double k__ = 1e-4;
              double hk__ = h__*k__;

              geo::MatrixXr f_h(3,diffSize), f_k(3,diffSize);

              for(short int iter = 0 ; iter < diffSize ; iter++){
                  //! h variation
                  c_aux = c;
                  c_aux(iter) += h__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeCenterOfMass(false, false);

                  f_h.col(iter) = robotDynamics->getRobotCoM();

                  //! k variation
                  c_aux = c;
                  c_aux(iter) += k__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeCenterOfMass(false, false);

                  f_k.col(iter) = robotDynamics->getRobotCoM();
                }

              geo::MatrixXr temporalD_CoM(3, diffSize);
              numericDD_CoM = geo::MatrixXr::Zero(3, diffSize*diffSize);

              for(int i__ = 0 ; i__ < diffSize ; i__++){
                  //! h variation
                  auto c_aux = c;
                  c_aux(i__) += h__;

                  for(int j__ = 0 ; j__ < diffSize ; j__++){
                      //! k variation
                      c_aux_in = c_aux;
                      c_aux_in(j__) += k__;

                      //! f(x+h+k) evaluation
                      q_aux   = D_q*c_aux_in;
                      dq_aux  = D_dq*c_aux_in;
                      ddq_aux = D_ddq*c_aux_in;
                      robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                      robotDynamics->computeCenterOfMass(false, false);

                      geo::VectorXr f_h_k = robotDynamics->getRobotCoM();

                      temporalD_CoM.col(j__) = (f_h_k - f_h.col(i__) - f_k.col(j__) + CoM_exa)/hk__;
                    }
                  numericDD_CoM.middleCols(i__*diffSize, diffSize) = temporalD_CoM;
                }
            }
        }
    }
}


//! Time loop for centroidal momentum
//!------------------------------------------------------------------------------!//
void loop_Mu ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeCentroidalMomentum(false, false);

      if (firstDerivative) {
          auto Mu_exa = robotDynamics->getCentroidalMomentum();
          numericD_Mu = geo::MatrixXr::Zero(6, diffSize);

          for(short int iter = 0 ; iter < diffSize ; iter++){
              //! wrt c
              auto c_aux = c;
              c_aux(iter) += inc_c;

              auto q_aux   = D_q*c_aux;
              auto dq_aux  = D_dq*c_aux;
              auto ddq_aux = D_ddq*c_aux;
              robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
              robotDynamics->computeCentroidalMomentum(false, false);

              auto temporalMu = robotDynamics->getCentroidalMomentum();

              numericD_Mu.col(iter) = (temporalMu-Mu_exa)/inc_c;
            }

          if (secondDerivative) {

              //! f(x+h) and f(x+k) evaluation
              double h__ = 1e-4;
              double k__ = 1e-4;
              double hk__ = h__*k__;

              geo::MatrixXr f_h(6,diffSize), f_k(6,diffSize);

              for(short int iter = 0 ; iter < diffSize ; iter++){
                  //! h variation
                  c_aux = c;
                  c_aux(iter) += h__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeCentroidalMomentum(false, false);

                  f_h.col(iter) = robotDynamics->getCentroidalMomentum();

                  //! k variation
                  c_aux = c;
                  c_aux(iter) += k__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeCentroidalMomentum(false, false);

                  f_k.col(iter) = robotDynamics->getCentroidalMomentum();
                }

              geo::MatrixXr temporalD_Mu(6, diffSize);
              numericDD_Mu = geo::MatrixXr::Zero(6, diffSize*diffSize);

              for(int i__ = 0 ; i__ < diffSize ; i__++){
                  //! h variation
                  auto c_aux = c;
                  c_aux(i__) += h__;

                  for(int j__ = 0 ; j__ < diffSize ; j__++){
                      //! k variation
                      c_aux_in = c_aux;
                      c_aux_in(j__) += k__;

                      //! f(x+h+k) evaluation
                      q_aux   = D_q*c_aux_in;
                      dq_aux  = D_dq*c_aux_in;
                      ddq_aux = D_ddq*c_aux_in;
                      robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                      robotDynamics->computeCentroidalMomentum(false, false);

                      geo::VectorXr f_h_k = robotDynamics->getCentroidalMomentum();

                      temporalD_Mu.col(j__) = (f_h_k - f_h.col(i__) - f_k.col(j__) + Mu_exa)/hk__;
                    }
                  numericDD_Mu.middleCols(i__*diffSize, diffSize) = temporalD_Mu;
                }
            }
        }
    }
}


//! Time loop for inverse dynamics
//!------------------------------------------------------------------------------!//
void loop_InvDyn ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeInverseDynamics(false, false);

      if (firstDerivative) {
          auto Tau_exa = robotDynamics->getGeneralizedTorques();
          numericD_Tau = geo::MatrixXr::Zero(n, diffSize);

          for(short int iter = 0 ; iter < diffSize ; iter++){
              //! wrt c
              auto c_aux = c;
              c_aux(iter) += inc_c;

              auto q_aux   = D_q*c_aux;
              auto dq_aux  = D_dq*c_aux;
              auto ddq_aux = D_ddq*c_aux;
              robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
              robotDynamics->computeInverseDynamics(false, false);

              auto temporalTau = robotDynamics->getGeneralizedTorques();

              numericD_Tau.col(iter) = (temporalTau-Tau_exa)/inc_c;
            }

          if (secondDerivative) {

              //! f(x+h) and f(x+k) evaluation
              double h__ = 1e-4;
              double k__ = 1e-4;
              double hk__ = h__*k__;

              geo::MatrixXr f_h(n,diffSize), f_k(n,diffSize);

              for(short int iter = 0 ; iter < diffSize ; iter++){
                  //! h variation
                  c_aux = c;
                  c_aux(iter) += h__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeInverseDynamics(false, false);

                  f_h.col(iter) = robotDynamics->getGeneralizedTorques();

                  //! k variation
                  c_aux = c;
                  c_aux(iter) += k__;

                  //! functions evaluation
                  q_aux   = D_q*c_aux;
                  dq_aux  = D_dq*c_aux;
                  ddq_aux = D_ddq*c_aux;
                  robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                  robotDynamics->computeInverseDynamics(false, false);

                  f_k.col(iter) = robotDynamics->getGeneralizedTorques();
                }

              geo::MatrixXr temporalD_Tau(n, diffSize);
              numericDD_Tau = geo::MatrixXr::Zero(n, diffSize*diffSize);

              for(int i__ = 0 ; i__ < diffSize ; i__++){
                  //! h variation
                  auto c_aux = c;
                  c_aux(i__) += h__;

                  for(int j__ = 0 ; j__ < diffSize ; j__++){
                      //! k variation
                      c_aux_in = c_aux;
                      c_aux_in(j__) += k__;

                      //! f(x+h+k) evaluation
                      q_aux   = D_q*c_aux_in;
                      dq_aux  = D_dq*c_aux_in;
                      ddq_aux = D_ddq*c_aux_in;
                      robotDynamics->setGeneralizedCoordinates(q_aux, dq_aux, ddq_aux);
                      robotDynamics->computeInverseDynamics(false, false);

                      geo::VectorXr f_h_k = robotDynamics->getGeneralizedTorques();

                      temporalD_Tau.col(j__) = (f_h_k - f_h.col(i__) - f_k.col(j__) + Tau_exa)/hk__;
                    }
                  numericDD_Tau.middleCols(i__*diffSize, diffSize) = temporalD_Tau;
                }
            }
        }
    }
}


int main( int argc, char** argv ){
  if (argc > 1) {
      urdf_dir.append( argv[1] );
    } else {
      urdf_dir.append( "nao_inertial_python.urdf" );
    }

  cout << "//!------------------------------------------------------------------------------------!//" << endl;
  cout << "Unit test of dynamic-objects time and their numeric (finite differences) derivatives wrt c" << endl;
  cout << "//!------------------------------------------------------------------------------------!//" << endl;

  //! Build multibody
  //!------------------------------------------------------------------------------!//
  geo::World world = geo::World();
  int robotId = world.getRobotsVector()->size();
  world.loadMultiBodyURDF(urdf_dir,robotId, geo::kNAO);

  std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);
  n = robot->getDoF();

  //! Optimization settings pointer
  //!------------------------------------------------------------------------------!//
  int numberControlPoints = 2;
  diffSize = n*numberControlPoints;
  robot->setDifferentiationSize( diffSize );

  //! Build control-points vectors
  //!------------------------------------------------------------------------------!//
  c = geo::VectorXr::LinSpaced(diffSize, 1, 2*3.1416);

  //! Dynamics object pointer
  //!------------------------------------------------------------------------------!//
  auto robotDynamics = std::make_shared< geo::InverseDynamics >( robot );

  //! Build random basis functions
  //!------------------------------------------------------------------------------!//
  D_q   = geo::MatrixXr::Random(n, diffSize);
  D_dq  = geo::MatrixXr::Random(n, diffSize);
  D_ddq = geo::MatrixXr::Random(n, diffSize);

  //! Set generalized coordinates and derivatives
  //!------------------------------------------------------------------------------!//
  q   = D_q*c;
  dq  = D_dq*c;
  ddq = D_ddq*c;
  robotDynamics->setGeneralizedCoordinates(q, dq, ddq);
  robotDynamics->setGeneralizedCoordinatesDifferentiation(D_q, D_dq, D_ddq);

  DD_q   = geo::MatrixXr::Random(n, diffSize*diffSize);;
  robotDynamics->setGeneralizedCoordinatesSecondDifferentiation( DD_q );

  //! First calls
  //!------------------------------------------------------------------------------!//
  robotDynamics->computeCenterOfMass(true, true);
  auto CoM = robotDynamics->getRobotCoM();
  auto D_CoM = robotDynamics->getRobotD_CoM();
  auto DD_CoM = robotDynamics->getRobotDD_CoM();

  robotDynamics->computeCentroidalMomentum(true, true);
  auto Mu = robotDynamics->getCentroidalMomentum();
  auto D_Mu = robotDynamics->getD_CentroidalMomentum();
  auto DD_Mu = robotDynamics->getDD_CentroidalMomentum();

  robotDynamics->computeInverseDynamics(true, true);
  auto Tau = robotDynamics->getGeneralizedTorques();
  auto D_Tau = robotDynamics->getGeneralizedTorquesFirstDifferentiation();
  auto DD_Tau = robotDynamics->getGeneralizedTorquesSecondDifferentiation();

  //! Time settings
  //!------------------------------------------------------------------------------!//
  int t_total = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  auto t2 = std::chrono::high_resolution_clock::now();
  typedef std::chrono::microseconds time_preci;

  //!------------------------------------------------------------------------------!//
  //!                               Center of Mass                                 !//
  //!------------------------------------------------------------------------------!//

  //! Perform the center of mass loop
  //!------------------------------------------------------------------------------!//
  firstDerivative = false;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_CoM( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << endl << "Center of mass = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the center of mass loop enabling the first derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_CoM( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Center of mass + D = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the center of mass loop enabling the second derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = true;
  t1 = std::chrono::high_resolution_clock::now();
  loop_CoM( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Center of mass + D + DD = " << t_total/M << " microseconds" << endl << endl;

  //!------------------------------------------------------------------------------!//
  //!                            Centroidal Momentum                               !//
  //!------------------------------------------------------------------------------!//

  //! Perform the centroidal momentum loop
  //!------------------------------------------------------------------------------!//
  firstDerivative = false;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_Mu( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Centroidal momentum = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the centroidal momentum loop enabling the first derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_Mu( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Centroidal momentum + D = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the centroidal momentum loop enabling the second derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = true;
  t1 = std::chrono::high_resolution_clock::now();
  loop_Mu( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Centroidal momentum + D + DD = " << t_total/M << " microseconds" << endl << endl;

  //!------------------------------------------------------------------------------!//
  //!                              Inverse Dynamics                                !//
  //!------------------------------------------------------------------------------!//

  //! Perform the inverse dynamics loop
  //!------------------------------------------------------------------------------!//
  firstDerivative = false;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_InvDyn( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Inverse dynamics = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the inverse dynamics loop enabling the first derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = false;
  t1 = std::chrono::high_resolution_clock::now();
  loop_InvDyn( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Inverse dynamics + D = " << t_total/M << " microseconds" << endl;

  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//

  //! Perform the inverse dynamics loop enabling the second derivative
  //!------------------------------------------------------------------------------!//
  firstDerivative = true;   secondDerivative = true;
  t1 = std::chrono::high_resolution_clock::now();
  loop_InvDyn( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  //! Time casting
  //!------------------------------------------------------------------------------!//
  t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
  cout << "Inverse dynamics + D + DD = " << t_total/M << " microseconds" << endl << endl << endl;


  return 0;
}
