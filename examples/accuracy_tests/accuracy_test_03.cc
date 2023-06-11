/**
 *	\file examples/example_07.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2021
 *
 *	Sparsity exploiting test
 */

#define EIGEN_NO_DEBUG

#include "geombd/core.h"
#include "geombd/dynamics.h"

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
const int M = 1e2;    // sample size 1e2
int k = 0;

bool firstDerivative, secondDerivative;

geo::VectorXr q, dq, ddq;


//! Time loop for center of mass
//!------------------------------------------------------------------------------!//
void loop_CoM ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeCenterOfMassII(firstDerivative, secondDerivative);
    }
}


//! Time loop for contracted center of mass
//!------------------------------------------------------------------------------!//
void loop_CoM_contracted ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeContractedCenterOfMass(firstDerivative, secondDerivative);
    }
}


//! Time loop for centroidal momentum
//!------------------------------------------------------------------------------!//
void loop_Mu ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeCentroidalMomentumII(firstDerivative, secondDerivative);
    }
}


//! Time loop for centroidal momentum
//!------------------------------------------------------------------------------!//
void loop_Mu_contracted ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeContractedCentroidalMomentum(firstDerivative, secondDerivative);
    }
}


//! Time loop for inverse dynamics
//!------------------------------------------------------------------------------!//
void loop_InvDyn_contracted ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeContractedInverseDynamics(firstDerivative, secondDerivative);
    }
}


//! Time loop for inverse dynamics
//!------------------------------------------------------------------------------!//
void loop_InvDyn ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeInverseDynamics(firstDerivative, secondDerivative);
    }
}


int main( int argc, char** argv ){
    if (argc > 1) {
        urdf_dir.append( argv[1] );
      } else {
        urdf_dir.append( "nao_inertial_python.urdf" );  // delete urdf once verifyed
      }

    cout << "//!-----------------------------------------------------------------------------!//" << endl;
    cout << "Unit test for second-order derivative of dynamic-objects when symmetry is exploited" << endl;
    cout << "//!-----------------------------------------------------------------------------!//" << endl;

    //! Build multibody
    //!------------------------------------------------------------------------------!//
    geo::World world = geo::World();
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(urdf_dir,robotId, geo::kNAO);

    std::shared_ptr<geo::MultiBody> robot = world.getRobot(0);
    int n = robot->getDoF();

    //! Build input vectors
    //!------------------------------------------------------------------------------!//
    q   = geo::VectorXr::LinSpaced(n, 1, 2*3.1416);
    dq  = geo::VectorXr::LinSpaced(n, 1, 2*3.1416);
    ddq = geo::VectorXr::LinSpaced(n, 1, 2*3.1416);

    //! Optimization settings pointer
    //!------------------------------------------------------------------------------!//
    int numberControlPoints = 2;
    int diffSize = n*numberControlPoints;
    robot->setDifferentiationSize( diffSize );

    //! Dynamics object pointer
    //!------------------------------------------------------------------------------!//
    auto robotDynamics = std::make_shared< geo::InverseDynamics >( robot );
    robotDynamics->setGeneralizedCoordinates(q, dq, ddq);

    //! Build random basis functions
    //!------------------------------------------------------------------------------!//
    geo::MatrixXr D_q, D_dq, D_ddq;
    D_q   = Eigen::kroneckerProduct(geo::RowVectorXr::LinSpaced(numberControlPoints, 1, 1*3.1416),
                                    geo::MatrixXr::Identity(n,n));
    D_dq  = Eigen::kroneckerProduct(geo::RowVectorXr::LinSpaced(numberControlPoints, 1, 2*3.1416),
                                    geo::MatrixXr::Identity(n,n));
    D_ddq = Eigen::kroneckerProduct(geo::RowVectorXr::LinSpaced(numberControlPoints, 1, 3*3.1416),
                                    geo::MatrixXr::Identity(n,n));
    robotDynamics->setGeneralizedCoordinatesDifferentiation(D_q, D_dq, D_ddq);

    geo::MatrixXr DD_q(n, diffSize*diffSize), DD_q_Van(n, diffSize*(diffSize+1)/2);
    for(int _x_ = 0; _x_ < n; _x_++) {
        DD_q.row(_x_) = Eigen::kroneckerProduct(D_q.row(_x_), D_q.row(_x_));
        DD_q_Van.row(_x_) = Eigen::vanishProduct(D_q.row(_x_), D_q.row(_x_));
      }

    robotDynamics->setGeneralizedCoordinatesSecondDifferentiation( DD_q );
    robotDynamics->setContractedDD_q( DD_q_Van );

    //! Time settings
    //!------------------------------------------------------------------------------!//
    int t_total = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    typedef std::chrono::microseconds time_preci;

    //!------------------------------------------------------------------------------!//
    //!                               Center of Mass                                 !//
    //!------------------------------------------------------------------------------!//

    //! Perform the center of mass loop enabling the second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_CoM( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_CoM = robotDynamics->getRobotDD_CoM();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << endl;
    cout << "Second-order derivative of Center of Mass                             =  " << t_total/M << " microseconds" << endl;
    double first = t_total/M;

    //! Perform the center of mass loop enabling the CONTRACTED second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_CoM_contracted( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_CoM_contracted = robotDynamics->getRobotDD_ContractedCoM();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Second-order derivative of Center of Mass when symmetry is exploited  =  " << t_total/M << " microseconds" << endl;
    double second = t_total/M;

    cout << "Improvement rate                                                      =  " << 100 - (second*100)/first << " % faster" << endl << endl;


    //!------------------------------------------------------------------------------!//
    //!                            Centroidal Momentum                               !//
    //!------------------------------------------------------------------------------!//

    //! Perform the centroidal momentum loop enabling the second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_Mu( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_Mu = robotDynamics->getDD_CentroidalMomentum();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Second-order derivative of Centroidal Momentum                             =  " << t_total/M << " microseconds" << endl;
    first = t_total/M;

    //! Perform the centroidal momentum loop enabling the second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_Mu_contracted( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_Mu_contracted = robotDynamics->getDD_ContractedCentroidalMomentum();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Second-order derivative of Centroidal Momentum when symmetry is exploited  =  " << t_total/M << " microseconds" << endl;
    second = t_total/M;

    cout << "Improvement rate                                                           =  " << 100 - (second*100)/first << " % faster" << endl << endl;


    //!------------------------------------------------------------------------------!//
    //!                              Inverse Dynamics                                !//
    //!------------------------------------------------------------------------------!//

    //! Perform the inverse dynamics loop enabling the second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_InvDyn( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_Tau = robotDynamics->getGeneralizedTorquesSecondDifferentiation();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Second-order derivative of Inverse Dynamics                             =  " << t_total/M << " microseconds" << endl;
    first = t_total/M;

    //! Perform the inverse dynamics loop enabling the second derivative
    //!------------------------------------------------------------------------------!//
    firstDerivative = true;   secondDerivative = true;
    t1 = std::chrono::high_resolution_clock::now();
    loop_InvDyn_contracted( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto DD_Tau_contracted = robotDynamics->getContractedDDTorque();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Second-order derivative of Inverse Dynamics when symmetry is exploited  =  " << t_total/M << " microseconds" << endl;
    second = t_total/M;

    cout << "Improvement rate                                                        =  " << 100 - (second*100)/first << " % faster" << endl << endl << endl;


//    geo::DirectCollocation* robot_nlp;
//    robot_nlp{nullptr};
//    robot_nlp = new geo::DirectCollocation( robot, optSettings );


return 0;
}
