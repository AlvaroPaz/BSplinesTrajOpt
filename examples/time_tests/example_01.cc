/**
 *	\file examples/example_01.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2021
 *
 *	Example to test the dynamic objects and their two first partial derivatives wrt to control points with random c
 */

#define EIGEN_NO_DEBUG

#include "geombd/core.h"
#include "geombd/dynamics.h"

//std::string naoFile = "../../data/nao_inertial_python.urdf";


#define __FU_PATH_PREFIX__ "../../data/TROmodels/"
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
      robotDynamics->computeContractedCenterOfMass(firstDerivative, secondDerivative);
    }
}


//! Time loop for centroidal momentum
//!------------------------------------------------------------------------------!//
void loop_Mu ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeContractedCentroidalMomentum(firstDerivative, secondDerivative);
    }
}


//! Time loop for inverse dynamics
//!------------------------------------------------------------------------------!//
void loop_InvDyn ( std::shared_ptr< geo::InverseDynamics > robotDynamics ) {
  for (k = 0 ; k < M ; k++) {
      robotDynamics->computeContractedInverseDynamics(firstDerivative, secondDerivative);
    }
}


int main( int argc, char** argv ){
    if (argc > 1) {
        urdf_dir.append( argv[1] );
      } else {
        urdf_dir.append( "nao_inertial_python.urdf" );  // delete urdf once verifyed
      }

    cout << "//!----------------------------------------------------------------------!//" << endl;
    cout << "Example to test the dynamic objects time and their partial derivatives wrt c" << endl;
    cout << "//!----------------------------------------------------------------------!//" << endl;

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

    //! Perform the center of mass loop
    //!------------------------------------------------------------------------------!//
    firstDerivative = false;   secondDerivative = false;
    t1 = std::chrono::high_resolution_clock::now();
    loop_CoM( robotDynamics );
    t2 = std::chrono::high_resolution_clock::now();

    //! Time casting
    //!------------------------------------------------------------------------------!//
    auto CoM = robotDynamics->getRobotCoM();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Center of mass = " << t_total/M << endl;

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
    auto D_CoM = robotDynamics->getRobotD_CoM();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Center of mass + D = " << t_total/M << endl;

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
    auto DD_CoM = robotDynamics->getRobotDD_ContractedCoM();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Center of mass + D + DD = " << t_total/M << endl;

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
    auto Mu = robotDynamics->getCentroidalMomentum();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Centroidal momentum = " << t_total/M << endl;

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
    auto D_Mu = robotDynamics->getD_CentroidalMomentum();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Centroidal momentum + D = " << t_total/M << endl;

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
    auto DD_Mu = robotDynamics->getDD_ContractedCentroidalMomentum();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Centroidal momentum + D + DD = " << t_total/M << endl;

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
    auto Tau = robotDynamics->getGeneralizedTorques();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Inverse dynamics = " << t_total/M << endl;

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
    auto D_Tau = robotDynamics->getGeneralizedTorquesFirstDifferentiation();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Inverse dynamics + D = " << t_total/M << endl;

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
    auto DD_Tau = robotDynamics->getContractedDDTorque();
    t_total = (int) std::chrono::duration_cast< time_preci >( t2 - t1 ).count();
    cout << "Inverse dynamics + D + DD = " << t_total/M << endl;

return 0;
}
