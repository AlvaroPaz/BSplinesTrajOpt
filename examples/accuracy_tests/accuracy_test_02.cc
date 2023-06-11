/**
 *	\file examples/example_02.cc
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Example to test the cost function, its gradient and its Hessian wrt control points
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
double inc_c = 1e-8;    // sensitivity

geo::VectorXr c, q, dq, ddq;
geo::MatrixXr D_q, D_dq, D_ddq, DD_q;
geo::MatrixXr numericD_CoM, numericDD_CoM, numericD_Mu, numericDD_Mu, numericD_Tau, numericDD_Tau;


int main( int argc, char** argv ){

  if (argc > 1) {
      urdf_dir.append( argv[1] );
    } else {
//      urdf_dir.append( "nao_inertial_python.urdf" );
      urdf_dir.append( "atlas_cpp_reduced.urdf" );
    }

  cout << "//!--------------------------------------------------------------------------!//" << endl;
  cout << "Unit test for accuracy in cost function and constraints derivatives wrt random c" << endl;
  cout << "//!--------------------------------------------------------------------------!//" << endl;


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


  //! Optimization parameters
  //!------------------------------------------------------------------------------!//
  std::shared_ptr< geo::robotSettingsTrajectoryOptimization > optSettings(new geo::robotSettingsTrajectoryOptimization);

  optSettings->n =  robot->getDoF();
  optSettings->numberControlPoints = 4;
  optSettings->numberPartitions = 20;
  optSettings->si = 0.0;
  optSettings->sf = 1.0;
  optSettings->S = geo::VectorXr::LinSpaced(optSettings->numberPartitions+1, optSettings->si, optSettings->sf);
  optSettings->DifferentiationWRT = geo::wrt_controlPoints;
  robot->setDifferentiationSize( optSettings->n*optSettings->numberControlPoints );

  geo::VectorXr weights( optSettings->n );
  weights.setOnes();
  weights *= 0.01;


  //! Fill up stack of constraints
  //!------------------------------------------------------------------------------!//
  optSettings->StackConstraints.clear();

  optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
  optSettings->initialConfiguration = geo::VectorXr::Zero(optSettings->n,1);

  optSettings->StackConstraints.push_back(geo::constraint_finalConfiguration);
  optSettings->finalConfiguration = geo::VectorXr::Ones(optSettings->n,1);

  optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
  optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

  optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
  optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

  optSettings->StackConstraints.push_back(geo::constraint_centerOfMass);

  optSettings->StackConstraints.push_back(geo::constraint_centroidalMomentum);


  //! Trajectory optimization instantiation
  //!------------------------------------------------------------------------------!//
  geo::DirectCollocation robotNonlinearProblem(robot, optSettings);


  //! Build basis functions
  //!------------------------------------------------------------------------------!//
  robotNonlinearProblem.buildBasisFunctions(  );


  //! Random control points
  //!------------------------------------------------------------------------------!//
  geo::VectorXr c(optSettings->n*optSettings->numberControlPoints);
  c.setRandom();


  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Objective Function Gradient wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//
  robotNonlinearProblem.computeObjectiveFunction(c, weights, true, true);

  auto cost = robotNonlinearProblem.getCost();
  auto D_cost = robotNonlinearProblem.getCostGradient();
  auto DD_cost = robotNonlinearProblem.getCostHessian();


  //! Numeric verification through finite differences
  //!------------------------------------------------------------------------------!//
  geo::MatrixXr numericD_Cost(D_cost.rows(),D_cost.cols());
  geo::real_t inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      robotNonlinearProblem.computeObjectiveFunction(c_aux1, weights, false, false);

      auto temporalCost = robotNonlinearProblem.getCost();

      numericD_Cost(iter) = (temporalCost-cost)/inc_s;
    }

  cout << endl << "//!----------------------------->> Cost Function Gradient Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_cost.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_Cost.cwiseAbs().sum() << endl;
  auto error_in_D_Cost = D_cost - numericD_Cost;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_Cost.cwiseAbs().sum() << endl << endl;


  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Objective Function Hessian wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//


  //! Numeric verification through finite differences
  geo::MatrixXr numericDD_Cost(DD_cost.rows(),DD_cost.cols());
  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      robotNonlinearProblem.computeObjectiveFunction(c_aux1, weights, true, false);

      auto temporalD_Cost = robotNonlinearProblem.getCostGradient();

      numericDD_Cost.row(iter) = (temporalD_Cost-D_cost)/inc_s;
    }

  cout << "//!----------------------------->> Cost Function Hessian Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_cost.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_Cost.cwiseAbs().sum() << endl;
  auto error_in_DD_Cost = DD_cost - numericDD_Cost;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_Cost.cwiseAbs().sum() << endl << endl;


  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Constraints Jacobian wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//


  //! Constraints Jacobian wrt Control Points Verification

  robotNonlinearProblem.computeConstraints(c, true, true);

  auto Constraints = robotNonlinearProblem.getConstraints();
  auto D_constraints = robotNonlinearProblem.getConstraintsJacobian();
  auto DD_constraints = robotNonlinearProblem.getConstraintsHessian();

  //! Numeric verification through finite differences //!
  geo::MatrixXr numericD_constraints(D_constraints.rows(),D_constraints.cols());
  inc_s = pow(2,-24);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      robotNonlinearProblem.computeConstraints(c_aux1, false, false);

      auto temporalConstraints = robotNonlinearProblem.getConstraints();

      numericD_constraints.col(iter) = (temporalConstraints-Constraints)/inc_s;
    }

  cout << "//!----------------------------->> Constraints Jacobian Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << D_constraints.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericD_constraints.cwiseAbs().sum() << endl;
  auto error_in_D_constraints = D_constraints - numericD_constraints;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_D_constraints.cwiseAbs().sum() << endl << endl;


  //! ---------------------------------------------------------------------------------------------------------- !//
  //! Constraints Hessian wrt Control Points Verification
  //! ---------------------------------------------------------------------------------------------------------- !//


  //! Numeric verification through finite differences //!
  geo::MatrixXr numericDD_constraints(DD_constraints.rows(),DD_constraints.cols());
  inc_s = 1e-8;//pow(2,-8);

  for(short int iter = 0 ; iter < c.size() ; iter++){
      //! wrt c
      auto c_aux1 = c;
      c_aux1(iter) += inc_s;

      robotNonlinearProblem.computeConstraints(c_aux1, true, false);

      auto temporalD_Constraints = robotNonlinearProblem.getConstraintsJacobian();

      numericDD_constraints.middleCols(iter*c.size(),c.size()) = (temporalD_Constraints-D_constraints)/inc_s;
    }


  cout << "//!----------------------------->> Constraints Second-Order Derivative Accuracy" << endl;
  cout << "Analytic  ==>  sum(abs( Analytic ))            =  " << DD_constraints.cwiseAbs().sum() << endl;
  cout << "Numeric   ==>  sum(abs( Numeric ))             =  " << numericDD_constraints.cwiseAbs().sum() << endl;
  auto error_in_DD_constraints = DD_constraints - numericDD_constraints;
  cout << "Error     ==>  sum(abs( Analytic - Numeric ))  =  " << error_in_DD_constraints.cwiseAbs().sum() << endl << endl;



  ////////////////////////////////////////////////////////////////////


//    //!   MODIFY THE CONSTRAINTS   !//

//    //! Fill up stack of constraints
//    optSettings->StackConstraints.clear();

//    optSettings->StackConstraints.push_back(geo::constraint_finalGeneralizedVelocity);
//    optSettings->finalGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);

//    optSettings->StackConstraints.push_back(geo::constraint_initialConfiguration);
//    optSettings->initialConfiguration = geo::VectorXr::Zero(optSettings->n,1);

//    optSettings->StackConstraints.push_back(geo::constraint_initialGeneralizedVelocity);
//    optSettings->initialGeneralizedVelocity = geo::VectorXr::Zero(optSettings->n,1);


//    //! Update last changes
//    robotNonlinearProblem.updateSettings();

return 0;
}
