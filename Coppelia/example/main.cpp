//! By WhiteRose from Dark Army

#include "openhrc/core.h"
#include "openhrc/dynamics.h"
#include "openhrc/trajectoryOptimization.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <iostream>
#include <fstream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>


//! VRep headers files
extern "C" {
    #include "remoteApi/extApi.h"
    #include "remoteApi/extApiPlatform.h"
    #include "remoteApi/simConst.h"
//    #include "remoteApi/extApiInternal.h"
//    #include "remoteApi/shared_memory.h"
}


std::string naoFile = "../../OpenHRC_Geo/src/data/NAOURDF_inertial.xml";

using std::cout;
using std::endl;
using namespace Eigen;




//USING_NAMESPACE_QPOASES
#define infty 100000

std::ofstream  Qfile;

int portNb = 19997;
int clientID = simxStart("127.0.0.1", portNb, true, true, 2500, 5);

std::vector< const simxChar* > CharHandle{"HeadYaw","HeadPitch","LHipYawPitch3","LHipRoll3","LHipPitch3","LKneePitch3","LAnklePitch3","LAnkleRoll3","RHipYawPitch3","RHipRoll3","RHipPitch3","RKneePitch3","RAnklePitch3","RAnkleRoll3","LShoulderPitch3","LShoulderRoll3","LElbowYaw3","LElbowRoll3","LWristYaw3","RShoulderPitch3","RShoulderRoll3","RElbowYaw3","RElbowRoll3","RWristYaw3"};
int n = 24;

std::vector< int > IntHandle(n);
std::vector< int > FeasibleHandle(n);


int main(){

    //! Build multibody
    hr::core::World world = hr::core::World();
    std::string sFile = naoFile;
    int robotId = world.getRobotsVector()->size();
    world.loadMultiBodyURDF(sFile,robotId, hr::core::kNAO);


    std::shared_ptr<hr::core::MultiBody> robot = world.getRobot(0);

    hr::VectorXr q(robot->getDoF());    hr::VectorXr dq(robot->getDoF());    hr::VectorXr ddq(robot->getDoF());

    q.setZero();
    dq.setZero();
    ddq.setZero();

    q(15)=0.3;     q(20)=-0.3;



    for(short int ID = 0 ; ID < n ; ID++){
        FeasibleHandle.at(ID) = simxGetObjectHandle(clientID, CharHandle.at(ID), &IntHandle.at(ID), simx_opmode_blocking);
    }



    //! Get initial configuration and compute forward kinematics
    robot->setConfiguration(q);
    robot->setGeneralizedVelocity(dq);

    //! Send data init
    for(short int ID = 0 ; ID < n ; ID++){
        simxSetJointTargetPosition(clientID, IntHandle.at(ID), q(ID),  simx_opmode_blocking);
    }


    //! Initial configuration before iterating
    hr::real_t integrationStepSize = 0.01;
    hr::real_t t = 0.0;
    hr::real_t T = 5.0;
    int N = T/integrationStepSize;


    //! Dynamic variables
    Eigen::VectorXd Q=Eigen::VectorXd::Zero(30);
    Eigen::VectorXd QDot = Eigen::VectorXd::Zero(30);
    Eigen::VectorXd QDDot = Eigen::VectorXd::Zero(30);
    Eigen::VectorXd Solution;



    if (clientID == -1){
        std::cout<<std::endl << "Coppelia connection failed" <<std::endl <<std::endl;
        simxFinish(clientID);  }
    else { std::cout << std::endl << "Connected to Coppelia remote API server" << std::endl << std::endl;  }


    for(int i = 0; i != N; i++){
        cout << "here " << i << " ::: " << i*0.004 << endl;

        Q.setOnes();
        Q *= i*0.004;
        cout<<"value = "<<Q(1)<<endl;

        //! Send data
        for(short int ID = 0 ; ID < n ; ID++){
            simxSetJointTargetPosition(clientID, IntHandle.at(ID), Q(ID), simx_opmode_oneshot_wait);
        }

//        extApi_sleepMs(2000);

    }

}
