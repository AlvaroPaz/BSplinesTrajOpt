#ifndef JOHNSON_H
#define JOHNSON_H

#include "geombd/core.h"


namespace geo{
namespace variational{


class Johnson : public Kinematics
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

public:

    Johnson(std::shared_ptr< MultiBody > robot, float alpha, float dt): Kinematics(robot){

        this->robot = robot;
        this->n = robot->getDoF();
        this->alpha = alpha;
        this->dt = dt;
    }

    ~Johnson(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    public:

    //! Dof
    short int n;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;

    //! Algorithm parameter
    double alpha;

    //! Time step
    double dt;

    VectorXr dL_q, dL_dq;
    MatrixXr dL_qq, dL_qdq, dL_dqdq;




    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //!
    void setDiscreteConfiguration(const VectorXr &qk, const VectorXr &qk1);

    //!
    SpatialVector dv_q(const int &n, const int &k);

    //!
    SpatialVector dv_dq(const int &n, const int &k);

    MatrixXr dG(const int &n, const int &k);

    SpatialVector dv_qi_qj(const int &n, const int &i, const int &j);

    SpatialVector dv_qj_dqi(const int &n, const int &i, const int &j);

    MatrixXr dG_qjqi(const int &n, const int &i, const int &j);

    void computeFirstLagragianDerivative();
    void computeFirstLagragianDerivative1();

    void computeSecondLagrangianDerivative();

    //! Compute the first slot Derivatives that are used to compute the dynamic equations of motion
         /*! \param configuration vectors in instant of time "k", "k+1" and "k-1"
         *   \param Slot Derivative Vector with respect to the first variable
         *   \param Slot Derivative Vector with respect to the second variable
         * \return void
         */
    void computeSlotDerivatives(const VectorXr &q0, const VectorXr &q1, const  VectorXr &q2, VectorXr &D1, VectorXr &D2);

    //! Compute the second slot Derivatives that are used to compute the derivative of  dynamic equations of motion
         /*! \param configuration vectors in instant of time "k", "k+1" and "k-1"
         * \return void
         */
    void computeSecondSlotDerivatives(const VectorXr &q0, const VectorXr &q1, const  VectorXr &q2, MatrixXr &D1D1, MatrixXr &D2D1, MatrixXr &D1D2, MatrixXr &D2D2);


};  //End Johnson class


} // End namespace variational
} // End namespace geo

#endif // JOHNSON_H

