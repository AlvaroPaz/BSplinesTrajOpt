// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

/**
 *	\file ipopt_interface_01.hpp
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2020
 *
 *	Header class to include the Trajectory Optimization Problem into the Ipopt Solver
 */

#ifndef IPOPT_INTERFACE_HPP
#define IPOPT_INTERFACE_HPP

#include <IpTNLP.hpp>


#include "geombd/types.h"
#include "geombd/core.h"
#include "geombd/dynamics.h"
#include "geombd/trajectoryOptimization.h"


// using namespace Ipopt;

class PracticeNLP : public Ipopt::TNLP
{
public:

    /** default constructor */
    PracticeNLP( std::shared_ptr< hr::core::MultiBody > robot, std::shared_ptr< hr::core::robotSettingsTrajectoryOptimization > robotSettings );

    /** default destructor */
    virtual ~PracticeNLP();

    // --------------------------------------------
    // Members
    // --------------------------------------------

protected:

    //! MultiBody
    std::shared_ptr< hr::core::MultiBody > robot;

    //! Multibody-Nonlinear-Problem object pointer
    hr::core::DirectCollocation* robotNonlinearProblem;

    //! Smart pointer to robot settings
    std::shared_ptr< hr::core::robotSettingsTrajectoryOptimization > robotSettings;

    //! Stack of constraints
    hr::core::ConstraintsStack StackConstraints;

    //! Number of degrees of freedom
    short int nDoF;

    //! Number of control points
    short int numberControlPoints;

    //! Number of partitions
    int numberPartitions;

    //! Number of constraints
    int numberConstraints;

    //! Time vector
    hr::VectorXr S;

    //! Initial and final times
    hr::real_t si, sf;

    //! Control points
    hr::VectorXr controlPoints;

    //! Initial configuration
    hr::VectorXr initialConfiguration;

    //! Final configuration
    hr::VectorXr finalConfiguration;

    //! Initial generalized velocity
    hr::VectorXr initialGeneralizedVelocity;

    //! Final generalized velocity
    hr::VectorXr finalGeneralizedVelocity;

    // --------------------------------------------
    // Methods
    // --------------------------------------------

public:

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Ipopt::Index& n,
                              Ipopt::Index& m,
                              Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag,
                              Ipopt::TNLP::IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Ipopt::Index n,
                                 Ipopt::Number* x_l,
                                 Ipopt::Number* x_u,
                                 Ipopt::Index m,
                                 Ipopt::Number* g_l,
                                 Ipopt::Number* g_u);

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Ipopt::Index n,
                                    bool init_x,
                                    Ipopt::Number* x,
                                    bool init_z,
                                    Ipopt::Number* z_L,
                                    Ipopt::Number* z_U,
                                    Ipopt::Index m,
                                    bool init_lambda,
                                    Ipopt::Number* lambda);

    /** Method to return the objective value */
    virtual bool eval_f(Ipopt::Index n,
                        const Ipopt::Number* x,
                        bool new_x,
                        Ipopt::Number& obj_value);

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Ipopt::Index n,
                             const Ipopt::Number* x,
                             bool new_x,
                             Ipopt::Number* grad_f);

    /** Method to return the constraint residuals */
    virtual bool eval_g(Ipopt::Index n,
                        const Ipopt::Number* x,
                        bool new_x,
                        Ipopt::Index m,
                        Ipopt::Number* g);

    /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
    virtual bool eval_jac_g(Ipopt::Index n,
                            const Ipopt::Number* x,
                            bool new_x,
                            Ipopt::Index m,
                            Ipopt::Index nele_jac,
                            Ipopt::Index* iRow,
                            Ipopt::Index *jCol,
                            Ipopt::Number* values);

    /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
    virtual bool eval_h(Ipopt::Index n,
                        const Ipopt::Number* x,
                        bool new_x,
                        Ipopt::Number obj_factor,
                        Ipopt::Index m,
                        const Ipopt::Number* lambda,
                        bool new_lambda,
                        Ipopt::Index nele_hess,
                        Ipopt::Index* iRow,
                        Ipopt::Index* jCol,
                        Ipopt::Number* values);

    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(Ipopt::SolverReturn status,
                                   Ipopt::Index n,
                                   const Ipopt::Number* x,
                                   const Ipopt::Number* z_L,
                                   const Ipopt::Number* z_U,
                                   Ipopt::Index m,
                                   const Ipopt::Number* g,
                                   const Ipopt::Number* lambda,
                                   Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq);
    //@}

private:
    /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
    //@{
    //  PracticeNLP();
    PracticeNLP(const PracticeNLP&);
    PracticeNLP& operator=(const PracticeNLP&);
    //@}
};



#endif
