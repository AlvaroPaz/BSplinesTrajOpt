// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Alvaro Paz

#pragma once

#ifndef VANISH_PRODUCT_H
#define VANISH_PRODUCT_H

inline int Vanisher_i;
inline int Vanisher_n;
inline int Vanisher_m;

namespace Eigen {

/*!
 * \ingroup VanishProduct_Module
 *
 * \brief The base class of dense and sparse Vanish product.
 *
 * \tparam Derived is the derived type.
 */
template<typename Derived>
class VanishProductBase : public ReturnByValue<Derived>
{
  private:
    typedef typename internal::traits<Derived> Traits;
    typedef typename Traits::Scalar Scalar;

  protected:
    typedef typename Traits::Lhs Lhs;
    typedef typename Traits::Rhs Rhs;

  public:
    /*! \brief Constructor. */
    VanishProductBase(const Lhs& A, const Rhs& B)
      : m_A(A), m_B(B)
    {}

    inline Index rows() const { return m_A.rows() * m_B.rows(); }
    inline Index cols() const { return m_A.cols() * (m_B.cols()+1) / 2; }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index row, Index col) const
    {
      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
             m_B.coeff(row % m_B.rows(), col % m_B.cols());
    }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index i) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
    }

  protected:
    typename Lhs::Nested m_A;
    typename Rhs::Nested m_B;
};

//!-------------------------------------------------------------------------

/*!
 * \ingroup VanishProduct_Module
 *
 * \brief Vanish tensor product helper class for dense matrices
 *
 * This class is the return value of VanishProduct(MatrixBase,
 * MatrixBase). Use the function rather than construct this class
 * directly to avoid specifying template prarameters.
 *
 * \tparam Lhs  Type of the left-hand side, a matrix expression.
 * \tparam Rhs  Type of the rignt-hand side, a matrix expression.
 */
template<typename Lhs, typename Rhs>
class VanishProduct : public VanishProductBase<VanishProduct<Lhs,Rhs> >
{
  private:
    typedef VanishProductBase<VanishProduct> Base;
    using Base::m_A;
    using Base::m_B;

  public:
    /*! \brief Constructor. */
    VanishProduct(const Lhs& A, const Rhs& B)
      : Base(A, B)
    {}

    /*! \brief Evaluate the Vanish tensor product. */
    template<typename Dest> void evalTo(Dest& dst) const;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
template<typename Dest>
void VanishProduct<Lhs,Rhs>::evalTo(Dest& dst) const
{
  const int BlockRows = Rhs::RowsAtCompileTime,
      BlockCols = Rhs::ColsAtCompileTime;
  const Index Br = m_B.rows(),
      Bc = m_B.cols();
  Index h, k;  h = 0;  k = 0;
  for (Index j=0; j < m_A.cols(); ++j) {
      h += k;  k = Bc-j;
      for (Index i=0; i < m_A.rows(); ++i)
        Block<Dest,BlockRows,BlockCols>(dst,i*Br,h,Br,k) = m_A.coeff(i,j) * m_B.block(0,j,Br,Bc-j);
    }
}

//!-------------------------------------------------------------------------

namespace internal {

template<typename _Lhs, typename _Rhs>
struct traits<VanishProduct<_Lhs,_Rhs> >
{
  typedef typename remove_all<_Lhs>::type Lhs;
  typedef typename remove_all<_Rhs>::type Rhs;
  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

  enum {
    Rows = size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
    MaxRows = size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
  };

  typedef Matrix<Scalar,Rows,Cols> ReturnType;
};

} // end namespace internal

//!-------------------------------------------------------------------------

/*!
 * \ingroup VanishProduct_Module
 *
 * Computes Vanish tensor product of two dense matrices
 *
 * \warning If you want to replace a matrix by its Vanish product
 *          with some matrix, do \b NOT do this:
 * \code
 * A = VanishProduct(A,B); // bug!!! caused by aliasing effect
 * \endcode
 * instead, use eval() to work around this:
 * \code
 * A = VanishProduct(A,B).eval();
 * \endcode
 *
 * \param a  Dense matrix a
 * \param b  Dense matrix b
 * \return   Vanish tensor product of a and b
 */
template<typename A, typename B>
VanishProduct<A,B> vanishProduct(const MatrixBase<A>& a, const MatrixBase<B>& b)
{
  return VanishProduct<A, B>(a.derived(), b.derived());
}


//!-------------------------------------------------------------------------
//!------------------------------NEW CLASS----------------------------------
//!-------------------------------------------------------------------------

template<typename Derived>
class VanishedKCMBase : public ReturnByValue<Derived>
{
  private:
    typedef typename internal::traits<Derived> Traits;
    typedef typename Traits::Scalar Scalar;

  protected:
    typedef typename Traits::Lhs Lhs;
    typedef typename Traits::Rhs Rhs;

  public:
    /*! \brief Constructor. */
    VanishedKCMBase(const Lhs& A, const Rhs& B)
      : m_A(A), m_B(B)
    {}

    inline Index rows() const { return m_A.rows(); }
    inline Index cols() const { return Vanisher_m*(Vanisher_m+1)/2; }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index row, Index col) const
    {
      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
             m_B.coeff(row % m_B.rows(), col % m_B.cols());
    }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index i) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
    }

  protected:
    typename Lhs::Nested m_A;
    typename Rhs::Nested m_B;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
class Vanished_KCM : public VanishedKCMBase<Vanished_KCM<Lhs,Rhs> >
{
  private:
    typedef VanishedKCMBase<Vanished_KCM> Base;
    using Base::m_A;
    using Base::m_B;

  public:
    /*! \brief Constructor. */
    Vanished_KCM(const Lhs& A, const Rhs& B)
      : Base(A, B)
    {}

    /*! \brief Evaluate the Vanish tensor product. */
    template<typename Dest> void evalTo(Dest& dst) const;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
template<typename Dest>
void Vanished_KCM<Lhs,Rhs>::evalTo(Dest& dst) const
{
  const int BlockRows = Rhs::RowsAtCompileTime,
            BlockCols = Rhs::ColsAtCompileTime;
  const Index Ar = m_A.rows(),
              Ac = m_A.cols();
  Index k = 0, z = 0, w = 0, y = 0;
  for (Index j=0; j < Vanisher_m*(Vanisher_m+1)/2; j++) {
      Block<Dest,6,1>(dst,0,j,6,1) = m_A.col(z) + m_A.col(y);
      z++;  w++;
      if(w == Vanisher_m-k){ w = 0;  k++;  z += k; }
      y += Vanisher_m;
      if(y >= Ac){ y = k + k*Vanisher_m; }
    }
}

//!-------------------------------------------------------------------------

namespace internal {

template<typename _Lhs, typename _Rhs>
struct traits<Vanished_KCM<_Lhs,_Rhs> >
{
  typedef typename remove_all<_Lhs>::type Lhs;
  typedef typename remove_all<_Rhs>::type Rhs;
  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

  enum {
    Rows = size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
    MaxRows = size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
  };

  typedef Matrix<Scalar,Rows,Cols> ReturnType;
};

} // end namespace internal

//!-------------------------------------------------------------------------

template<typename A, typename B>
Vanished_KCM<A,B> vanished_KCM(const MatrixBase<A>& a, const MatrixBase<B>& b)
{
  return Vanished_KCM<A, B>(a.derived(), b.derived());
}

//!-------------------------------------------------------------------------
//!------------------------------NEW CLASS----------------------------------
//!-------------------------------------------------------------------------


template<typename Derived>
class LeftKroneckerProductBase : public ReturnByValue<Derived>
{
  private:
    typedef typename internal::traits<Derived> Traits;
    typedef typename Traits::Scalar Scalar;

  protected:
    typedef typename Traits::Lhs Lhs;
    typedef typename Traits::Rhs Rhs;

  public:
    /*! \brief Constructor. */
    LeftKroneckerProductBase(const Lhs& A, const Rhs& B)
      : m_A(A), m_B(B)
    {}

    inline Index rows() const { return m_A.rows() * m_B.rows(); }
    inline Index cols() const { return m_A.cols() * m_B.cols(); }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index row, Index col) const
    {
      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
             m_B.coeff(row % m_B.rows(), col % m_B.cols());
    }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index i) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
    }

  protected:
    typename Lhs::Nested m_A;
    typename Rhs::Nested m_B;
};

//!-------------------------------------------------------------------------

/*!
 * \ingroup KroneckerProduct_Module
 *
 * \brief Kronecker tensor product helper class for dense matrices
 *
 * This class is the return value of kroneckerProduct(MatrixBase,
 * MatrixBase). Use the function rather than construct this class
 * directly to avoid specifying template prarameters.
 *
 * \tparam Lhs  Type of the left-hand side, a matrix expression.
 * \tparam Rhs  Type of the rignt-hand side, a matrix expression.
 */
template<typename Lhs, typename Rhs>
class LeftKroneckerProduct : public LeftKroneckerProductBase<LeftKroneckerProduct<Lhs,Rhs> >
{
  private:
    typedef LeftKroneckerProductBase<LeftKroneckerProduct> Base;
    using Base::m_A;
    using Base::m_B;

  public:
    /*! \brief Constructor. */
    LeftKroneckerProduct(const Lhs& A, const Rhs& B)
      : Base(A, B)
    {}

    /*! \brief Evaluate the Kronecker tensor product. */
    template<typename Dest> void evalTo(Dest& dst) const;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
template<typename Dest>
void LeftKroneckerProduct<Lhs,Rhs>::evalTo(Dest& dst) const
{
  const int BlockRows = Rhs::RowsAtCompileTime,
            BlockCols = Rhs::ColsAtCompileTime;
  const Index Br = m_B.rows(),
              Bc = m_B.cols();

  dst.setZero();

  for (Index i=0; i < m_A.rows(); ++i)
    for (Index j = Vanisher_i; j < m_A.cols(); j += Vanisher_n)
      Block<Dest,BlockRows,BlockCols>(dst,i*Br,j*Bc,Br,Bc) = m_A.coeff(i,j) * m_B;
}

//!-------------------------------------------------------------------------

namespace internal {

template<typename _Lhs, typename _Rhs>
struct traits<LeftKroneckerProduct<_Lhs,_Rhs> >
{
  typedef typename remove_all<_Lhs>::type Lhs;
  typedef typename remove_all<_Rhs>::type Rhs;
  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

  enum {
    Rows = size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
    MaxRows = size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
  };

  typedef Matrix<Scalar,Rows,Cols> ReturnType;
};

} // end namespace internal

//!-------------------------------------------------------------------------

template<typename A, typename B>
LeftKroneckerProduct<A,B> leftKroneckerProduct(const MatrixBase<A>& a, const MatrixBase<B>& b)
{
  return LeftKroneckerProduct<A, B>(a.derived(), b.derived());
}

} // end namespace Eigen

#endif // Vanish_TENSOR_PRODUCT_H
