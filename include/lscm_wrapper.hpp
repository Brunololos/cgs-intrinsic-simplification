// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2023 Alec Jacobson
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// Code taken from https://github.com/libigl/libigl/blob/main/include/igl/lscm.cpp
//        and https://github.com/libigl/libigl/blob/main/include/igl/lscm_hessian.cpp
// For some reason cmake always installed a version of lscm that doesn't have this signature
#ifndef LSCM_WRAPPER_H
#define LSCM_WRAPPER_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/lscm.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include <igl/cotmatrix.h>
#include <igl/vector_area_matrix.h>
#include <igl/min_quad_with_fixed.h>

template <typename DerivedV, typename DerivedF, typename QScalar>
void lscm_hessian(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    Eigen::SparseMatrix<QScalar> & Q)
{
    // Assemble the area matrix (note that A is #Vx2 by #Vx2)
    Eigen::SparseMatrix<QScalar> A;
    igl::vector_area_matrix(F,A);
    // Assemble the cotan laplacian matrix
    Eigen::SparseMatrix<QScalar> L;
    igl::cotmatrix(V,F,L);
    Eigen::SparseMatrix<QScalar> L_flat;
    igl::repdiag(L,2,L_flat);
    Q = -L_flat - 2.*A;
}

template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedV_uv>
bool lscm(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedV_uv> & V_uv)
{
    using Scalar = typename DerivedV_uv::Scalar;
    Eigen::SparseMatrix<Scalar> Q;
    lscm_hessian(V,F,Q);
    Eigen::SparseMatrix<Scalar> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    Eigen::SparseMatrix<Scalar> M2;
    igl::repdiag(M,2,M2);
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> U;
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> S;
    bool success = igl::eigs(Q,M2,3,igl::EIGS_TYPE_SM,U,S);
    if(!success)
    {
    return false;
    }
    V_uv.resize(V.rows(),2);
    V_uv<< U.col(0).head(V.rows()),U.col(0).tail(V.rows());
    return true;
}

#endif