/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 * Ted Ralphs. All Rights Reserved.                                          *
 *===========================================================================*/


#ifndef DcoCbfIO_hpp_
#define DcoCbfIO_hpp_

class CoinPackedMatrix;

enum CONES {
  // F
  FREE_RANGE,
  // positive orthant, L+
  POSITIVE_ORT,
  // negative orthant, L-
  NEGATIVE_ORT,
  // fixpoint zero, L=
  FIXPOINT_ZERO,
  // quadratic cone, Q
  QUAD_CONE,
  // rotated quadratic cone
  RQUAD_CONE
};

/*!
  DisCO's CBF (Conic Benchmark Format) reader.

 */

class DcoCbfIO {
  /// cbf version
  int version_;
  /// direction of optimization, 1 for min, -1 for max
  int sense_;
  int num_cols_;
  int num_col_domains_;
  CONES * col_domains_;
  int * col_domain_size_;
  int num_int_;
  int * integers_;
  int num_rows_;
  int num_row_domains_;
  CONES * row_domains_;
  int * row_domain_size_;
  double * obj_coef_;
  /// coefficient matrix number of nonzeros
  int num_nz_;
  /// coefficient matrix row coordinates
  int * row_coord_;
  /// coefficient matrix column coordinates
  int * col_coord_;
  /// values of matrix coefficients
  double * coef_;
  double * fixed_term_;
public:
  DcoCbfIO();
  ~DcoCbfIO();
  void readCbf(char const * cbf_file);
  /// returns nonzero if rows are in quadratic or rotated quadratic domains
  int check_row_domains() const;
  /// get number of columns
  int getNumCols() const { return num_cols_; }
  /// get number of rows
  int getNumRows() const { return num_rows_; }
  int objSense() const { return sense_; }
  double const * objCoef() const { return obj_coef_; }
  ///@name Getting integrality
  int getNumInteger() const { return num_int_; }
  int const * integerCols() const { return integers_; }
  //@}
  //@name Getting linear constraints
  //@{
  int const * rowCoord() const { return row_coord_; }
  int const * colCoord() const { return col_coord_; }
  double const * matCoef() const { return coef_; }
  //@}
  ///@name Getting domains
  //@{
  /// get number of column domains
  int numColDomains() const { return num_col_domains_; }
  /// get type of column domains
  CONES const * colDomains() const { return col_domains_; }
  /// get number of row domains
  int numRowDomains() const { return num_row_domains_; }
  /// get type of row domains
  CONES const * rowDomains() const { return row_domains_; }
  //@}
  void getProblem(double *& colLB, double *& colUB,
                  double *& rowLB, double *& rowUB,
                  CoinPackedMatrix *& matrix,
                  int & numCones, int *& coneStart,
                  int *& coneMembers, int *& coneType) const;
  double getInfinity() const;
private:
  DcoCbfIO(DcoCbfIO const &);
  DcoCbfIO & operator=(DcoCbfIO const &);
};

#endif
