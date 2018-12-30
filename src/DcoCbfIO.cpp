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
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "DcoCbfIO.hpp"
#include "Dco.hpp"
#include <CoinPackedMatrix.hpp>
#include <CoinFinite.hpp>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <exception>


DcoCbfIO::DcoCbfIO() {
  col_domains_ = NULL;
  col_domain_size_ = NULL;
  integers_ = NULL;
  row_domains_ = NULL;
  row_domain_size_ = NULL;
  obj_coef_ = NULL;
  row_coord_ = NULL;
  col_coord_ = NULL;
  coef_ = NULL;
  fixed_term_ = NULL;
}

DcoCbfIO::DcoCbfIO(int sense,
                   int num_cols, int num_col_domains,
                   CONES const * col_domains, int const * col_domain_size,
                   int num_int, int const * integers,
                   int num_rows, int num_row_domains,
                   CONES const * row_domains, int const * row_domain_size,
                   double const * obj_coef,
                   int num_nz, int const * row_coord,
                   int const * col_coord, double const * coef,
                   double const * fixed_term) {
  version_ = 1;
  if (sense!=-1 && sense!=1) {
    std::cerr << "Sense should be either -1 (max) or 1 (min)." << std::endl;
    throw std::exception();
  }
  sense_ = sense;
  num_cols_ = num_cols;
  num_col_domains_ = num_col_domains;
  col_domains_ = new CONES[num_col_domains];
  std::copy(col_domains, col_domains+num_col_domains, col_domains_);
  col_domain_size_ = new int[num_col_domains];
  std::copy(col_domain_size, col_domain_size+num_col_domains, col_domain_size_);
  num_int_ = num_int;
  integers_ = new int[num_int];
  std::copy(integers, integers+num_int, integers_);
  num_rows_ = num_rows;
  num_row_domains_ = num_row_domains;
  row_domains_ = new CONES[num_row_domains];
  std::copy(row_domains, row_domains+num_row_domains, row_domains_);
  row_domain_size_ = new int[num_row_domains];
  std::copy(row_domain_size, row_domain_size+num_row_domains, row_domain_size_);
  obj_coef_ = new double[num_cols];
  std::copy(obj_coef, obj_coef+num_cols, obj_coef_);
  num_nz_ = num_nz;
  row_coord_ = new int[num_nz];
  std::copy(row_coord, row_coord+num_nz, row_coord_);
  col_coord_ = new int[num_nz];
  std::copy(col_coord, col_coord+num_nz, col_coord_);
  coef_ = new double[num_nz];
  std::copy(coef, coef+num_nz, coef_);
  fixed_term_ = new double[num_rows];
  std::copy(fixed_term, fixed_term+num_rows, fixed_term_);
}

void DcoCbfIO::readCbf(char const * prob_file_path) {
  // open file
  std::ifstream prob_file;
  prob_file.open(prob_file_path);
  std::string line;
  while (std::getline(prob_file, line)) {
    if (!line.compare("VER")) {
      // read VER block
      prob_file >> version_;
      if (version_ != 1) {
        std::cerr << "Only version 1 is supported." << std::endl;
        throw std::exception();
      }
    }
    else if (!line.compare("OBJSENSE")) {
      // read objective sense
      std::string sense_str;
      prob_file >> sense_str;
      sense_ = !sense_str.compare("MAX") ? -1 : 1;
    }
    else if (!line.compare("VAR")) {
      // read var
      prob_file >> num_cols_ >> num_col_domains_;
      col_domains_ = new CONES[num_col_domains_];
      col_domain_size_ = new int[num_col_domains_];
      for (int i=0; i<num_col_domains_; ++i) {
        std::string dom;
        prob_file >> dom >> col_domain_size_[i];
        if (!dom.compare("F")) {
          col_domains_[i] = FREE_RANGE;
        }
        else if (!dom.compare("L+")) {
          col_domains_[i] = POSITIVE_ORT;
        }
        else if (!dom.compare("L-")) {
          col_domains_[i] = NEGATIVE_ORT;
        }
        else if (!dom.compare("L=")) {
          col_domains_[i] = FIXPOINT_ZERO;
        }
        else if (!dom.compare("Q")) {
          col_domains_[i] = QUAD_CONE;
        }
        else if (!dom.compare("QR")) {
          col_domains_[i] = RQUAD_CONE;
        }
        else {
          std::cerr << "Unknown domain!" << std::endl;
          throw std::exception();
        }
      }
    }
    else if (!line.compare("INT")) {
      // read integrality info
      prob_file >> num_int_;
      integers_ = new int[num_int_];
      for (int i=0; i<num_int_; ++i) {
        prob_file >> integers_[i];
      }
    }
    else if (!line.compare("CON")) {
      // read constraint info
      prob_file >> num_rows_ >> num_row_domains_;
      row_domains_ = new CONES[num_row_domains_];
      row_domain_size_ = new int[num_row_domains_];
      for (int i=0; i<num_row_domains_; ++i) {
        std::string dom;
        prob_file >> dom >> row_domain_size_[i];
        if (!dom.compare("F")) {
          row_domains_[i] = FREE_RANGE;
        }
        else if (!dom.compare("L+")) {
          row_domains_[i] = POSITIVE_ORT;
        }
        else if (!dom.compare("L-")) {
          row_domains_[i] = NEGATIVE_ORT;
        }
        else if (!dom.compare("L=")) {
          row_domains_[i] = FIXPOINT_ZERO;
        }
        else if (!dom.compare("Q")) {
          row_domains_[i] = QUAD_CONE;
        }
        else if (!dom.compare("QR")) {
          row_domains_[i] = RQUAD_CONE;
        }
      }
    }
    else if (!line.compare("OBJACOORD")) {
      // read objective coef
      obj_coef_ = new double[num_cols_]();
      int num_coef;
      prob_file >> num_coef;
      for (int i=0; i<num_coef; ++i) {
        int index;
        double value;
        prob_file >> index >> value;
        obj_coef_[index] = value;
      }
    }
    else if (!line.compare("ACOORD")) {
      // read constraint coefficient
      prob_file >> num_nz_;
      row_coord_ = new int[num_nz_];
      col_coord_ = new int[num_nz_];
      coef_ = new double[num_nz_];
      for (int i=0; i<num_nz_; ++i) {
        prob_file >> row_coord_[i]
                  >> col_coord_[i]
                  >> coef_[i];
      }
    }
    else if (!line.compare("BCOORD")) {
      // read constant term
      int nonzero_rhs;
      fixed_term_ = new double[num_rows_]();
      prob_file >> nonzero_rhs;
      for (int i=0; i<nonzero_rhs; ++i) {
        int index;
        double value;
        prob_file >> index >> value;
        fixed_term_[index] = value;
      }
    }
  }
  prob_file.close();
}

void DcoCbfIO::writeCbf(std::stringstream & problem_stream) const {
  problem_stream << "# Problem written by COIN-OR DisCO." << std::endl;
  problem_stream << std::endl;
  // version
  problem_stream << "VER" << std::endl;
  problem_stream << version_ << std::endl;
  problem_stream << std::endl;
  // objective sense
  problem_stream << "OBJSENSE" << std::endl;
  if (sense_==-1) {
    problem_stream << "MAX" << std::endl;
  }
  else {
    problem_stream << "MIN" << std::endl;
  }
  problem_stream << std::endl;
  // var
  problem_stream << "VAR" << std::endl;
  problem_stream << num_cols_ << " " << num_col_domains_ << std::endl;
  for (int i=0; i<num_col_domains_; ++i) {
    std::string domain;
    if (col_domains_[i]==FREE_RANGE) {
      domain = "F";
    }
    else if (col_domains_[i]==POSITIVE_ORT) {
      domain = "L+";
    }
    else if (col_domains_[i]==NEGATIVE_ORT) {
      domain = "L-";
    }
    else if (col_domains_[i]==FIXPOINT_ZERO) {
      domain = "L=";
    }
    else if (col_domains_[i]==QUAD_CONE) {
      domain = "Q";
    }
    else if (col_domains_[i]==RQUAD_CONE) {
      domain = "QR";
    }
    problem_stream << domain << " " << col_domain_size_[i];
  }
  problem_stream << std::endl;
  problem_stream << std::endl;
  // int
  problem_stream << "INT" << std::endl;
  problem_stream << num_int_ << std::endl;
  for (int i=0; i<num_int_; ++i) {
    problem_stream << integers_[i] << std::endl;
  }
  problem_stream << std::endl;
  // con
  problem_stream << "CON" << std::endl;
  problem_stream << num_rows_ << " " << num_row_domains_ << std::endl;
  for (int i=0; i<num_row_domains_; ++i) {
    std::string domain;
    if (row_domains_[i]==FREE_RANGE) {
      domain = "F";
    }
    else if (row_domains_[i]==POSITIVE_ORT) {
      domain = "L+";
    }
    else if (row_domains_[i]==NEGATIVE_ORT) {
      domain = "L-";
    }
    else if (row_domains_[i]==FIXPOINT_ZERO) {
      domain = "L=";
    }
    else if (row_domains_[i]==QUAD_CONE) {
      domain = "Q";
    }
    else if (row_domains_[i]==RQUAD_CONE) {
      domain = "QR";
    }
    problem_stream << domain << " " << row_domain_size_[i] << std::endl;
  }
  problem_stream << std::endl;
  // objacoord
  int obj_num_nz = 0;
  for (int i=0; i<num_cols_; ++i) {
    if (obj_coef_[i]!=0.0) {
      obj_num_nz++;
    }
  }
  // find number of nonzero
  problem_stream << "OBJACOORD" << std::endl;
  problem_stream << obj_num_nz << std::endl;
  for (int i=0; i<num_cols_; ++i) {
    if (obj_coef_[i]!=0.0) {
      problem_stream << i << " " << obj_coef_[i] << std::endl;
    }
  }
  problem_stream << std::endl;
  // acoord
  problem_stream << "ACOORD" << std::endl;
  problem_stream << num_nz_ << std::endl;
  for (int i=0; i<num_nz_; ++i) {
    problem_stream << row_coord_[i] << " "
                   << col_coord_[i] << " "
                   << coef_[i]
                   << std::endl;
  }
  problem_stream << std::endl;
  // bcoord
  int ft_num_nz = 0;
  for (int i=0; i<num_rows_; ++i) {
    if (fixed_term_[i]!=0.0) {
      ft_num_nz++;
    }
  }
  problem_stream << "BCOORD" << std::endl;
  problem_stream << ft_num_nz << std::endl;
  for (int i=0; i<num_rows_; ++i) {
    if (fixed_term_[i]!=0.0) {
      problem_stream << i << " " << fixed_term_[i] << std::endl;
    }
  }
}

int DcoCbfIO::check_row_domains() const {
  for (int i=0; i<num_row_domains_; ++i) {
    if (row_domains_[i]==QUAD_CONE or
        row_domains_[i]==RQUAD_CONE) {
      return 1;
    }
  }
  return 0;
}

DcoCbfIO::~DcoCbfIO() {
  if (col_domains_) {
    delete[] col_domains_;
  }
  if (col_domain_size_) {
    delete[] col_domain_size_;
  }
  if (row_domains_) {
    delete[] row_domains_;
  }
  if (row_domain_size_) {
    delete[] row_domain_size_;
  }
  if (obj_coef_) {
    delete[] obj_coef_;
  }
  if (row_coord_) {
    delete[] row_coord_;
  }
  if (col_coord_) {
    delete[] col_coord_;
  }
  if (coef_) {
    delete[] coef_;
  }
  if (fixed_term_) {
    delete[] fixed_term_;
  }
}


// get problem in following form, i.e.,
//
// rowLB <= Ax <= rowUB
// colLB <= x  <= colUB
// x^i in L for i=1, ..., k
//
// this function is called after problem is read from file.
void DcoCbfIO::getProblem(double *& colLB, double *& colUB,
                double *& rowLB, double *& rowUB,
                CoinPackedMatrix *& matrix,
                int & numCones, int *& coneStart,
                int *& coneMembers, int *& coneType) const {
  // regular columns
  int numCols = num_cols_;
  // add columns from lifting, reduce row domains to column domains (Lorentz cones)
  // convert ax + alpha in L
  // to -y + ax + alpha = 0 and y in L
  for (int i=0; i<num_row_domains_; ++i) {
    if (row_domains_[i]==QUAD_CONE or row_domains_[i]==RQUAD_CONE) {
      numCols += row_domain_size_[i];
    }
  }
  colLB = new double[numCols];
  colUB = new double[numCols];
  std::vector<int*> cMembers;
  std::vector<int> cType;
  std::vector<int> cStart;
  cStart.push_back(0);
  // set upper and lower bounds to inf and -inf
  // no column bounds initially
  std::fill_n(colLB, numCols, -getInfinity());
  std::fill_n(colUB, numCols, getInfinity());
  // iterate over col domains and set column bounds and get cone info
  int col_index = 0;
  for (int i=0; i<num_col_domains_; ++i) {
    if (col_domains_[i] == FREE_RANGE) {
      // free domain, do nothing
    }
    else if (col_domains_[i] == POSITIVE_ORT) {
      // set lower bound of domain size many elements to 0
      std::fill_n(colLB+col_index, col_domain_size_[i], 0.0);
    }
    else if (col_domains_[i] == NEGATIVE_ORT) {
      // set upper bound of domain size many elements to 0
      std::fill_n(colUB+col_index, col_domain_size_[i], 0.0);
    }
    else if (col_domains_[i] == FIXPOINT_ZERO) {
      // set lower bound of domain size many elements to 0
      std::fill_n(colLB+col_index, col_domain_size_[i], 0.0);
      // set upper bound of domain size many elements to 0
      std::fill_n(colUB+col_index, col_domain_size_[i], 0.0);
    }
    else if (col_domains_[i] == QUAD_CONE) {
      // variables are in a lorentz cone. Add a conic constraint.
      int * members = new int[col_domain_size_[i]];
      for (int j=0; j<col_domain_size_[i]; ++j) members[j] = j+col_index;
      cMembers.push_back(members);
      cType.push_back(1);
      cStart.push_back(cStart.back()+col_domain_size_[i]);
      // leading variable is nonnegative.
      colLB[members[0]] = 0.0;
      members = 0;
    }
    else if (col_domains_[i] == RQUAD_CONE) {
      // variables are in a lorentz cone. Add a conic constraint.
      int * members = new int[col_domain_size_[i]];
      for (int j=0; j<col_domain_size_[i]; ++j) members[j] = j+col_index;
      cMembers.push_back(members);
      cType.push_back(2);
      cStart.push_back(cStart.back()+col_domain_size_[i]);
      // leading variable is nonnegative.
      colLB[members[0]] = 0.0;
      colLB[members[1]] = 0.0;
      members = 0;
    }
    col_index += col_domain_size_[i];
  }

  double * rhs = new double[num_rows_];
  for(int i=0; i<num_rows_; ++i) {
    rhs[i] = -fixed_term_[i];
  }

  // iterate over row domains and set column bounds and get cone info
  rowLB = new double[num_rows_];
  rowUB = new double[num_rows_];
  // no row bounds initially
  std::fill_n(rowLB, num_rows_, -getInfinity());
  std::fill_n(rowUB, num_rows_, getInfinity());
  // iterate over row domains and set row bounds and get cone info
  for (int i=0, row_index=0; i<num_row_domains_; ++i) {
    if (row_domains_[i]==FREE_RANGE) {
      // free domain, do nothing
    }
    else if (row_domains_[i]==POSITIVE_ORT) {
      // set lower bound of domain size many elements
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i],
                rowLB+row_index);
    }
    else if (row_domains_[i]==NEGATIVE_ORT) {
      // set upper bound of domain size many elements
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i],
                rowUB+row_index);
    }
    else if (row_domains_[i]==FIXPOINT_ZERO) {
      // set lower bound of domain size many elements
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i],
                rowLB+row_index);
      // set upper bound of domain size many elements
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i],
                rowUB+row_index);
    }
    else if (row_domains_[i]==QUAD_CONE) {
      // variables are in a lorentz cone. Add a conic constraint.
      int * members = new int[row_domain_size_[i]];
      for (int j=0; j<row_domain_size_[i]; ++j) members[j] = j+col_index;
      cMembers.push_back(members);
      cType.push_back(1);
      cStart.push_back(cStart.back()+row_domain_size_[i]);
      // leading row is nonnegative.
      colLB[members[0]] = 0.0;
      members = 0;
      col_index += row_domain_size_[i];
      // set rowlb and rowub to rhs
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i], rowLB+row_index);
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i], rowUB+row_index);
    }
    else if (row_domains_[i]==RQUAD_CONE) {
      // variables are in a lorentz cone. Add a conic constraint.
      int * members = new int[row_domain_size_[i]];
      for (int j=0; j<row_domain_size_[i]; ++j) members[j] = j+col_index;
      cMembers.push_back(members);
      cType.push_back(2);
      cStart.push_back(cStart.back()+row_domain_size_[i]);
      // leading variable is nonnegative.
      colLB[members[0]] = 0.0;
      colLB[members[1]] = 0.0;
      members = 0;
      col_index += row_domain_size_[i];
      // set rowlb and rowub to rhs
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i], rowLB+row_index);
      std::copy(rhs+row_index, rhs+row_index+row_domain_size_[i], rowUB+row_index);
    }
    row_index += row_domain_size_[i];
  }

  // compute cone info from cMembers and cType
  numCones = static_cast<int>(cMembers.size());
  coneStart = new int[numCones+1];
  coneType = new int[numCones];
  std::copy(cType.begin(), cType.end(), coneType);
  std::copy(cStart.begin(), cStart.end(), coneStart);
  int total_in_cone = 0;
  for (int i=0; i<numCones; ++i) {
    total_in_cone += coneStart[i+1] - coneStart[i];
  }
  coneMembers = new int[total_in_cone];
  for (int i=0; i<numCones; ++i) {
    int size = coneStart[i+1] - coneStart[i];
    std::copy(cMembers[i], cMembers[i]+size, coneMembers+coneStart[i]);
  }

  // get rows in CoinPackedMatrix format

  matrix = new CoinPackedMatrix(0, row_coord_, col_coord_, coef_, num_nz_);
  {
    // convert Ax + b in L
    // to -y + Ax + b = 0 and y in L
    int numcols = 0;


    // rows is always row_index, row_index+1, ...,
    // row_index+row_domain_size[i]-1
    std::vector<int> rows;

    for (int i=0, row_index=0; i<num_row_domains_; ++i) {
      if (row_domains_[i]==QUAD_CONE or row_domains_[i]==RQUAD_CONE) {
        for (int j=0; j<row_domain_size_[i]; ++j)
          rows.push_back(row_index+j);
        numcols += row_domain_size_[i];
      }
      row_index += row_domain_size_[i];
    }

    if (numcols) {
      // colstarts is always col_index, col_index+1, ...,
      // col_index+row_domain_size[i]-1
      std::vector<int> colstarts;
      for (int j=0; j<numcols; ++j)
        colstarts.push_back(j);
      colstarts.push_back(numcols);

      double * elements = new double[numcols];
      std::fill_n(elements, numcols, -1.0);
      // colstarts, rows
      int * colstarts_arr = new int[numcols+1];
      std::copy(colstarts.begin(), colstarts.end(), colstarts_arr);
      int * rows_arr = new int[numcols];
      std::copy(rows.begin(), rows.end(), rows_arr);
      matrix->appendCols(numcols, colstarts_arr, rows_arr, elements);
      delete[] rows_arr;
      delete[] colstarts_arr;
      delete[] elements;
      colstarts.clear();
    }
    rows.clear();
  }
}


double DcoCbfIO::getInfinity() const {
  return COIN_DBL_MAX;
}
