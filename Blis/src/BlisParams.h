/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BlisParams_h
#define BlisParams_h

#include "AlpsKnowledge.h"
#include "AlpsParameterBase.h"


//#############################################################################
//#############################################################################

//** Parameters used in Blis. */
class BlisParams : public AlpsParameterSet {
 public:
  /** Character parameters. All of these variable are used as booleans
      (ture = 1, false = 0). */
  enum chrParams{
    /// Whether generate cuts at all. Default: true.
    useCons,
    /// Whether call heuristics. Default: true.
    useHeuristics,
    /// Whether generate cuts during rampup. Default: false.
    cutDuringRampup,
    //
    endOfChrParams
  };
  
  /** Integer paramters. */
  enum intParams{
      /** The number of candidate used in strong branching. Default: 10. */
      strongCandSize,

      /// -1: disable, 0: default, 1: verbose. Default: 0
      ///logLevel,

      /** Branching strategy.
          0: max infeasibilty, 
          1: pseudocost, 
          2: relibility, 
          3: strong branching.
      */
      branchStrategy,

      /** Heuristics control.
          -1: disable, 
          0: default, 
          1: moderate, 
          2: aggressive. 
      */
      heurRound,

      /** Cut generators control. 
          -1: disable, 
          0: default, 
          1: moderate, 
          2: aggressive. 
      */
      cutClique,
      cutGomory,
      cutFlowCover,
      cutKnapsack,
      cutMir,
      cutOddHole,
      cutProbing,
      cutTwoMir,
      /** The relibility of pseudocost. */
      pseudoRelibility,
      /** The look ahead of pseudocost. */
      lookAhead,
      //
      endOfIntParams
  };

  /** Double parameters. */
  enum dblParams{
    /** Tolrence to treat as an integer. Default: 1.0e-5 */
    integerTol,

    /** The value added to relaxation value when deciding fathom. 
	Default:1.0e-6 */
    cutoffInc,
    
    /** If the relative gap between best feasible and best relaxed fall into
	this gap, search stops. Default: 1.0e-6 */
    optimalRelGap,

    /** If the absolute gap between best feasible and best relaxed fall into
	this gap, search stops. Default: 1.0e-4 */
    optimalAbsGap,

    /** Weight used to calculate pseudocost. */
    pseudoWeight,

    /** Limit the max number cuts applied at a node. 
        maxNumCons = (CutFactor - 1) * numCoreConstraints. */
    cutFactor,

    /** Dense constraint factor.*/    
    denseConFactor,

    /** Scaling indicator of a constraint.*/    
    scaleConFactor,
    //
    endOfDblParams
  };

  /** String parameters. */
  enum strParams{
    strDummy,
    //
    endOfStrParams
  };

  /** There are no string array parameters. */
  enum strArrayParams{
    strArrayDummy,
    //
    endOfStrArrayParams
  };

 public:
  /**@name Constructors. */
  /*@{*/
  /** The default constructor creates a parameter set with from the template
      argument structure. The keyword list is created and the defaults are
      set. */
  BlisParams() :
    AlpsParameterSet(
		     static_cast<int>(endOfChrParams),
		     static_cast<int>(endOfIntParams),
		     static_cast<int>(endOfDblParams),
		     static_cast<int>(endOfStrParams),
		     static_cast<int>(endOfStrArrayParams)
		     )
    {
      createKeywordList();
      setDefaultEntries();
    }
  /*@}*/

  /** Method for creating the list of keyword looked for in the parameter
      file. */
  virtual void createKeywordList();
  /** Method for setting the default values for the parameters. */
  virtual void setDefaultEntries();
  /*@}*/


 public:
  //===========================================================================
  /** For user application: 
   *   Following code are do NOT need to change. 
   *   The reason can not put following functions in base class 
   *   <CODE> AlpsParameterSet </CODE> is that <CODE> chrParams </CODE>
   *   and <CODE> endOfChrParams </CODE> etc., are NOT the same as those
   *   declared in base class.
   */
  //===========================================================================

  
  /**@name Query methods 
      
     The members of the parameter set can be queried for using the overloaded
     entry() method. Using the example in the class
     documentation the user can get a parameter with the
     "<code>param.entry(USER_par::parameter_name)</code>" expression.
  */
  /*@{*/
  ///
  inline char entry(const chrParams key) const { return cpar_[key]; }
  ///
  inline int entry(const intParams key) const { return ipar_[key]; }
  ///
  inline double  entry(const dblParams key) const { return dpar_[key]; }
  ///
  inline const std::string&
    entry(const strParams key) const { return spar_[key]; }
  ///
  inline const std::vector<std::string>&
    entry(const strArrayParams key) const { return sapar_[key]; }
  /*@}*/

  //---------------------------------------------------------------------------
  /// char* is true(1) or false(0), not used
  void setEntry(const chrParams key, const char * val) {
    cpar_[key] = atoi(val); }
  /// char is true(1) or false(0), not used
  void setEntry(const chrParams key, const char val) {
    cpar_[key] = val; }
  /// This method is the one that ever been used.
  void setEntry(const chrParams key, const bool val) {
    cpar_[key] = val; }
  ///
  void setEntry(const intParams key, const char * val) {
    ipar_[key] = atoi(val); }
  ///
  void setEntry(const intParams key, const int val) {
    ipar_[key] = val; }
  ///
  void setEntry(const dblParams key, const char * val) {
    dpar_[key] = atof(val); }
  ///
  void setEntry(const dblParams key, const double val) {
    dpar_[key] = val; }
  ///
  void setEntry(const strParams key, const char * val) {
    spar_[key] = val; }
  ///
  void setEntry(const strArrayParams key, const char *val) {
    sapar_[key].push_back(val); }

  //---------------------------------------------------------------------------

  /**@name Packing/unpacking methods */
  /*@{*/
  /** Pack the parameter set into the buffer (AlpsEncoded is used 
      as buffer Here). */
  void pack(AlpsEncoded& buf) {
    buf.writeRep(cpar_, endOfChrParams)
      .writeRep(ipar_, endOfIntParams)
      .writeRep(dpar_, endOfDblParams);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.writeRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      buf.writeRep(sapar_[i].size());
      for (size_t j = 0; j < sapar_[i].size(); ++j)
	buf.writeRep(sapar_[i][j]);
    }
  }
  /** Unpack the parameter set from the buffer. */
  void unpack(AlpsEncoded& buf) {
    int dummy;
    // No need to allocate the arrays, they are of fixed length
    dummy = static_cast<int>(endOfChrParams);
    buf.readRep(cpar_, dummy, false);
    dummy = static_cast<int>(endOfIntParams);
    buf.readRep(ipar_, dummy, false);
    dummy = static_cast<int>(endOfDblParams);
    buf.readRep(dpar_, dummy, false);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.readRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      size_t str_size;
      buf.readRep(str_size);
      sapar_[i].reserve(str_size);
      for (size_t j = 0; j < str_size; ++j){
	//	sapar_[i].unchecked_push_back(std::string());
	sapar_[i].push_back(std::string());
	buf.readRep(sapar_[i].back());
      }
    }
  }
  /*@}*/

};

#endif