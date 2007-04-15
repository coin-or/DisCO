/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <vector>

#include "VrpModel.h"
#include "VrpConstants.h"

//#############################################################################

/** 1) Read in the instance data
 *  2) Set colMatrix_, varLB_, varUB_, conLB_, conUB
 *     numCols_, numRows_, objSense_, objCoef_.
 *  3) Set numIntObjects_ and intColIndices_.
 */
void
VrpModel::readInstance(const char* dataFile)
{
   static char keywords[KEY_NUM][22] = {
      "NAME", 
      "NAME:",                 /* This section lists the names of the */
      "TYPE",                  /* possible fields in the data file    */
      "TYPE:",
      "COMMENT",
      "COMMENT:",
      "DIMENSION",
      "DIMENSION:",
      "CAPACITY",
      "CAPACITY:",
      "EDGE_WEIGHT_TYPE",
      "EDGE_WEIGHT_TYPE:",
      "EDGE_WEIGHT_FORMAT", 
      "EDGE_WEIGHT_FORMAT:", 
      "DISPLAY_DATA_TYPE",
      "DISPLAY_DATA_TYPE:",
      "EDGE_WEIGHT_SECTION", 
      "EDGE_WEIGHT_SECTION:", 
      "DISPLAY_DATA_SECTION", 
      "DISPLAY_DATA_SECTION:",
      "NODE_COORD_SECTION",
      "NODE_COORD_SECTION:",
      "NODE_COORD_TYPE",
      "NODE_COORD_TYPE:",
      "DEPOT_SECTION",
      "DEPOT_SECTION:",
      "CAPACITY_VOL",
      "CAPACITY_VOL:",
      "DEMAND_SECTION",
      "DEMAND_SECTION:",
      "TIME_WINDOW_SECTION",
      "TIME_WINDOW_SECTION:",
      "STANDTIME_SECTION",
      "STANDTIME_SECTION:",
      "PICKUP_SECTION",
      "PICKUP_SECTION:",
      "EOF",
      "EOF.",
      "NUMBER_OF_TRUCKS",
      "NUMBER_OF_TRUCKS:",
      "",
      "",
      "NO_MORE_TYPE"
   };
   
#define NCTYPE_NUM 3
   
   static char nctypes[NCTYPE_NUM][14] = {
      "TWOD_COORDS",
      "THREED_COORDS",     /*This section lists the possible node*/
      "NO_COORDS"          /*coordinate data types               */
   };
   
#define WTYPE_NUM 10
   
   static char wtypes[WTYPE_NUM][9] = {
      "EXPLICIT",
      "EUC_2D",            /*This is a list of the possible data types for */
      "EUC_3D",            /*edge weights                                  */
      "MAX_2D",
      "MAX_3D",
      "MAN_2D",
      "MAN_3D",
      "CEIL_2D",
      "GEO",
      "ATT"
   };
   
#define WFORMAT_NUM 9
   
   static char wformats[WFORMAT_NUM][20] = {
      "UPPER_ROW",
      "LOWER_ROW",          /*This is a list of the possible formats that*/
      "UPPER_DIAG_ROW",     /*the edge weight matrix could be given in   */
      "LOWER_DIAG_ROW",     /*if it is given explicitly                  */
      "UPPER_COL",
      "LOWER_COL",
      "UPPER_DIAG_COL",
      "LOWER_DIAG_COL",
      "FULL_MATRIX"
   };
   
#define DTYPE_NUM 3
   
   static char dtypes[DTYPE_NUM][14] = {
      "COORD_DISPLAY",
      "TWOD_DISPLAY",     /*This is a list of the various display data*/
      "NO_DISPLAY"        /*types                                     */
   };
   
   char line[LENGTH], line1[LENGTH], key[30], tmp[80];
   int wformat=-1, dtype=-1, nctype=-1;
   double fdummy;
   int i, j = 0;
   int l, m;
   FILE *f;
   int node;
   double deg, min, coord_x, coord_y, coord_z;
   double x, y;
   int capacity_vol = false;
   int k;
   
   if (!strcmp(dataFile, "")){
      printf("\nVrp I/O: No problem data file specified\n\n");
      exit(1);
   }
   
   if ((f = fopen(dataFile, "r")) == NULL){
      fprintf(stderr, "Vrp I/O: file '%s' can't be opened\n", dataFile);
      exit(1);
   }
   
   /*This loop reads in the next line of the data file and compares it
     to the list of possible keywords to determine what data will follow.
     It then reads the data into the appropriate field and iterates */
   
   while(NULL != fgets( line1, LENGTH, f)){
      strcpy(key,"");
      sscanf(line1,"%s",key); /*read in next keyword*/
      
      for (k = 0; k < KEY_NUM; k++) /*which data field comes next?*/
	 if (strcmp(keywords[k], key) == 0) break;
      
      if (k == KEY_NUM){
	 continue;
	 fprintf(stderr, "Unknown keyword! bye.\n");
	 exit(1); /*error check for acceptable data field*/
      }
      
      k >>= 1; /* This is a bit shift operation that divides k by 2    */
      /* since in the list of keywords, there are two possible*/
      /* formats for the keyword                              */
      
      if (strchr(line1,':')){
	 strcpy(line, strchr(line1, ':')+1);
      }
      
      switch (k){
	 
       case 0: /* NAME */
	 if (!sscanf(line, "%s", name_))
	    fprintf(stderr, "\nVrp I/O: error reading NAME\n\n");
	 printf("PROBLEM NAME: \t\t%s\n", name_);
	 break;
       case 1 : /*TYPE*/
	 sscanf(line, "%s", tmp);
	 if (strcmp("CVRP", tmp) != 0){
	    fprintf(stderr, "This is not a recognized file type!\n");
	    exit(1);
	 }
	 printf("TYPE: \t\t\t%s\n", tmp);
       case 2 : /*COMMENT*/
#if 0
	 if (!strncpy(tmp, line, 80))
	    fprintf(stderr, "\nVrp I/O: error reading COMMENT\n\n");
	 printf("DESCRIPTION: \t\t%s\n", tmp);
#endif
	 break;
       case 3 : /* DIMENSION */
	 if (!sscanf(line, "%i", &k)){
	    fprintf(stderr, "Vrp I/O: error reading DIMENSION\n\n");
	    exit(1);
	 }
	 vertnum_ = (int) k;
	 edgenum_ = (int) vertnum_ * (vertnum_ - 1)/2;
	 printf("DIMENSION: \t\t%i\n", k);
	 break;
       case 4 : /*CAPACITY*/
	 if (!sscanf(line, "%i", &k)){
	    fprintf(stderr, "Vrp I/O: error reading CAPACITY\n\n");
	    exit(1);
	 }
	 capacity_ = (int) k;
	 break;
       case 5 : /* EDGE_WEIGHT_TYPE */
	 sscanf(line, "%s", tmp);
	 for (wtype_ = 0; wtype_ < WTYPE_NUM; (wtype_)++)
	    if (strcmp(wtypes[wtype_], tmp) == 0) break;
	 if (wtype_ == WTYPE_NUM) {
	    fprintf(stderr, "Unknown weight type : %s !!!\n", tmp);
	    exit(1);
	 }
	 break;
       case 6 : /* EDGE_WEIGHT_FORMAT */
	 sscanf(line, "%s", tmp);
	 for (wformat = 0; wformat < WFORMAT_NUM; wformat++)
	    if (strcmp(wformats[wformat], tmp) == 0) break;
	 if (wformat == WFORMAT_NUM) {
	    fprintf(stderr, "Unknown weight type : %s !!!\n", tmp);
	    exit(1);
	 }
	 break;
       case 7 : /* DISPLAY_DATA_TYPE */
	 sscanf(line, "%s", tmp);
	 for (dtype = 0; dtype < DTYPE_NUM; dtype++)
	    if (strcmp(dtypes[dtype], tmp) == 0) break;
	 if (dtype == DTYPE_NUM) {
	    fprintf(stderr, "Unknown display type : %s !!!\n", tmp);
	    exit(1);
	 }
	 break;
       case 8: /* EDGE_WEIGHT_SECTION */
	 /*------------------------break if not EXPLICIT -*/
	 if (wtype_ != _EXPLICIT) break; 
	 switch (wformat){
	  case 1 : /* LOWER_ROW */
	  case 4 : /* UPPER_COL */
	  case 3 : /* LOWER_DIAG_ROW */
	  case 6 : /* UPPER_DIAG_COL */
	    for (i = 0; i < vertnum_; i++){
	       for (j = 0; j < i; j++){
		  if (!fscanf(f,"%lf", &fdummy)){
		     fprintf(stderr, "Not enough data -- DIMENSION or "
			     "EDGE_WEIGHT_TYPE declared wrong\n");
		     exit(1);
		  } else {
		     edges_.push_back(new VrpVariable(i, j, (int) fdummy));
		  }
	       }
	       if ((wformat==3 || wformat==6) && 
		   !fscanf(f,"%lf", &fdummy)){
		  fprintf(stderr, "Not enough data -- DIMENSION or "
			  "EDGE_WEIGHT_TYPE declared wrong\n");
		  exit(1);
	       }
	    }
	    if (fscanf(f,"%lf", &fdummy)){
	       fprintf(stderr, "Too much data -- DIMENSION or "
		       "EDGE_WEIGHT_TYPE declared wrong\n");
	       exit(1);
	    }
	    break;
	  case 0 : /* UPPER_ROW */
	  case 5 : /* LOWER_COL */
	  case 2 : /* UPPER_DIAG_ROW */
	  case 7 : /* LOWER_DIAG_COL */
	    for (i = 0; i < vertnum_; i++){
	       if (wformat==2 || wformat==7) 
		  if (!fscanf(f,"%lf", &fdummy)){
		     fprintf(stderr, "Not enough data -- DIMENSION or "
			     "EDGE_WEIGHT_TYPE declared wrong");
		     exit(1);
		  }
	       for (j= i + 1; j < vertnum_; j++){
		  if (!fscanf(f,"%lf", &fdummy)){
		     fprintf(stderr, "Not enough data -- DIMENSION or "
			     "EDGE_WEIGHT_TYPE declared wrong");
		     exit(1);
		  } else {
		     edges_[index(i, j)] = new VrpVariable(i, j, (int) fdummy);
		  }
	       }
	    }
	    if (fscanf(f,"%lf", &fdummy)){
	       fprintf(stderr, "Too much data -- DIMENSION or "
		       "EDGE_WEIGHT_TYPE declared wrong\n");
	       exit(1);
	    }
	    break;
	  case 8 : /* FULL_MATRIX */
	    for (i = 0; i < vertnum_; i++){
	       for (j = 0; j <= i; j++)
		  if(!fscanf(f,"%lf", &fdummy)){
		     fprintf(stderr, "Not enough data -- DIMENSION or "
			     "EDGE_WEIGHT_TYPE declared wrong");
		     exit(1);
		  }
	       for (j = i + 1; j < vertnum_; j++){
		  if(!fscanf(f,"%lf", &fdummy)){
		     fprintf(stderr, "Not enough data -- DIMENSION or "
			     "EDGE_WEIGHT_TYPE declared wrong");
		     exit(1);
		  }
		  edges_[index(i, j)] = new VrpVariable(i, j, (int) fdummy);
	       }
	    }
	    if (fscanf(f,"%lf", &fdummy)){
	       fprintf(stderr, "Too much data -- DIMENSION or "
		       "EDGE_WEIGHT_TYPE declared wrong\n");
	       exit(1);
	    }
	    break;
	 }
	 break;
       case 9 : /* DISPLAY_DATA_SECTION */
	 /*--------------------- break if NO_DISPLAY -*/
	 if (dtype != 1){
	    fprintf(stderr, "DISPLAY_DATA_SECTION exists"
		    "but not TWOD_DISPLAY!\n");
	    exit(1);
	 }
	 /* posx, posy -*/
	 posx_ = new int[vertnum_];
	 posy_ = new int[vertnum_];
	 for (i = 0; i < vertnum_; i++){
	    if ((k = fscanf(f,"%i%lf%lf", &node, &x, &y)) != 3){
	       fprintf(stderr, "\nVrp I/O: error reading DISPLAY_DATA\n");
	       break;
	    }
	    posx_[node-1] = (int)(x + 0.5);
	    posy_[node-1] = (int)(y + 0.5);
	 }
	 if (fscanf(f,"%lf", &fdummy)){
	    fprintf(stderr, "\nVrp I/O: too much display data\n");
	    break;
	 }
	 break;
       case 10 : /* NODE_COORD_SECTION */
	 if (nctype == -1) nctype = 0;  /*if not given: TWOD_COORDS*/
	 if (dtype == -1 && ((wtype_ == _EUC_2D) || /*display type*/
			     (wtype_ == _MAX_2D) ||	 /*not defd yet*/
			     (wtype_ == _MAN_2D)   ))/*&& can disp.*/
	    dtype = 0;				    /* COORD_DISPLAY */
	 if (dtype == 0){
	    posx_ = new int[vertnum_];
	    posy_ = new int[vertnum_];
	 }
	 coordx_ = new double[vertnum_];
	 coordy_ = new double[vertnum_];
	 if (nctype == 1)
	    coordz_ = new double[vertnum_];
	 for (i=0; i<vertnum_; i++){
	    if (nctype == 0)	     /* TWOD_COORDS */
	       if (fscanf(f,"%i%lf%lf", &node, &coord_x, &coord_y) != 3){
		  fprintf(stderr, "\nVrp I/O: error reading NODE_COORD\n\n");
		  exit(1);
	       }
	    if (nctype == 1)	     /* THREED_COORDS */
	       if (fscanf(f,"%i%lf%lf%lf", &node, &coord_x, &coord_y,
			  &coord_z) != 4){
		  fprintf(stderr, "\nVrp I/O: error reading NODE_COORD\n\n");
		  exit(1);
	       }
	    coordx_[node-1] = coord_x;
	    coordy_[node-1] = coord_y;
	    /*since position is an integer and coord is a double, I must
	      round off here if dtype is EXPLICIT*/
	    if (dtype == 0){
	       posx_[node-1] = (int)coord_x;
	       posy_[node-1] = (int)coord_y;
	    }
	    if (nctype == 1) coordz_[node-1] = coord_z;
	    if (wtype_ == _GEO){ /* GEO */
	       /*--- latitude & longitude for node ------------*/
	       deg = (int)(coordx_[node-1]);
	       min = coordx_[node-1] - deg;
	       coordx_[node-1] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	       deg = floor(coordy_[node-1]);
	       min = coordy_[node-1] - deg;
	       coordy_[node-1] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	    }
	 }
	 if (fscanf(f,"%i%lf%lf%lf", &node, &coord_x, &coord_y, &coord_z)){
	    fprintf(stderr, "\nVrp I/O: too much data in NODE_COORD\n\n");
	    exit(1);
	 }
	 break;
       case 11: /* NODE_COORD_TYPE */
	 sscanf(line, "%s", tmp);
	 for (nctype = 0; nctype < NCTYPE_NUM; nctype++)
	    if (strcmp(nctypes[nctype], tmp) == 0) break;
	 if (nctype == NCTYPE_NUM) {
	    fprintf(stderr, "Unknown node_coord_type : %s !!!\n", tmp);
	    exit(1);
	 }
	 break;
       case 12: /*DEPOT_SECTION*/
	 fscanf(f, "%i", &k);
	 if (k != 1){
	    fprintf(stderr, "Error in data: depot must be node 1");
	    exit(1);
	 }
	 depot_ = k - 1;
	 while (-1 != k) fscanf(f, "%i", &k);
	 break;
       case 13: /*CAPACITY_VOL*/
	 sscanf(line, "%i", &k);
	 capacity_vol = true;
	 break;
       case 14: /*DEMAND_SECTION*/
	 demand_ = new int[vertnum_];
	 for (i = 0; i < vertnum_; i++){
	    if (capacity_vol){
	       if (fscanf(f, "%i%i%i", &k, &l, &m) != 3){
		  fprintf(stderr,"\nVrp I/O: error reading DEMAND_SECTION\n\n");
		  exit(1);
	       }
	    }
	    else if (fscanf(f, "%i%i", &k, &l) != 2){
	       fprintf(stderr, "\nVrp I/O: error reading DEMAND_SECTION\n\n");
	       exit(1);
	    }
	    demand_[k-1] = l;
	    demand_[0] += l;
	 }
	 if (fscanf(f, "%i%i", &k, &l)){
	    fprintf(stderr, "\nVrp I/O: too much data in DEMAND_SECTION\n\n");
	    exit(1);
	 }
	 break;
       case 15: /*TIME_WINDOW_SECTION*/	/*These sections are not used*/
	 while (fscanf(f, "%d %*d:%*d %*d:%*d", &k));
	 break;
       case 16: /*STANDTIME_SECTION*/
	 while (fscanf(f, "%d%*d", &k));
	 break;
       case 17: /*PICKUP_SECTION*/	
	 while (fscanf(f, "%d%*d%*d", &k));
	 break;
       case 18: /*  EOF	*/
	 break;
       case 19: /*  NUMBER_OF_TRUCKS  */
	 if (!sscanf(line, "%i", &k)){
	    fprintf(stderr, "Vrp I/O: error reading NO_OF_TRUCKS\n\n");
	    exit(1);
	 }
	 numroutes_ = (int) k;
       default:
	 break;
      }
   }
   
   if (f != stdin)
      fclose(f);
   
   /*calculate all distances explicitly and then use distance type EXPLICIT*/
   if (wtype_ != _EXPLICIT){
      for (i = 1, k = 0; i < vertnum_; i++){
	 for (j = 0; j < i; j++){
	    edges_.push_back(new VrpVariable(i, j, compute_cost(i, j)));
	 }
      }
      wtype_ = _EXPLICIT;
   }

   setProblem(); 
}

//#############################################################################

int 
VrpModel::compute_cost(int v0, int v1){

   if (wtype_ == _EXPLICIT){
      return((int)edges_[index(v0, v1)]->getObjCoef());
   }

   double q1, q2, q3, dx, dy, dz;
   int cost = 0;
   
   if (wtype_ == _GEO){
      q1 = cos( coordy_[v0] - coordy_[v1] );
      q2 = cos( coordx_[v0] - coordx_[v1] );
      q3 = cos( coordx_[v0] + coordx_[v1] );
      cost = (int) (RRR*acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0);
   }else{
      dx = coordx_[v0] - coordx_[v1];
      dy = coordy_[v0] - coordy_[v1];
      switch (wtype_){
       case _EUC_2D : cost = (int) floor( sqrt( dx*dx + dy*dy ) + 0.5);
	 break;
       case _EUC_3D : dz = coordz_[v0] - coordz_[v1];
	 cost = (int) floor( sqrt( dx*dx + dy*dy + dz*dz) + 0.5);
	 break;
       case _MAX_2D : cost = (int) fabs(dx);
	 if (cost < fabs(dy)) cost = (int) fabs(dy);
	 break;
       case _MAX_3D : dz = coordz_[v0] - coordz_[v1];
	 cost = (int) fabs(dx);
	 if (cost < fabs(dy)) cost = (int) fabs(dy);
	 if (cost < fabs(dz)) cost = (int) fabs(dz);
	 break;
       case _MAN_2D : cost = (int) floor( dx+dy+0.5 );
	 break;
       case _MAN_3D : dz = coordz_[v0] - coordz_[v1];
	 cost = (int) floor( dx+dy+dz+0.5 );
	 break;
       case _CEIL_2D : cost = (int)ceil( sqrt( dx*dx + dy*dy ) + 0.5);
	 break;
       case _ATT     : cost = (int)( sqrt( (dx*dx + dy*dy ) / 10 ) + 1);
	 break;
      }
   }
   return(cost);
}

//#############################################################################

void 
VrpModel::setProblem()
{
   CoinBigIndex i, j, numNonzeros = 0, numInt = 0;
   int size;
   
   int* varIndices = 0;
   double* varValues = 0;

   double* collb = new double [edgenum_];
   double* colub = new double [edgenum_];
   double* conUpper = new double[vertnum_];
   double* conLower = new double[vertnum_];
   double* objCoef = new double[edgenum_];

   int* intVars = new int[numInt];

   CoinBigIndex * start = new CoinBigIndex [edgenum_+1];
   int* indices = new int[numNonzeros];
   double* values = new double[numNonzeros];

   for (i = 0; i < edgenum_; i++){
      objCoef[i] = edges_[i]->getObjCoef();
   }

   conLower[0] = conUpper[0] = numroutes_;
   for (i = 1; i < vertnum_; i++){
      conUpper[i] = conLower[i] = 2.0;
   }
   
   for (i = 0; i < edgenum_; ++i) {
      numNonzeros += edges_[i]->getSize();
      if (edges_[i]->getIntType() != 'C'){
	 numInt++;
      }
   }
   
   // Get collb, colub, obj, and matrix from variables
   for (numInt = 0, numNonzeros = 0, i = 0; i < edgenum_; ++i) {
      collb[i] = edges_[i]->getLbHard();
      colub[i] = edges_[i]->getUbHard();
      start[i] = numNonzeros;
      if (edges_[i]->getIntType() != 'C'){
	 intVars[numInt++] = i;
      }
      varValues = edges_[i]->getValues();
      varIndices = edges_[i]->getIndices();
      size = edges_[i]->getSize();
      for (j = 0; j < size; ++j, ++numNonzeros){
	 indices[numNonzeros] = varIndices[j];
	 values[numNonzeros] = varValues[j];
      }
   }
   start[i] = numNonzeros;

   // Load to lp solver
   setColMatrix(new CoinPackedMatrix(true, vertnum_, edgenum_, numNonzeros, 
				     values, indices, start, 0));
   
   setNumCons(vertnum_);
   setNumVars(edgenum_);
   setConLb(conLower);
   setConUb(conUpper);
   setVarLb(collb);
   setVarUb(colub);
   setObjSense(1.0);
   setObjCoef(objCoef);
   setInteger(intVars, numInt);    
   
   if (numInt == 0) {
      if (broker_->getMsgLevel() > 0) {
	 bcpsMessageHandler_->message(BLIS_W_LP, blisMessages())
	    << CoinMessageEol;
      }
   }

   delete [] start;
   delete [] indices;
   delete [] values;
}

//#############################################################################


