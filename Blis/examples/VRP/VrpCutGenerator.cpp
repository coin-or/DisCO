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

#include "VrpModel.h"
#include "VrpCutGenerator.h"
#include "VrpMacros.h"
#include "VrpParams.h"

#define DELETE_POWER 3
#define DELETE_AND 0x07

/*===========================================================================*/

VrpCutGenerator::VrpCutGenerator(VrpModel *vrp, int vertnum)
{
   model_ = vrp;
   if (vertnum){
      ref_ = new int[model_->vertnum_];
      cutVal_ = new double[model_->vertnum_];
      cutList_ = new char[(model_->vertnum_ >> DELETE_POWER) + 1];
      inSet_ = new char[model_->vertnum_];
   }else{
      ref_ = 0;
      cutVal_ = 0;
      cutList_ = 0;
      inSet_ = 0;
   }
   SRANDOM(1);
   setName("VRP");
}

/*===========================================================================*/

bool 
VrpCutGenerator::generateCons(OsiCuts &cs, bool fullScan)
{
   int vertnum = model_->vertnum_;
   int comp_num = 0, rcnt, cur_bins = 0, i, k, max_node;
   double cur_slack = 0.0, node_cut, max_node_cut;
   int cut_size = (vertnum >> DELETE_POWER) + 1, num_cuts = 0;
   elist *cur_edge = NULL;
   VrpParams *par = model_->VrpPar_;
   int which_connected_routine = par->entry(VrpParams::whichConnectedRoutine);
   bool do_greedy = par->entry(VrpParams::doGreedy);
   double etol = model_->etol_;

   int *demand = model_->demand_;
   elist *cur_edge1 = NULL, *cur_edge2 = NULL;
   int node1 = 0, node2 = 0;
   int total_demand = demand[0], num_trials = 0;
   int type, rhs, capacity = model_->capacity_;
   bool found_cut = false;
   VrpNetwork *n = model_->n_;

   if (n->isIntegral_){
      /* if the network is integral, check for connectivity */
      return connectivityCuts(cs);
   }

   vertex *verts = n->verts_;

   if (which_connected_routine == BOTH) which_connected_routine = CONNECTED;
      
   int *compnodes_copy = new int[vertnum + 1];
   int *compnodes = n->compNodes_;
   double *compcuts = n->compCuts_;
   int *compdemands = n->compDemands_;

   do{
      memset(compnodes, 0, (vertnum + 1)*sizeof(int));
      memset(compcuts, 0, (vertnum + 1)*sizeof(double));
      memset(compdemands, 0, (vertnum + 1)*sizeof(int));
      
      /*------------------------------------------------------------------*\
       * Get the connected components of the solution graph without the
       * depot and see if the number of components is more than one
       \*------------------------------------------------------------------*/
      rcnt = which_connected_routine == BICONNECTED ?
	     n->biconnected() : n->connected();

      /* copy the arrays as they will be needed later */
      if (!which_connected_routine && do_greedy){
	 compnodes_copy = (int *) memcpy((char *)compnodes_copy, 
					 (char*)compnodes,
					 (vertnum + 1)*sizeof(int));
	 compnodes = compnodes_copy;
	 comp_num = rcnt;
      }
      if (rcnt > 1){
	 /*---------------------------------------------------------------*\
	  * If the number of components is more then one, then check each
	  * component to see if it violates a capacity constraint
	  \*---------------------------------------------------------------*/
	 
	 char **coef_list = new char *[rcnt];
	 memset(coef_list, 0, rcnt*sizeof(char *));
	 coef_list[0] = new char[rcnt * cut_size];
	 memset(coef_list[0], 0, rcnt*cut_size*sizeof(char));
	 for(i = 1; i < rcnt; i++)
	    coef_list[i] = coef_list[0]+i*cut_size;
	 
	 for(i = 1; i < vertnum; i++)
	    (coef_list[(verts[i].comp)-1][i >> DELETE_POWER]) |=
	       (1 << (i & DELETE_AND));
	 
	 for (i = 0; i < rcnt; i++){
	    if (compnodes[i+1] < 2) continue;
	    /*check ith component to see if it violates a constraint*/
	    if (which_connected_routine == BOTH &&
		which_connected_routine == BICONNECTED && compcuts[i+1]==0)
	       continue;
	    if (compcuts[i+1] < 2*BINS(compdemands[i+1], capacity)-etol){
	       /*the constraint is violated so impose it*/
	       type = (compnodes[i+1] < vertnum/2 ?
		       SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs = (type == SUBTOUR_ELIM_SIDE ?
		      RHS(compnodes[i+1],compdemands[i+1],
			  capacity): 2*BINS(compdemands[i+1],capacity));
	       found_cut = addCut(cs, coef_list[i], rhs, type);
	       num_cuts++;
	    }
	    else{/*if the constraint is not violated, then try generating a
		   violated constraint by deleting customers that don't
		   change the number of trucks required by the customers in
		   the component but decrease the value of the cut*/
	       cur_bins = BINS(compdemands[i+1], capacity);/*the current
							     number of trucks 
							     required*/
	       /*current slack in the constraint*/
	       cur_slack = (compcuts[i+1] - 2*cur_bins);
	       while (compnodes[i+1]){/*while there are still nodes in the
					component*/
		  for (max_node = 0, max_node_cut = 0, k = 1;
		       k < vertnum; k++){
		     if (verts[k].comp == i+1){
			if (BINS(compdemands[i+1]-verts[k].demand, capacity)
			    == cur_bins){
			   /*if the number of trucks doesn't decrease upon
			     deleting this customer*/
			   for (node_cut = 0, cur_edge = verts[k].first;
				cur_edge; cur_edge = cur_edge->next_edge){
			      node_cut += (cur_edge->other_end ?
					   -cur_edge->data->weight :
					   cur_edge->data->weight);
			   }
			   if (node_cut > max_node_cut){/*check whether the
							  value of the cut 
							  decrease is the best
							  seen so far*/
			      max_node = k;
			      max_node_cut = node_cut;
			   }
			}
		     }
		  }
		  if (!max_node){
		     break;
		  }
		  /*delete the customer that exhibited the greatest
		    decrease in cut value*/
		  compnodes[i+1]--;
		  compdemands[i+1] -= verts[max_node].demand;
		  compcuts[i+1] -= max_node_cut;
		  cur_slack -= max_node_cut;
		  verts[max_node].comp = 0;
		  coef_list[i][max_node >> DELETE_POWER] ^=
		     (1 << (max_node & DELETE_AND));
		  if (cur_slack < 0){/*if the cut is now violated, impose
				       it*/
		     type = (compnodes[i+1] < vertnum/2 ?
				      SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
		     rhs = (type == SUBTOUR_ELIM_SIDE ?
			    RHS(compnodes[i+1], compdemands[i+1],
				capacity): 2*cur_bins);
		     found_cut = addCut(cs, coef_list[i], rhs, type);
		     num_cuts++;
		     break;
		  }
	       }
	    }
	 }
	 delete [] coef_list[0];
	 delete [] coef_list;
      }
      which_connected_routine++;
   }while(!num_cuts && which_connected_routine == BOTH &&
	  which_connected_routine < 2);

   compnodes = n->compNodes_;
   
   if (!do_greedy){
      return found_cut;
   }

   if (num_cuts < 10 && do_greedy){
      int numroutes = model_->numroutes_;
      char *coef = new char[cut_size];
      for (cur_edge=verts[0].first; cur_edge; cur_edge=cur_edge->next_edge){
	 for (cur_edge1 = cur_edge->other->first; cur_edge1;
	      cur_edge1 = cur_edge1->next_edge){
	    if (cur_edge1->data->weight + cur_edge->data->weight < 1 - etol)
	       continue; 
	    node1 = cur_edge->other_end; 
	    node2 = cur_edge1->other_end;
	    for (cur_edge2 = verts[node2].first; cur_edge2;
		 cur_edge2 = cur_edge2->next_edge){
	       if (!(cur_edge2->other_end) && node2){
		  if ((BINS(total_demand - demand[node1] - demand[node2],
			    capacity) > numroutes -1) &&
		      (cur_edge1->data->weight + cur_edge->data->weight +
		       cur_edge2->data->weight>2+etol)){
		     type = SUBTOUR_ELIM_ACROSS;
		     rhs =2*BINS(total_demand - demand[node1] -
				 demand[node2],capacity);
		     memset(coef, 0, cut_size);
		     for (i = 1; i <vertnum ; i++)
			if ((i != node1) && (i != node2))
			   (coef[i >> DELETE_POWER]) |= (1 << (i&DELETE_AND));
		     found_cut = addCut(cs, coef, rhs, type);
		  }
		  break; 
	       }
	    }
	 }
      }
      delete [] coef;
   }

#if 0
   if (*num_cuts < 10 && do_greedy){
      memcpy((char *)newDemand, (char *)demand, vertnum*ISIZE);
      reduce_graph();
      if (comp_num > 1){
	 greedy_shrinking1(compnodes_copy, comp_num);
      }else{
	 greedy_shrinking1_one();
      }
   }

   if (*num_cuts < 10 && do_greedy){
      if (par->entry(VrpParams::doExtraInRoot))
	 num_trials = level ? par->entry(VrpParams::greedyNumTrials):
	    2*par->entry(VrpParams::greedyNumTrials);
      else
	 num_trials = par->entry(VrpParams::greedyNumTrials);
      if (comp_num){
	 greedy_shrinking6(compnodes_copy, comp_num, num_cuts ? num_trials :
			   2 * num_trials, 10.5);
      }else{
	 greedy_shrinking6_one(num_cuts ? num_trials : 2 * num_trials, 10.5); 
      }
   }
#endif

   delete[] compnodes_copy;

   return found_cut;
}

/*===========================================================================*/

bool 
VrpCutGenerator::connectivityCuts(OsiCuts &cs)
{
   int vertnum = model_->vertnum_;
   elist *cur_route_start;
   edge *edge_data;
   int weight = 0, reduced_weight, *route;
   int cur_vert = 0, prev_vert, cust_num = 0, cur_route, rcnt;
   int i, reduced_cust_num, vert1, vert2, type, rhs;
   int cut_size = (vertnum >> DELETE_POWER) +1;
   int capacity = model_->capacity_;
   VrpNetwork *n = model_->n_;
   vertex *verts = n->verts_;
   double *compcuts = n->compCuts_;
   int *compnodes = n->compNodes_;
   int *compdemands = n->compDemands_;
   double etol = model_->etol_;
   bool found_cut = false;

   if (!n->isIntegral_) return false;
   
   memset(compnodes, 0, (vertnum + 1)*sizeof(int));
   memset(compcuts, 0, (vertnum + 1)*sizeof(double));
   memset(compdemands, 0, (vertnum + 1)*sizeof(int));
   /*get the components of the solution graph without the depot to check if the
     graph is connected or not*/
   /* rcnt = n->connected(); */ /* This was previously executed */
   rcnt = n->numComps_;
   char **coef_list = new char *[rcnt];
   memset(coef_list, 0, rcnt*sizeof(char *));
   coef_list[0] = new char[rcnt * cut_size];
   memset(coef_list[0], 0, rcnt*cut_size*sizeof(char));
   for(i = 1; i < rcnt; i++)
      coef_list[i] = coef_list[0]+i*cut_size;
   for(i = 1; i < rcnt; i++)
     coef_list[i] = coef_list[0]+i*cut_size;

  for(i = 1; i < vertnum; i++)
     (coef_list[(verts[i].comp)-1][i >> DELETE_POWER]) |=
	(1 << (i & DELETE_AND));
  
   /*------------------------------------------------------------------------*\
   | For each component check to see if the cut it induces is nonzero -- each |
   | component's cut value must be either 0 or 2 since we have integrality    |
   \*------------------------------------------------------------------------*/
  
  for (i = 0; i < rcnt; i++){
     if (compcuts[i+1] < etol){/*if the cut value is zero, the graph is
				 disconnected and we have a violated cut*/
	type = (compnodes[i+1] < vertnum/2 ?
		SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	rhs = (type == SUBTOUR_ELIM_SIDE ? 
	       RHS(compnodes[i+1], compdemands[i+1], capacity) :
	       2*BINS(compdemands[i+1], capacity));
	found_cut = addCut(cs, coef_list[i], rhs, type);
     }
  }

  delete [] coef_list[0];
  delete [] coef_list;

  /*-------------------------------------------------------------------------*\
  | if the graph is connected, check each route to see if it obeys the        |
  | capacity constraints                                                      |
  \*-------------------------------------------------------------------------*/

  int numroutes = model_->numroutes_;
  route = new int[vertnum];
  for (cur_route_start = verts[0].first, cur_route = 0,
	  edge_data = cur_route_start->data; cur_route < numroutes;
       cur_route++){
     edge_data = cur_route_start->data;
     edge_data->scanned = true;
     cur_vert = edge_data->v1;
     prev_vert = weight = cust_num = 0;
     
     char *coef = new char[cut_size];
     memset(coef, 0, cut_size*sizeof(char));

     route[0] = cur_vert;
     while (cur_vert){
	/*keep tracing around the route and whenever the addition
	  of the next customer causes a violation, impose the
	  constraint induced
	  by the set of customers seen so far on the route*/
	coef[cur_vert >> DELETE_POWER]|=(1 << (cur_vert & DELETE_AND));
	cust_num++;
	if ((weight += verts[cur_vert].demand) > capacity){
	   type = (cust_num < vertnum/2 ?
		   SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	   rhs = (type ==SUBTOUR_ELIM_SIDE ? RHS(cust_num, weight, capacity):
		  2*BINS(weight, capacity));
	   found_cut = addCut(cs, coef, rhs, type);
	   vert1 = route[0];
	   reduced_weight = weight;
	   reduced_cust_num = cust_num;
	   while (true){
	      if ((reduced_weight -= verts[vert1].demand) > capacity){
		 reduced_cust_num--;
		 coef[vert1 >> DELETE_POWER] &=
		    ~(1 << (vert1 & DELETE_AND));
		 type = (reduced_cust_num < vertnum/2 ?
			 SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
		 rhs = (type ==SUBTOUR_ELIM_SIDE ?
			RHS(reduced_cust_num, reduced_weight, capacity):
			2*BINS(reduced_weight, capacity));
		 found_cut = addCut(cs, coef, rhs, type);
		 vert1 = route[vert1];
	      }else{
		 break;
	      }
	   }
	   vert2 = route[0];
	   while (vert2 != vert1){
	      coef[vert2 >> DELETE_POWER] |= (1 << (vert2 & DELETE_AND));
	      vert2 = route[vert2];
	   }
	}
	if (verts[cur_vert].first->other_end != prev_vert){
	   prev_vert = cur_vert;
	   edge_data = verts[cur_vert].first->data;
	   cur_vert = verts[cur_vert].first->other_end;
	}
	else{
	   prev_vert = cur_vert;
	   edge_data = verts[cur_vert].last->data; /*This statement could
						     possibly be taken out to
						     speed things up a bit*/
	   cur_vert = verts[cur_vert].last->other_end;
	}
	route[prev_vert] = cur_vert;
     }
     edge_data->scanned = true;
     
     delete [] coef;
     
     while (cur_route_start->data->scanned){/*find the next edge leading out of
					      the depot which has not yet been
					      traversed to start the next
					      route*/
	if (!(cur_route_start = cur_route_start->next_edge)) break;
     }
  }
  
  delete [] route;
  
  for (cur_route_start = verts[0].first; cur_route_start;
       cur_route_start = cur_route_start->next_edge)
     cur_route_start->data->scanned = false;

  return found_cut;
}

/*===========================================================================*/

bool
VrpCutGenerator::addCut(OsiCuts &cs, char *coef, int rhs, int type)
{
   int i, nzcnt = 0, nzcnt_side = 0, nzcnt_across = 0;
   int v0, v1; 
   std::vector<VrpVariable *>edges = model_->getEdgeList();
   int *matind = NULL, *matind_across, *matind_side;
   double *matval = NULL;
   int edgenum = model_->edgenum_;
   double infinity = model_->solver()->getInfinity();
   char sense = 'A';

   switch (type){
      /*-------------------------------------------------------------------*\
       * The subtour elimination constraints are stored as a vector of
       * bits indicating which side of the cut each customer is on
       \*-------------------------------------------------------------------*/
      
    case SUBTOUR_ELIM:
      matind_side = new int[edgenum];
      matind_across = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if (coef[v0 >> DELETE_POWER] &
	     (1 << (v0 & DELETE_AND)) &&
	     (coef[v1 >> DELETE_POWER]) &
	     (1 << (v1 & DELETE_AND))){
	    matind_side[nzcnt_side++] = i;
	 }else if ((coef[v0 >> DELETE_POWER] >>
		    (v0 & DELETE_AND) & 1) ^
		   (coef[v1 >> DELETE_POWER] >>
		    (v1 & DELETE_AND) & 1)){
	    matind_across[nzcnt_across++] = i;
	 }
      }
      type = nzcnt_side < nzcnt_across ? SUBTOUR_ELIM_SIDE :
	 SUBTOUR_ELIM_ACROSS;
      switch (type){
       case SUBTOUR_ELIM_SIDE:
	 nzcnt = nzcnt_side;
	 matind = matind_side;
	 rhs = 0; /*RHS(compnodes[i+1],compdemands[i+1], capacity)*/
	 sense = 'L';
	 delete [] matind_across;
	 break;
	 
       case SUBTOUR_ELIM_ACROSS:
	 nzcnt = nzcnt_across;
	 matind = matind_across;
	 rhs = 0; /*2*BINS(compdemands[i+1], capacity)*/
	 sense = 'G';
	 delete [] matind_side;
	 break;
      }
      
      break;
      
      
    case SUBTOUR_ELIM_SIDE:
      matind = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if (coef[v0 >> DELETE_POWER] &
	     (1 << (v0 & DELETE_AND)) &&
	     (coef[v1 >> DELETE_POWER]) &
	     (1 << (v1 & DELETE_AND))){
	    matind[nzcnt++] = i;
	 }
      }
      sense = 'L';
      break;
      
    case SUBTOUR_ELIM_ACROSS:
      matind = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if ((coef[v0 >> DELETE_POWER] >>
	      (v0 & DELETE_AND) & 1) ^
	     (coef[v1 >> DELETE_POWER] >>
	      (v1 & DELETE_AND) & 1)){
	    matind[nzcnt++] = i;
	 }
      }
      sense = 'G';
      break;
      
    default:
      break;
      
   }
   
   matval = new double[nzcnt];
   for (i = nzcnt-1; i >= 0; i--)
      matval[i] = 1;
   OsiRowCut *cut;
   if (sense == 'L'){
      cut = new OsiRowCut(infinity, rhs, edgenum, nzcnt, matind, matval);
   }else if (sense == 'G'){
      cut = new OsiRowCut(rhs, infinity, edgenum, nzcnt, matind, matval);
   }else{
      return false;
   }
      
   cs.insert(cut);
 
   return true;
}

/*===========================================================================*/


