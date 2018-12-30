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


#include "DcoMessage.hpp"

#include <cstring>
#include <map>

// Define grumpy messages
std::map<DISCO_Grumpy_Msg_Type, char const *> grumpyMessage;
std::map<DcoNodeBranchDir, char> grumpyDirection;

typedef struct {
    DISCO_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Dco_message;

// From CoinMessages documentation
// <3000 are informational ('I')
// <6000 warnings ('W')
// <9000 non-fatal errors ('E')
// >=9000 aborts the program (after printing the message) ('S')

// From CoinMessageHandler documentation
// log (detail) levels explained.
//
// If the log level is equal or greater than the detail level of a message, the
// message will be printed. A rough convention for the amount of output
// expected is
//
// 0 - none
// 1 - minimal
// 2 - normal low
// 3 - normal high
// 4 - verbose
//
// Please assign log levels to messages accordingly. Log levels of 8 and above
// (8,16,32, etc.) are intended for selective debugging. The logical AND of the
// log level specified in the message and the current log level is used to
// determine if the message is printed. (In other words, you're using
// individual bits to determine which messages are printed.)

// Log messages external numbers
//
// 1 DisCO welcome message
// 100-199 DcoModel::readInstance information messages
// 200-299 DcoTreeNode information messages
// 300-399 Constraint generation information messages
// 400-499 Relaxation solver information messages
// 500-599 Heuristics messages
// 600-649 General messages
// 700-800 Solution process information
// 900     Grumpy messages
//
// 6200-6299 DcoTreeNode warnings
// 6300-6399 Constraint generation warnings
// 6400-6499 Relaxation solver warnings
// 6500-6599 Heuristics warnings
// 9200-9299 DcoTreeNode error messages
// 9900-9999 general error messages
// 9300-9399 Constraint generation errors
// 9400-9499 Relaxation solver errors
// 9500-9599 Heuristics errors
// 9600-9649 Parallelization related error messages
//
//

// Debug levels, assume 32 bit integers, first 3 bits are already in use
// by log levels, we will skip 32th bit, this leaves us with 28 bits.
// We can decide 28 different debug issues, branching, cutting, etc.
//
// Bit number   use
// 0,1,2        reserved for log levels
// 3            debug branch
// 4            debug cut generation
// 5            debug node process
//
//
//

static Dco_message us_english[]=
{
    {DISCO_CUTOFF_INC, 43, 1, "Objective coefficients are multiples of %g"},
    {DISCO_CUT_STATS_FINAL, 53, 1, "Called %s cut generator %d times, generated %d cuts, used %d, CPU time %.4f seconds, current strategy %d"},
    {DISCO_CUT_STATS_NODE, 55, 1, "Node %d, called %s cut generator %d times, generated %d cuts, used %d, CPU time %.4f seconds, current strategy %d"},
    {DISCO_CUT_GENERATED, 56, DISCO_DLOG_CUT, "[%d] Cut generator %s generated %d cuts."},
    {DISCO_GAP_NO, 57, 1, "Relative optimality gap is infinity because no solution was found"},
    {DISCO_GAP_YES, 58, 1, "Relative optimality gap is %.2f%%"},
    {DISCO_ROOT_PROCESS, 30, 1, "Processing the root node (%d rows, %d columns)"},
    {DISCO_ROOT_TIME, 35, 1, "Processing the first root relaxation took CPU time %.4f seconds"},
    // solution process information
    {DISCO_NODE_LOG_HEADER, 701, 2, "Nodes Processed  Nodes Left     Lower Bound    Upper Bound    Gap(%%)  CPU Time"},
    {DISCO_NODE_LOG, 702, 2, "%-16d %-14d %s %s %s  %-d"},
    {DISCO_NODE_LOG_NO_SOL, 703, 2, "%-16d %-14d %s NA             NA      %-d"},
    // reading mps files
    {DISCO_READ_NOINTS, 20, 1, "Problem does not have integer variables"},
    {DISCO_READ_NOCONES, 21, 1, "Problem does not have conic constraints."},
    {DISCO_READ_MPSERROR, 9001, 1, "Reading conic mps file failed with code %d." },
    {DISCO_READ_MPSCBFFILEONLY,9002, 1, "Input should be in CBLIB's cbf or Mosek's conic mps format with extension cbf or mps."},
    {DISCO_READ_CONEERROR, 9003, 1, "Invalid cone type."},
    {DISCO_READ_ROTATEDCONESIZE, 9004, 1, "Rotated cones should have at least 3 members."},
    {DISCO_READ_CONESTATS1, 101, 3, "Problem has %d cones."},
    {DISCO_READ_CONESTATS2, 102, 3, "Cone %d has %d entries (type %d)"},
    {DISCO_PROBLEM_INFO, 103, 1, "Problem info.\n"
     "  Problem name: %s\n"
     "  Objective sense: %s\n"
     "  Number of variables: %d\n"
     "  Number of linear constraints: %d\n"
     "  Number of nonzero in coefficient matrix: %d\n"
     "  Number of conic constraints: %d\n"
     "  Number of integer variables: %d"},
    // tree node
    {DISCO_NODE_BRANCHONINT, 9201, 1, "[%d] Branched on integer variable. Variable index %d."},
    {DISCO_NODE_UNEXPECTEDSTATUS,9202,1, "[%d] Unexpected node status %d"},
    {DISCO_NODE_FATHOM_PARENTQ, 203, DISCO_DLOG_PROCESS,
     "[%d] Node %d fathomed due to parent quality, abs gap %f, relative gap %f."},
    {DISCO_NODE_FATHOM, 204, DISCO_DLOG_PROCESS,
     "[%d] Node %d fathomed due to quality, abs gap %f, relative gap %f."},
    {DISCO_NODE_BCP_DECISION, 205, DISCO_DLOG_PROCESS,
     "[%d] BCP function decided to keep bounding %d, branch %d, generate cons %d."},
    {DISCO_NODE_BRANCH, 206, DISCO_DLOG_BRANCH,
     "[%d] Branching node %d, variable %d, value %f, score %f."},
    {DISCO_NODE_ENCODED, 207, DISCO_DLOG_MPI, "[%d] Node %d encoded."},
    {DISCO_NODE_DECODED, 208, DISCO_DLOG_MPI, "[%d] Node decoded into %d."},
    // constraint generation
    {DISCO_INVALID_CUT_FREQUENCY,9301,1, "%d is not a valid cut frequency, changed it to %d."},
    {DISCO_INEFFECTIVE_CUT, 302, DISCO_DLOG_CUT, "[%d] Node %d, cut is ignored since the activity is low."},
    {DISCO_CUTS_ADDED, 303, DISCO_DLOG_CUT, "[%d] Node %d, %d out of %d cuts added to the solver."},
    // relaxation solver messages
    {DISCO_SOLVER_UNKNOWN_STATUS,9401, 1, "[%d] Unknown relaxation solver status."},
    {DISCO_SOLVER_FAILED,9402, 1, "[%d] Relaxation solver failed in node %d."},
    {DISCO_SOLVER_STATUS, 403, DISCO_DLOG_PROCESS,
     "[%d] Subproblem solved, status %d, obj value %f, quality %f, estimate %f."},
    {DISCO_SOLVER_INFEASIBLE, 404, DISCO_DLOG_PROCESS, "[%d] Node %d is infeasible. Status set to fathom."},
    {DISCO_FAILED_WARM_START, 9405, DISCO_DLOG_PROCESS, "[%d] Node %d, setting warm start failed."},
    {DISCO_SOLVER_ITERATIONS, 406, 1, "Total number of relaxation solver iterations %d"},
    // heuristics
    {DISCO_HEUR_BEFORE_ROOT, 501, 4, "%s heuristic found a solution; quality is %g"},
    {DISCO_HEUR_STATS_FINAL, 502, 1, "Called %s heuristic %d times, found %d solutions, CPU time %.4f seconds, current strategy %d"},
    {DISCO_HEUR_STATS_NODE, 503, DISCO_DLOG_HEURISTIC, "Node %d, called %s heuristic %d times, found %d solutions, CPU time %.4f seconds, current strategy %d"},
    {DISCO_INVALID_HEUR_FREQUENCY, 9501, 1, "%d is not a valid heuristic frequency, changed it to %d."},
    {DISCO_HEUR_SOL_FOUND, 504, DISCO_DLOG_HEURISTIC, "[%d] %s heuristic found solution, quality %f."},
    {DISCO_HEUR_NOSOL_FOUND, 505, DISCO_DLOG_HEURISTIC, "[%d] %s heuristic is called and no solution is found."},
    // branch strategies
    {DISCO_PSEUDO_REPORT, 551, DISCO_DLOG_BRANCH, "[%d] Pseudocost score of variable %d is %f."},
    {DISCO_PSEUDO_DUP, 552, DISCO_DLOG_BRANCH, "[%d] Updating down pseudocost of %d from %f to %f, frac value %f."},
    {DISCO_PSEUDO_UUP, 553, DISCO_DLOG_BRANCH, "[%d] Updating up pseudocost of %d from %f to %f, frac value %f."},
    {DISCO_STRONG_REPORT, 554, DISCO_DLOG_BRANCH, "[%d] Strong score of variable %d is %f."},

    // grumpy messages
    // time, node status, node id, parent id, branch direction, obj val [,sum
    // of column infeasibilities, count of infeasible cols]
    {DISCO_GRUMPY_MESSAGE_LONG, 900, DISCO_DLOG_GRUMPY, "[%d] %.6f %s %d %d %c %.6f %.6f %d"},
    {DISCO_GRUMPY_MESSAGE_MED, 900, DISCO_DLOG_GRUMPY, "[%d] %.6f %s %d %d %c %.6f"},
    {DISCO_GRUMPY_MESSAGE_SHORT, 900, DISCO_DLOG_GRUMPY, "[%d] %.6f %s %d %d %c"},
    // Parallelization related messages
    {DISCO_UNEXPECTED_ENCODE_STATUS, 9601, 0, "Unexpected encode return value, file: %s, line: %d."},
    {DISCO_UNEXPECTED_DECODE_STATUS, 9602, 0, "Unexpected decode return value, file: %s, line: %d."},
    // general messages
    {DISCO_INFEAS_REPORT, 603, DISCO_DLOG_PROCESS, "[%d] Column infeas %f, row infeas %f."},
    {DISCO_SOL_FOUND, 604, DISCO_DLOG_PROCESS, "[%d] Solution found, quality %f."},
    {DISCO_SOL_INT_FEAS_REPORT, 605, 1, "Integrality maximum violation %f."},
    {DISCO_SOL_CONE_FEAS_REPORT, 606, 1, "Conic constraints maximum violation %f."},
    // welcome message
    {DISCO_WELCOME, 1, 0,
     "\nThis program contains DisCO, a library for solving mixed integer second order\n"
     "cone optimization problems. Copyright 2014-2018 Lehigh University and others,\n"
     "all rights reserved. DisCO is distributed with Eclipse Public License 1.0.\n"
     "For more information visit https://github.com/coin-or/DisCO.\n"
     "Version: %s\n"
     "Build Date: %s\n"},
    {DISCO_OUT_OF_MEMORY,9901,1, "[%d] Out of memory, file: %s, line: %d."},
    {DISCO_NOT_IMPLEMENTED,9902,1, "[%d] Not implemented yet, file: %s, line: %d."},
    {DISCO_UNKNOWN_CONETYPE,9903,1, "Unknown cone type %d"},
    {DISCO_UNKNOWN_BRANCHSTRATEGY,9904,1, "Unknown branch strategy %d"},
    {DISCO_UNKNOWN_CUTSTRATEGY,9905,1, "[%d] Unknown cut strategy %d"},
    {DISCO_SHOULD_NOT_HAPPEN, 9906, 0, "[%d] Node %d, this should not happen."},
    {DISCO_DUMMY_END, 9999, 0, ""}
};

/* Constructor */
DcoMessage::DcoMessage(Language language):
  CoinMessages(sizeof(us_english) / sizeof(Dco_message)) {
  language_ = language;
  strcpy(source_, "Dco");
  Dco_message * message = us_english;
  while (message->internalNumber != DISCO_DUMMY_END) {
    CoinOneMessage oneMessage(message->externalNumber, message->detail,
                              message->message);
    addMessage(message->internalNumber, oneMessage);
    message++;
  }
  // now override any language ones
  switch (language) {
  default:
    message = NULL;
    break;
  }
  // replace if any found
  if (message) {
    while (message->internalNumber != DISCO_DUMMY_END) {
      replaceMessage(message->internalNumber, message->message);
      message++;
    }
  }
  // initialize grumpy message node status map
  grumpyMessage[DISCO_GRUMPY_BRANCHED] = "branched";
  grumpyMessage[DISCO_GRUMPY_CANDIDATE] = "candidate";
  grumpyMessage[DISCO_GRUMPY_HEURISTIC] = "heuristic";
  grumpyMessage[DISCO_GRUMPY_INTEGER] = "integer";
  grumpyMessage[DISCO_GRUMPY_FATHOMED] = "fathomed";
  grumpyMessage[DISCO_GRUMPY_PREGNANT] = "pregnant";
  grumpyMessage[DISCO_GRUMPY_INFEASIBLE] = "infeasible";
  // initialize grumpy direction map
  grumpyDirection[DcoNodeBranchDirectionDown] = 'L';
  grumpyDirection[DcoNodeBranchDirectionUp] = 'R';
}
