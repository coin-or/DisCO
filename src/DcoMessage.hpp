#ifndef DcoMessage_hpp_
#define DcoMessage_hpp_

#include <Dco.hpp>

//#############################################################################

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Dco messages. */
#include "CoinMessageHandler.hpp"

//#############################################################################

//todo(aykut) order these messages properly before publishing in github.
enum DISCO_Message {
    DISCO_CUTOFF_INC,
    DISCO_CUT_STATS_FINAL,
    DISCO_CUT_STATS_NODE,
    DISCO_CUT_GENERATED,
    DISCO_GAP_NO,
    DISCO_GAP_YES,
    DISCO_ROOT_PROCESS,
    DISCO_ROOT_TIME,
    // solution process information
    DISCO_NODE_LOG_HEADER,
    DISCO_NODE_LOG,
    DISCO_NODE_LOG_NO_SOL,
    // reading mps files
    DISCO_READ_NOINTS,
    DISCO_READ_NOCONES,
    DISCO_READ_MPSERROR,
    DISCO_READ_MPSFILEONLY,
    DISCO_READ_CONEERROR,
    DISCO_READ_ROTATEDCONESIZE,
    DISCO_READ_CONESTATS1,
    DISCO_READ_CONESTATS2,
    DISCO_PROBLEM_INFO,
    // tree node
    DISCO_NODE_BRANCHONINT,
    DISCO_NODE_UNEXPECTEDSTATUS,
    DISCO_NODE_FATHOM_PARENTQ,
    DISCO_NODE_FATHOM,
    DISCO_NODE_BCP_DECISION,
    DISCO_NODE_BRANCH,
    DISCO_NODE_ENCODED,
    DISCO_NODE_DECODED,
    // constraint generation
    DISCO_INVALID_CUT_FREQUENCY,
    DISCO_INEFFECTIVE_CUT,
    DISCO_CUTS_ADDED,
    // relaxation solver messages
    DISCO_SOLVER_UNKNOWN_STATUS,
    DISCO_SOLVER_FAILED,
    DISCO_SOLVER_STATUS,
    DISCO_SOLVER_INFEASIBLE,
    DISCO_FAILED_WARM_START,
    DISCO_SOLVER_ITERATIONS,
    // heuristics
    DISCO_HEUR_BEFORE_ROOT,
    DISCO_HEUR_STATS_FINAL,
    DISCO_HEUR_STATS_NODE,
    DISCO_INVALID_HEUR_FREQUENCY,
    DISCO_HEUR_SOL_FOUND,
    DISCO_HEUR_NOSOL_FOUND,
    // branch strategies
    DISCO_PSEUDO_REPORT,
    DISCO_PSEUDO_DUP,
    DISCO_PSEUDO_UUP,
    DISCO_STRONG_REPORT,
    // grumpy messages
    DISCO_GRUMPY_MESSAGE_LONG,
    DISCO_GRUMPY_MESSAGE_MED,
    DISCO_GRUMPY_MESSAGE_SHORT,
    // Parallelization related messages
    DISCO_UNEXPECTED_ENCODE_STATUS,
    DISCO_UNEXPECTED_DECODE_STATUS,
    // more general messages
    DISCO_INFEAS_REPORT,
    DISCO_SOL_INT_FEAS_REPORT,
    DISCO_SOL_CONE_FEAS_REPORT,
    DISCO_SOL_FOUND,
    // welcome message
    DISCO_WELCOME,
    // out of memory
    DISCO_OUT_OF_MEMORY,
    DISCO_NOT_IMPLEMENTED,
    DISCO_UNKNOWN_CONETYPE,
    DISCO_UNKNOWN_BRANCHSTRATEGY,
    DISCO_UNKNOWN_CUTSTRATEGY,
    DISCO_SHOULD_NOT_HAPPEN,
    ///
    DISCO_DUMMY_END
};

//todo(aykut) for now bit masking for debug level after 128
// does not work sice the bit masking is done with a char, 8 bits,
// having masking done with an int would make things wonderfull since
// DisCO has many levels of debug messages. Contact Coinutils guys
// to see how feasible is this.
enum DISCO_Debug_Level {
  DISCO_DLOG_BRANCH = 8,
  DISCO_DLOG_CUT = 16,
  DISCO_DLOG_PROCESS = 32,
  DISCO_DLOG_PRESOLVE = 64,
  DISCO_DLOG_MPI = 32,
  DISCO_DLOG_GRUMPY = 32, // 256
  DISCO_DLOG_HEURISTIC = 32 // 512

};

enum DISCO_Grumpy_Msg_Type {
  DISCO_GRUMPY_BRANCHED = 0,
  DISCO_GRUMPY_CANDIDATE,
  DISCO_GRUMPY_HEURISTIC,
  DISCO_GRUMPY_INTEGER,
  DISCO_GRUMPY_FATHOMED,
  DISCO_GRUMPY_PREGNANT,
  DISCO_GRUMPY_INFEASIBLE
};


//#############################################################################

class DcoMessage : public CoinMessages {
public:
    /**@name Constructors etc */
    //@{
    /** Constructor */
    DcoMessage(Language language=us_en);
    //@}
};

//#############################################################################
#endif
