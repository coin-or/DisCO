/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef DcoMessage_H_
#define DcoMessage_H_

//#############################################################################

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Dco messages. */
#include "CoinMessageHandler.hpp"

//#############################################################################

enum DISCO_Message
{
    DISCO_CUTOFF_INC,
    DISCO_CUT_STAT_FINAL,
    DISCO_CUT_STAT_NODE,
    DISCO_GAP_NO,
    DISCO_GAP_YES,
    DISCO_HEUR_BEFORE_ROOT,
    DISCO_HEUR_STAT_FINAL,
    DISCO_HEUR_STAT_NODE,
    DISCO_ROOT_PROCESS,
    DISCO_ROOT_TIME,
    DISCO_W_LP,
    ///
    DISCO_DUMMY_END
};

//#############################################################################

class DcoMessage : public CoinMessages 
{
public:
    
    /**@name Constructors etc */
    //@{
    /** Constructor */
    DcoMessage(Language language=us_en);
    //@}
};

//#############################################################################

#endif
