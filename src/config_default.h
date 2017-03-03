/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/* include the public project specific macros */
#include "config_dco_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to the debug sanity check level (0 is no test) */
#define COIN_BLIS_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_BLIS_VERBOSITY 0

/* Define to 1 if the ALPS package is used */
#define COIN_HAS_ALPS 1

/* Define to 1 if the BiCePS package is used */
#define COIN_HAS_BCPS 1

/* Define to 1 if the Blis package is used */
#define COIN_HAS_DISCO 1

/* Define to 1 if the CoinUtils package is used */
#define COIN_HAS_COINUTILS 1

/* Define to 1 if the Clp package is used */
#define COIN_HAS_CLP 1
