/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/* Mark tidy up Nov 2005 */
/* General purpose header file - rf 12/90 */

#ifndef _H_general
    #define _H_general

namespace clustalw
{

    /* Macintosh specific rf 12/9/94 */
    #ifdef OS_MAC

        //#define const           /* THINK C doesn't know about these identifiers */
        //#define signed
        //#define volatile
        //#define int long
        //#define pint short            /* cast ints in printf statements as pint*/
        //#define int int            /* cast ints for sequence lengths */
        //#define lint int            /* cast ints for profile scores */

    #else

        //#define pint int            /* cast ints in printf statements as pint */
        //#define int int            /* cast ints for sequence lengths */
        //#define lint int             /* cast ints for profile scores */

    #endif

    /* definitions for all machines */

    #define EOS '\0'                /* End-Of-String */
    #define MAXLINE 5000            /* Max. line length */

}
#endif /* ifndef _H_general */


