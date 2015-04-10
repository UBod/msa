/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <math.h>
#include "NJTree.h"

namespace clustalw
{

    /****************************************************************************
     * [ Improvement ideas in fast_nj_tree() ] by DDBJ & FUJITSU Limited.
     *                        written by Tadashi Koike
     *                        (takoike@genes.nig.ac.jp)
     *******************
     * <IMPROVEMENT 1> : Store the value of sum of the score to temporary array,
     *                   and use again and again.
     *
     *    In the main cycle, these are calculated again and again :
     *        diq = sum of distMat[n][ii]   (n:1 to lastSeq-firstSeq+1),
     *        djq = sum of distMat[n][jj]   (n:1 to lastSeq-firstSeq+1),
     *        dio = sum of distMat[n][mini] (n:1 to lastSeq-firstSeq+1),
     *        djq = sum of distMat[n][minj] (n:1 to lastSeq-firstSeq+1)
     *        // 'lastSeq' and 'firstSeq' are both constant values //
     *    and the result of above calculations is always same until
     *    a best pair of neighbour nodes is joined.
     *
     *    So, we change the logic to calculate the sum[i] (=sum of distMat[n][i]
     *    (n:1 to lastSeq-firstSeq+1)) and store it to array, before
     *    beginning to find a best pair of neighbour nodes, and after that
     *    we use them again and again.
     *
     *        tmat[i][j]
     *                  1   2   3   4   5
     *                +---+---+---+---+---+
     *              1 |   |   |   |   |   |
     *                +---+---+---+---+---+
     *              2 |   |   |   |   |   |  1) calculate sum of distMat[n][i]
     *                +---+---+---+---+---+        (n: 1 to lastSeq-firstSeq+1)
     *              3 |   |   |   |   |   |  2) store that sum value to sum[i]
     *                +---+---+---+---+---+
     *              4 |   |   |   |   |   |  3) use sum[i] during finding a best
     *                +---+---+---+---+---+     pair of neibour nodes.
     *              5 |   |   |   |   |   |
     *                +---+---+---+---+---+
     *                  |   |   |   |   |
     *                  V   V   V   V   V  Calculate sum , and store it to sum[i]
     *                +---+---+---+---+---+
     *         sum[i] |   |   |   |   |   |
     *                +---+---+---+---+---+
     *
     *    At this time, we thought that we use upper triangle of the matrix
     *    because distMat[i][j] is equal to distMat[j][i] and distMat[i][i] is equal
     *    to zero. Therefore, we prepared sum_rows[i] and sum_cols[i] instead
     *    of sum[i] for storing the sum value.
     *
     *        distMat[i][j]
     *                  1   2   3   4   5     sum_cols[i]
     *                +---+---+---+---+---+     +---+
     *              1     | # | # | # | # | --> |   | ... sum of distMat[1][2..5]
     *                + - +---+---+---+---+     +---+
     *              2         | # | # | # | --> |   | ... sum of distMat[2][3..5]
     *                + - + - +---+---+---+     +---+
     *              3             | # | # | --> |   | ... sum of distMat[3][4..5]
     *                + - + - + - +---+---+     +---+
     *              4                 | # | --> |   | ... sum of distMat[4][5]
     *                + - + - + - + - +---+     +---+
     *              5                     | --> |   | ... zero
     *                + - + - + - + - + - +     +---+
     *                  |   |   |   |   |
     *                  V   V   V   V   V  Calculate sum , sotre to sum[i]
     *                +---+---+---+---+---+
     *    sum_rows[i] |   |   |   |   |   |
     *                +---+---+---+---+---+
     *                  |   |   |   |   |
     *                  |   |   |   |   +----- sum of distMat[1..4][5]
     *                  |   |   |   +--------- sum of distMat[1..3][4]
     *                  |   |   +------------- sum of distMat[1..2][3]
     *                  |   +----------------- sum of distMat[1][2]
     *                  +--------------------- zero
     *
     *    And we use (sum_rows[i] + sum_cols[i]) instead of sum[i].
     *
     *******************
     * <IMPROVEMENT 2> : We manage valid nodes with chain list, instead of
     *                   tkill[i] flag array.
     *
     *    In original logic, invalid(killed?) nodes after nodes-joining
     *    are managed with tkill[i] flag array (set to 1 when killed).
     *    By this method, it is conspicuous to try next node but skip it
     *    at the latter of finding a best pair of neighbor nodes.
     *
     *    So, we thought that we managed valid nodes by using a chain list
     *    as below:
     *
     *    1) declare the list structure.
     *        struct {
     *            sint n;        // entry number of node.
     *            void *prev;        // pointer to previous entry.
     *            void *next;        // pointer to next entry.
     *        }
     *    2) construct a valid node list.
     *
     *       +-----+    +-----+    +-----+    +-----+        +-----+
     * NULL<-|prev |<---|prev |<---|prev |<---|prev |<- - - -|prev |
     *       |  0  |    |  1  |    |  2  |    |  3  |        |  n  |
     *       | next|--->| next|--->| next|--->| next|- - - ->| next|->NULL
     *       +-----+    +-----+    +-----+    +-----+        +-----+
     *
     *    3) when finding a best pair of neighbor nodes, we use
     *       this chain list as loop counter.
     *
     *    4) If an entry was killed by node-joining, this chain list is
     *       modified to remove that entry.
     *
     *       EX) remove the entry No 2.
     *       +-----+    +-----+               +-----+        +-----+
     * NULL<-|prev |<---|prev |<--------------|prev |<- - - -|prev |
     *       |  0  |    |  1  |               |  3  |        |  n  |
     *       | next|--->| next|-------------->| next|- - - ->| next|->NULL
     *       +-----+    +-----+               +-----+        +-----+
     *                             +-----+
     *                       NULL<-|prev |
     *                             |  2  |
     *                             | next|->NULL
     *                             +-----+
     *
     *    By this method, speed is up at the latter of finding a best pair of
     *    neighbor nodes.
     *
     *******************
     * <IMPROVEMENT 3> : Cut the frequency of division.
     *
     * At comparison between 'total' and 'tmin' in the main cycle, total is
     * divided by (2.0*fnseqs2) before comparison.  If N nodes are available,
     * that division happen (N*(N-1))/2 order.
     *
     * We thought that the comparison relation between tmin and total/(2.0*fnseqs2)
     * is equal to the comparison relation between (tmin*2.0*fnseqs2) and total.
     * Calculation of (tmin*2.0*fnseqs2) is only one time. so we stop dividing
     * a total value and multiply tmin and (tmin*2.0*fnseqs2) instead.
     *
     *******************
     * <IMPROVEMENT 4> : some transformation of the equation (to cut operations).
     *
     * We transform an equation of calculating 'total' in the main cycle.
     *
     */

void NJTree::generateTree(clustalw::PhyloTree* phyTree,
                          clustalw::DistMatrix* distMat,
                          clustalw::SeqInfo* seqInfo,
                          ofstream* log)
{
    if (log == NULL)
    {
        verbose = false;
    }
    
    register int i;
    int l[4], nude, k;
    int nc, mini, minj, j, ii, jj;
    double fnseqs, fnseqs2 = 0, sumd;
    double diq, djq, dij, dio, djo, da;
    double tmin, total, dmin;
    double bi, bj, b1, b2, b3, branch[4];
    int typei, typej; /* 0 = node; 1 = OTU */

    int firstSeq = seqInfo->firstSeq;
    int lastSeq = seqInfo->lastSeq;
    int numSeqs = seqInfo->numSeqs;
    
    /* IMPROVEMENT 1, STEP 0 : declare  variables */
    double *sumCols,  *sumRows,  *join;

    sumCols = new double[numSeqs + 1];
    sumRows = new double[numSeqs + 1];
    join = new double[numSeqs + 1];
    
    /* IMPROVEMENT 2, STEP 0 : declare  variables */
    int loop_limit;
    typedef struct _ValidNodeID
    {
        int n;
        struct _ValidNodeID *prev;
        struct _ValidNodeID *next;
    } ValidNodeID;
    ValidNodeID *tvalid,  *lpi,  *lpj,  *lpii,  *lpjj,  *lp_prev,  *lp_next;

    /*
     * correspondence of the loop counter variables.
     *   i .. lpi->n,    ii .. lpii->n
     *   j .. lpj->n,    jj .. lpjj->n
     */

    fnseqs = (double)lastSeq - firstSeq + 1;

    /*********************** First initialisation ***************************/

    if (verbose)
    {
        (*log) << "\n\n\t\t\tNeighbor-joining Method\n"
                << "\n Saitou, N. and Nei, M. (1987)" << " The Neighbor-joining Method:"
                << "\n A New Method for Reconstructing Phylogenetic Trees."
                << "\n Mol. Biol. Evol., 4(4), 406-425\n" << "\n\n This is an UNROOTED tree\n"
                << "\n Numbers in parentheses are branch lengths\n\n";
    }

    if (fnseqs == 2)
    {
        if (verbose)
        {
            (*log) << "Cycle   1     =  SEQ:   1 (" << setw(9) << setprecision(5) 
                    << (*distMat)(firstSeq, firstSeq + 1) 
                    << ") joins  SEQ:   2 ("
                    << setw(9) << setprecision(5) 
                    << (*distMat)(firstSeq, firstSeq + 1) << ")";
        }
        return ;
    }

    mini = minj = 0;

    /* IMPROVEMENT 1, STEP 1 : Allocate memory */
    /* IMPROVEMENT 1, STEP 2 : Initialize arrays to 0 */
    phyTree->leftBranch.resize(numSeqs + 2, 0.0);
    phyTree->rightBranch.resize(numSeqs + 2, 0.0);
    tkill.resize(numSeqs + 1, 0);
    av.resize(numSeqs + 1, 0.0);
    
    /* IMPROVEMENT 2, STEP 1 : Allocate memory */

    tvalid = new ValidNodeID[numSeqs + 1];
    
    /* tvalid[0] is special entry in array. it points a header of valid entry list */
    tvalid[0].n = 0;
    tvalid[0].prev = NULL;
    tvalid[0].next = &tvalid[1];

    /* IMPROVEMENT 2, STEP 2 : Construct and initialize the entry chain list */
    for (i = 1, loop_limit = lastSeq - firstSeq + 1, lpi = &tvalid[1],
        lp_prev = &tvalid[0], lp_next = &tvalid[2]; i <= loop_limit; ++i,
        ++lpi, ++lp_prev, ++lp_next)
    {
        (*distMat)(i, i) = av[i] = 0.0;
        tkill[i] = 0;
        lpi->n = i;
        lpi->prev = lp_prev;
        lpi->next = lp_next;
    }
    tvalid[loop_limit].next = NULL;

    /*
     * IMPROVEMENT 1, STEP 3 : Calculate the sum of score value that
     * is sequence[i] to others.
     */
    double matValue; 
    sumd = 0.0;
    for (lpj = tvalid[0].next; lpj != NULL; lpj = lpj->next)
    {
        double tmp_sum = 0.0;
        j = lpj->n;
        /* calculate sumRows[j] */
        for (lpi = tvalid[0].next; lpi->n < j; lpi = lpi->next)
        {
            i = lpi->n;
            matValue = (*distMat)(i, j);
            tmp_sum = tmp_sum + matValue;
        }
        sumRows[j] = tmp_sum;

        tmp_sum = 0.0;
        /* Set lpi to that lpi->n is greater than j */
        if ((lpi != NULL) && (lpi->n == j))
        {
            lpi = lpi->next;
        }
        /* calculate sumCols[j] */
        for (; lpi != NULL; lpi = lpi->next)
        {
            i = lpi->n;
            tmp_sum += (*distMat)(j, i);
        }
        sumCols[j] = tmp_sum;
    }

    /*********************** Enter The Main Cycle ***************************/

    for (nc = 1, loop_limit = (lastSeq - firstSeq + 1-3); nc <= loop_limit; ++nc)
    {
        sumd = 0.0;
        /* IMPROVEMENT 1, STEP 4 : use sum value */
        for (lpj = tvalid[0].next; lpj != NULL; lpj = lpj->next)
        {
            sumd += sumCols[lpj->n];
        }

        /* IMPROVEMENT 3, STEP 0 : multiply tmin and 2*fnseqs2 */
        fnseqs2 = fnseqs - 2.0; /* Set fnseqs2 at this point. */
        tmin = 99999.0 * 2.0 * fnseqs2;

        /*.................compute SMATij values and find the smallest one ........*/

        mini = minj = 0;

        /* jj must starts at least 2 */
        if ((tvalid[0].next != NULL) && (tvalid[0].next->n == 1))
        {
            lpjj = tvalid[0].next->next;
        }
        else
        {
            lpjj = tvalid[0].next;
        }

        for (; lpjj != NULL; lpjj = lpjj->next)
        {
            jj = lpjj->n;
            for (lpii = tvalid[0].next; lpii->n < jj; lpii = lpii->next)
            {
                ii = lpii->n;
                diq = djq = 0.0;

                /* IMPROVEMENT 1, STEP 4 : use sum value */
                diq = sumCols[ii] + sumRows[ii];
                djq = sumCols[jj] + sumRows[jj];
                /*
                 * always ii < jj in this point. Use upper
                 * triangle of score matrix.
                 */
                dij = (*distMat)(ii, jj);
                /*
                 * IMPROVEMENT 3, STEP 1 : fnseqs2 is
                 * already calculated.
                 */
                /* fnseqs2 = fnseqs - 2.0 */

                /* IMPROVEMENT 4 : transform the equation */
                /*-------------------------------------------------------------------*
                 * OPTIMIZE of expression 'total = d2r + fnseqs2*dij + dr*2.0'       *
                 * total = d2r + fnseq2*dij + 2.0*dr                                 *
                 *       = d2r + fnseq2*dij + 2(sumd - dij - d2r)                    *
                 *       = d2r + fnseq2*dij + 2*sumd - 2*dij - 2*d2r                 *
                 *       =       fnseq2*dij + 2*sumd - 2*dij - 2*d2r + d2r           *
                 *       = fnseq2*dij + 2*sumd - 2*dij - d2r                         *
                 *       = fnseq2*dij + 2*sumd - 2*dij - (diq + djq - 2*dij)         *
                 *       = fnseq2*dij + 2*sumd - 2*dij - diq - djq + 2*dij           *
                 *       = fnseq2*dij + 2*sumd - 2*dij + 2*dij - diq - djq           *
                 *       = fnseq2*dij + 2*sumd  - diq - djq                          *
                 *-------------------------------------------------------------------*/
                total = fnseqs2 * dij + 2.0 * sumd - diq - djq;
                /*
                 * IMPROVEMENT 3, STEP 2 : abbrevlate
                 * the division on comparison between
                 * total and tmin.
                 */
                /* total = total / (2.0*fnseqs2); */

                if (total < tmin)
                {
                    tmin = total;
                    mini = ii;
                    minj = jj;
                }
            }
        }

        /* MEMO: always ii < jj in avobe loop, so mini < minj */

        /*.................compute branch lengths and print the results ........*/


        dio = djo = 0.0;

        /* IMPROVEMENT 1, STEP 4 : use sum value */
        dio = sumCols[mini] + sumRows[mini];
        djo = sumCols[minj] + sumRows[minj];

        dmin = (*distMat)(mini, minj);
        dio = (dio - dmin) / fnseqs2;
        djo = (djo - dmin) / fnseqs2;
        bi = (dmin + dio - djo) *0.5;
        bj = dmin - bi;
        bi = bi - av[mini];
        bj = bj - av[minj];
        
#if 0
        (*log) << endl << "***  cycle " << nc << endl;
        (*log) << "dmin     = " << setw(9) << right << setprecision(9) << dmin << endl;
        (*log) << "dio      = " << setw(9) << right << setprecision(9) <<  dio << endl;
        (*log) << "djo      = " << setw(9) << right << setprecision(9) <<  djo << endl;
        (*log) << "bi       = " << setw(9) << right << setprecision(9) << bi << endl;
        (*log) << "bj       = " << setw(9) << right << setprecision(9) <<  bj << endl;
        (*log) << "mini     = " << setw(9) << mini << endl;
        (*log) << "minj     = " << setw(9) << minj << endl;
        (*log) << "av[minj] = " << setw(9) << right << setprecision(9) <<  av[minj] << endl;
        (*log) << "av[minj] = " << setw(9) << right << setprecision(9) << av[minj] << endl;
#endif
        if (av[mini] > 0.0)
        {
            typei = 0;
        }
        else
        {
            typei = 1;
        }
        if (av[minj] > 0.0)
        {
            typej = 0;
        }
        else
        {
            typej = 1;
        }

        if (verbose)
        {
            (*log) <<  "\n Cycle" << setw(4) << nc << "     = ";
        }

        /*
         set (tiny? (AW&FS)) negative branch lengths to zero.  Also set any tiny positive
         branch lengths to zero.
         */
        if (fabs(bi) < 0.0001)
        {
            bi = 0.0;
        }
        if (fabs(bj) < 0.0001)
        {
            bj = 0.0;
        }

        if (verbose)
        {
            if (typei == 0)
            {
                (*log) <<  "Node:" << setw(4) << mini << " (" << setw(9) << setprecision(5) 
                        << bi << ") joins ";
            }
            else
            {
                (*log) << " SEQ:" << setw(4) << mini << " (" << setw(9) << setprecision(5)
                        << bi << ") joins ";
            }

            if (typej == 0)
            {
                (*log) << "Node:" << setw(4) << minj << " (" << setw(9) << setprecision(5)
                        << bj << ")";
            }
            else
            {
                (*log) << " SEQ:" << setw(4) << minj << " (" << setw(9) << setprecision(5)
                        << bj << ")";
            }

            (*log) << "\n";
        }

        phyTree->leftBranch[nc] = bi;
        phyTree->rightBranch[nc] = bj;

        for (i = 1; i <= lastSeq - firstSeq + 1; i++)
        {
            phyTree->treeDesc[nc][i] = 0;
        }

        if (typei == 0)
        {
            for (i = nc - 1; i >= 1; i--)
                if (phyTree->treeDesc[i][mini] == 1)
                {
                    for (j = 1; j <= lastSeq - firstSeq + 1; j++)
                        if (phyTree->treeDesc[i][j] == 1)
                        {
                            phyTree->treeDesc[nc][j] = 1;
                        }
                    break;
                }
        }
        else
        {
            phyTree->treeDesc[nc][mini] = 1;
        }

        if (typej == 0)
        {
            for (i = nc - 1; i >= 1; i--)
                if (phyTree->treeDesc[i][minj] == 1)
                {
                    for (j = 1; j <= lastSeq - firstSeq + 1; j++)
                        if (phyTree->treeDesc[i][j] == 1)
                        {
                            phyTree->treeDesc[nc][j] = 1;
                        }
                    break;
                }
        }
        else
        {
            phyTree->treeDesc[nc][minj] = 1;
        }


        /*
        Here is where the -0.00005 branch lengths come from for 3 or more
        identical seqs.
         */
        /*        if(dmin <= 0.0) dmin = 0.0001; */
        if (dmin <= 0.0)
        {
            dmin = 0.000001;
        }
        av[mini] = dmin * 0.5;

        /*........................Re-initialisation................................*/

        fnseqs = fnseqs - 1.0;
        tkill[minj] = 1;

        /* IMPROVEMENT 2, STEP 3 : Remove tvalid[minj] from chain list. */
        /* [ Before ]
         *  +---------+        +---------+        +---------+
         *  |prev     |<-------|prev     |<-------|prev     |<---
         *  |    n    |        | n(=minj)|        |    n    |
         *  |     next|------->|     next|------->|     next|----
         *  +---------+        +---------+        +---------+
         *
         * [ After ]
         *  +---------+                           +---------+
         *  |prev     |<--------------------------|prev     |<---
         *  |    n    |                           |    n    |
         *  |     next|-------------------------->|     next|----
         *  +---------+                           +---------+
         *                     +---------+
         *              NULL---|prev     |
         *                     | n(=minj)|
         *                     |     next|---NULL
         *                     +---------+
         */
        (tvalid[minj].prev)->next = tvalid[minj].next;
        if (tvalid[minj].next != NULL)
        {
            (tvalid[minj].next)->prev = tvalid[minj].prev;
        }
        tvalid[minj].prev = tvalid[minj].next = NULL;

        /* IMPROVEMENT 1, STEP 5 : re-calculate sum values. */
        for (lpj = tvalid[0].next; lpj != NULL; lpj = lpj->next)
        {
            double tmp_di = 0.0;
            double tmp_dj = 0.0;
            j = lpj->n;

            /*
             * subtrace a score value related with 'minj' from
             * sum arrays .
             */
            if (j < minj)
            {
                tmp_dj = (*distMat)(j, minj);
                sumCols[j] -= tmp_dj;
            }
            else if (j > minj)
            {
                tmp_dj = (*distMat)(minj, j);
                sumRows[j] -= tmp_dj;
            } /* nothing to do when j is equal to minj. */


            /*
             * subtrace a score value related with 'mini' from
             * sum arrays .
             */
            if (j < mini)
            {
                tmp_di = (*distMat)(j, mini);
                sumCols[j] -= tmp_di;
            }
            else if (j > mini)
            {
                tmp_di = (*distMat)(mini, j);
                sumRows[j] -= tmp_di;
            } /* nothing to do when j is equal to mini. */

            /*
             * calculate a score value of the new inner node.
             * then, store it temporary to join[] array.
             */
            join[j] = (tmp_dj + tmp_di) *0.5;
        }

        /*
         * 1)
         * Set the score values (stored in join[]) into the matrix,
         * row/column position is 'mini'.
         * 2)
         * Add a score value of the new inner node to sum arrays.
         */
        for (lpj = tvalid[0].next; lpj != NULL; lpj = lpj->next)
        {
            j = lpj->n;
            if (j < mini)
            {
                distMat->SetAt(j, mini, join[j]);
                sumCols[j] += join[j];
            }
            else if (j > mini)
            {
                distMat->SetAt(mini, j, join[j]);
                sumRows[j] += join[j];
            } /* nothing to do when j is equal to mini. */
        }

        /* Re-calculate sumRows[mini],sumCols[mini]. */
        sumCols[mini] = sumRows[mini] = 0.0;

        /* calculate sumRows[mini] */
        da = 0.0;
        for (lpj = tvalid[0].next; lpj->n < mini; lpj = lpj->next)
        {
            da = da + join[lpj->n];
        }
        sumRows[mini] = da;

        /* skip if 'lpj->n' is equal to 'mini' */
        if ((lpj != NULL) && (lpj->n == mini))
        {
            lpj = lpj->next;
        }

        /* calculate sumCols[mini] */
        da = 0.0;
        for (; lpj != NULL; lpj = lpj->next)
        {
            da = da + join[lpj->n];
        }
        sumCols[mini] = da;

        /*
         * Clean up sumRows[minj], sumCols[minj] and score matrix
         * related with 'minj'.
         */
        sumCols[minj] = sumRows[minj] = 0.0;
        for (j = 1; j <= lastSeq - firstSeq + 1; ++j)
        {
            distMat->SetAt(minj, j, 0.0);
            distMat->SetAt(j, minj, 0.0);
            join[j] = 0.0;
        }

    } /** end main cycle **/

    /******************************Last Cycle (3 Seqs. left)********************/

    nude = 1;

    for (lpi = tvalid[0].next; lpi != NULL; lpi = lpi->next)
    {
        l[nude] = lpi->n;
        ++nude;
    }

    b1 = ((*distMat)(l[1], l[2]) + (*distMat)(l[1], l[3]) - (*distMat)(l[2], l[3])) *0.5;
    b2 = (*distMat)(l[1], l[2]) - b1;
    b3 = (*distMat)(l[1], l[3]) - b1;

    branch[1] = b1 - av[l[1]];
    branch[2] = b2 - av[l[2]];
    branch[3] = b3 - av[l[3]];

    /* Reset tiny negative and positive branch lengths to zero */
    if (fabs(branch[1]) < 0.0001)
    {
        branch[1] = 0.0;
    }
    if (fabs(branch[2]) < 0.0001)
    {
        branch[2] = 0.0;
    }
    if (fabs(branch[3]) < 0.0001)
    {
        branch[3] = 0.0;
    }

    phyTree->leftBranch[lastSeq - firstSeq + 1-2] = branch[1];
    phyTree->leftBranch[lastSeq - firstSeq + 1-1] = branch[2];
    phyTree->leftBranch[lastSeq - firstSeq + 1] = branch[3];

    for (i = 1; i <= lastSeq - firstSeq + 1; i++)
    {
        phyTree->treeDesc[lastSeq - firstSeq + 1-2][i] = 0;
    }

    if (verbose)
    {
        (*log) << "\n Cycle" << setw(4) << nc << " (Last cycle, trichotomy):\n";
    }

    for (i = 1; i <= 3; ++i)
    {
        if (av[l[i]] > 0.0)
        {
            if (verbose)
            {
                (*log) << "\n\t\t Node:" << setw(4) << l[i] <<" (" << setw(9) 
                        << setprecision(5) << branch[i] << ") ";
            }
            for (k = lastSeq - firstSeq + 1-3; k >= 1; k--)
                if (phyTree->treeDesc[k][l[i]] == 1)
                {
                    for (j = 1; j <= lastSeq - firstSeq + 1; j++)
                        if (phyTree->treeDesc[k][j] == 1)
                        {
                            phyTree->treeDesc[lastSeq - firstSeq + 1-2][j] = i;
                        }
                    break;
                }
        }
        else
        {
            if (verbose)
            {
                (*log) << "\n\t\t  SEQ:" << setw(4) << l[i] << " (" << setw(9) 
                        << setprecision(5) << branch[i] << ") ";
            }
            phyTree->treeDesc[lastSeq - firstSeq + 1-2][l[i]] = i;
        }
        if (i < 3)
        {
            if (verbose)
            {
                (*log) << "joins";
            }
        }
    }

    if (verbose)
    {
        (*log) << "\n";
    }

    /* IMPROVEMENT 2, STEP 4 : release memory area */

    delete [] tvalid;
    delete [] sumCols;
    delete [] sumRows;
    delete [] join;
    

}
}
