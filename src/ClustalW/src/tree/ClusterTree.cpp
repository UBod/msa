/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ClusterTree.h"
#include "../general/utils.h"
#include "dayhoff.h"
#include "RandomGenerator.h"
#include <math.h>
#include <sstream>
#include "../general/OutputFile.h"

namespace clustalw
{

ClusterTree::ClusterTree()
 : numSeqs(0),
   firstSeq(0),
   lastSeq(0)
{
    bootstrapPrompt = "\nEnter name for bootstrap output file  ";
    bootstrapFileTypeMsg = "Bootstrap output";
}



/** ****************************************************************************************
 *                          The rest are private functions.                                *
 *                                                                                         *
 *                                                                                         *
 *                                                                                         *
 *******************************************************************************************/
 
 
 
void ClusterTree::distanceMatrixOutput(ofstream* outFile, clustalw::DistMatrix* matToPrint,
                                        clustalw::Alignment *alignPtr)
{
    if(outFile == NULL || !outFile->is_open())
    {
        clustalw::utilityObject->error("Cannot output the distance matrix, file is not open\n");        
        return;
    }
     
    int i, j;
    int _maxNames = alignPtr->getMaxNames();
    (*outFile) << setw(6) << lastSeq - firstSeq + 1;
    
    for (i = 1; i <= lastSeq - firstSeq + 1; i++)
    {
        (*outFile) << "\n" << left << setw(_maxNames) << alignPtr->getName(i) << " ";
        for (j = 1; j <= lastSeq - firstSeq + 1; j++)
        {
            (*outFile) << " " << setw(6) << setprecision(3) << fixed << (*matToPrint)(i, j);
            if (j % 8 == 0)
            {
                if (j != lastSeq - firstSeq + 1)
                {
                    (*outFile) << "\n";
                }
                if (j != lastSeq - firstSeq + 1)
                {
                    (*outFile) << "          ";
                }
            }
        }
    }
}

void ClusterTree::overspillMessage(int overspill,int totalDists)
{
    std::ostringstream ssOverSpill;
    std::ostringstream ssTotalDists;
    string message;

    ssOverSpill << overspill;
    message += ssOverSpill.str();
    message += " of the distances out of a total of ";
    ssTotalDists << totalDists;
    message += ssTotalDists.str();
    message += "\n were out of range for the distance correction.\n"
      "\n SUGGESTIONS: 1) remove the most distant sequences"
      "\n           or 2) use the PHYLIP package"
      "\n           or 3) turn off the correction."
      "\n Note: Use option 3 with caution! With this degree"
      "\n of divergence you will have great difficulty"
      "\n getting robust and reliable trees.\n\n";

    clustalw::utilityObject->warning(message.c_str()); 
}

/*
 * The function treeGapDelete flags all the positions that have a gap in any sequence.
 *
 */
// NOTE there is something wrong with using _lenFirstSeq. But this is what old clustal does.
void ClusterTree::treeGapDelete(clustalw::Alignment *alignPtr)
{
    int seqn;
    int posn;
    int _maxAlnLength = alignPtr->getMaxAlnLength();
    int _lenFirstSeq = alignPtr->getSeqLength(firstSeq);
    int _gapPos1 = clustalw::userParameters->getGapPos1();
    int _gapPos2 = clustalw::userParameters->getGapPos2();
    
    treeGaps.resize(_maxAlnLength + 1);
    
    for (posn = 1; posn <= _lenFirstSeq; ++posn)
    {
        treeGaps[posn] = 0;
        for (seqn = 1; seqn <= lastSeq - firstSeq + 1; ++seqn)
        {
            const vector<int>* _seqM = alignPtr->getSequence(seqn + firstSeq - 1);
            
            if(posn > alignPtr->getSeqLength(seqn + firstSeq - 1))
            {
                break; // Dont read locations that cannot be read!
            }
            if (((*_seqM)[posn] == _gapPos1) || ((*_seqM)[posn] == _gapPos2))
            {
                treeGaps[posn] = 1;
                break;
            }
        }
    }
}

int ClusterTree::dnaDistanceMatrix(ofstream* treeFile, clustalw::Alignment *alignPtr)
{
    int m, n;
    int j, i;
    int res1, res2;
    int overspill = 0;
    double p, q, e, a, b, k;
    
    treeGapDelete(alignPtr); // flag positions with gaps (tree_gaps[i] = 1 ) 

    if (verbose)
    {
        (*treeFile) << "\n";
        (*treeFile) <<  "\n DIST   = percentage divergence (/100)";
        (*treeFile) << "\n p      = rate of transition (A <-> G; C <-> T)";
        (*treeFile) << "\n q      = rate of transversion";
        (*treeFile) << "\n Length = number of sites used in comparison";
        (*treeFile) << "\n";
        if (clustalw::userParameters->getTossGaps())
        {
            (*treeFile) << "\n All sites with gaps (in any sequence) deleted!";
            (*treeFile) << "\n";
        }
        if (clustalw::userParameters->getKimura())
        {
            (*treeFile) << "\n Distances corrected by Kimura's 2 parameter model:";
            (*treeFile) << "\n\n Kimura, M. (1980)";
            (*treeFile) << " A simple method for estimating evolutionary ";
            (*treeFile) << "rates of base";
            (*treeFile) << "\n substitutions through comparative studies of ";
            (*treeFile) << "nucleotide sequences.";
            (*treeFile) << "\n J. Mol. Evol., 16, 111-120.";
            (*treeFile) << "\n\n";
        }
    }
    
    int _numSeqs = alignPtr->getNumSeqs();
    quickDistMat.reset(new clustalw::DistMatrix(_numSeqs + 1));
    int _lenFirstSeq = alignPtr->getSeqLength(firstSeq);
    int _gapPos1 = clustalw::userParameters->getGapPos1();
    int _gapPos2 = clustalw::userParameters->getGapPos2();
    int lenSeqM, lenSeqN;
     // for every pair of sequence 
    for (m = 1; m < lastSeq - firstSeq + 1; ++m)
    {
        const vector<int>* _seqM = alignPtr->getSequence(m + firstSeq - 1);
        lenSeqM = alignPtr->getSeqLength(m + firstSeq - 1);
        
        for (n = m + 1; n <= lastSeq - firstSeq + 1; ++n)
        {
            const vector<int>* _seqN = alignPtr->getSequence(n + firstSeq - 1);
            lenSeqN = alignPtr->getSeqLength(n + firstSeq - 1);
            
            p = q = e = 0.0;
            quickDistMat->SetAt(m, n, 0.0);
            quickDistMat->SetAt(n ,m, 0.0);
            
            for (i = 1; i <= _lenFirstSeq; ++i)
            {
                j = bootPositions[i];
                if (clustalw::userParameters->getTossGaps() && (treeGaps[j] > 0))
                {
                    goto skip;
                }
                // gap position 
/** *******************************************************************************
 * BUG!!!!!!!      NOTE this was found for protein. Presuming the same here       *
 * NOTE: the following if statements were coded in so as to produce               *
 * the same distance results as the old clustal. Old clustal compares             *
 * up to the length of the first sequence. If this is longer than the             *
 * other sequences, then the -3 and 0's are compared at the end of the            *
 * array. These should not be compared, but I need to stick to this to            *
 * produce the same results as the old version for testing!                       *
 **********************************************************************************/
                if(j > lenSeqM)
                {
                    if(j == lenSeqM + 1)
                    {
                        res1 = -3;
                    }
                    else
                    {    
                        res1 = 0;
                    }
                }
                else
                { 
                    res1 = (*_seqM)[j];
                }
                if(j > lenSeqN)
                {
                    if(j == lenSeqN + 1)
                    {
                        res2 = -3;
                    }
                    else
                    {    
                        res2 = 0;
                    }
                }
                else
                {
                    res2 = (*_seqN)[j];
                }
                if ((res1 == _gapPos1) || (res1 == _gapPos2) || (res2 == _gapPos1)
                    || (res2 == _gapPos2))
                {
                    goto skip;
                }
                // gap in a seq
                if (!clustalw::userParameters->getUseAmbiguities())
                {
                    if (isAmbiguity(res1) || isAmbiguity(res2))
                    {
                        goto skip;
                    }
                }
                // ambiguity code in a seq
                e = e+1.0;
                if (res1 != res2)
                {
                    if (transition(res1, res2))
                    {
                        p = p + 1.0;
                    }
                    else
                    {
                        q = q + 1.0;
                    }
                }
                skip: ;
            }


            // Kimura's 2 parameter correction for multiple substitutions

            if (!clustalw::userParameters->getKimura())
            {
                if (e == 0)
                {
                    cerr << "\n WARNING: sequences " << m << " and " << n 
                         << " are non-overlapping\n";
                    k = 0.0;
                    p = 0.0;
                    q = 0.0;
                }
                else
                {
                    k = (p + q) / e;
                    if (p > 0.0)
                    {
                        p = p / e;
                    }
                    else
                    {
                        p = 0.0;
                    }
                    if (q > 0.0)
                    {
                        q = q / e;
                    }
                    else
                    {
                        q = 0.0;
                    }
                }
                quickDistMat->SetAt(m, n, k);
                quickDistMat->SetAt(n ,m, k);
                if (verbose) 
                {
                    (*treeFile) << setw(4) << m << " vs." << setw(4) << n << ":  DIST =  " 
                                << setw(4) << fixed << setprecision(4)
                                << k << "; p = " << fixed << setprecision(4) << p << "; q = "
                                << fixed << setprecision(4) << q << "; length = " << setw(6)
                                << fixed << setprecision(0) << e << "\n"; 
                }
            }
            else
            {
                if (e == 0)
                {
                    cerr << "\n WARNING: sequences " << m << " and " << n 
                         << " are non-overlapping\n";;
                    p = 0.0;
                    q = 0.0;
                }
                else
                {
                    if (p > 0.0)
                    {
                        p = p / e;
                    }
                    else
                    {
                        p = 0.0;
                    }
                    if (q > 0.0)
                    {
                        q = q / e;
                    }
                    else
                    {
                        q = 0.0;
                    }
                }

                if (((2.0 *p) + q) == 1.0)
                {
                    a = 0.0;
                }
                else
                {
                    a = 1.0 / (1.0 - (2.0 *p) - q);
                }

                if (q == 0.5)
                {
                    b = 0.0;
                }
                else
                {
                    b = 1.0 / (1.0 - (2.0 *q));
                }

                // watch for values going off the scale for the correction.
                if ((a <= 0.0) || (b <= 0.0))
                {
                    overspill++;
                    k = 3.5; // arbitrary high score
                }
                else
                {
                    k = 0.5 * log(a) + 0.25 * log(b);
                }
                quickDistMat->SetAt(m, n, k);
                quickDistMat->SetAt(n ,m, k);
                if (verbose)
                // if screen output 
                {
                    (*treeFile) << setw(4) << m << " vs." << setw(4) << n 
                                << ":  DIST =  " << fixed << setprecision(4)
                                << k << "; p = " << fixed <<  setprecision(4) << p << "; q = "
                                << fixed << setprecision(4) << q << "; length = " << setw(6)
                                << fixed << setprecision(0) << e << "\n"; 
                }

            }
        }
    }
    return overspill; // return the number of off-scale values
}

int ClusterTree::protDistanceMatrix(ofstream* treeFile, clustalw::Alignment *alignPtr)
{
    int m, n;
    int j, i;
    int res1, res2;
    int overspill = 0;
    double p, e, k, tableEntry;

    treeGapDelete(alignPtr); // flag positions with gaps (tree_gaps[i] = 1 ) 

    if (verbose)
    {
        (*treeFile) <<  "\n";
        (*treeFile) << "\n DIST   = percentage divergence (/100)";
        (*treeFile) << "\n Length = number of sites used in comparison";
        (*treeFile) << "\n\n";
        if (clustalw::userParameters->getTossGaps())
        {
            (*treeFile) << "\n All sites with gaps (in any sequence) deleted";
            (*treeFile) << "\n";
        }
        if (clustalw::userParameters->getKimura())
        {
            (*treeFile) <<
                "\n Distances up to 0.75 corrected by Kimura's empirical method:";
            (*treeFile) << "\n\n Kimura, M. (1983)";
            (*treeFile) << " The Neutral Theory of Molecular Evolution.";
            (*treeFile) <<
                "\n Page 75. Cambridge University Press, Cambridge, England.";
            (*treeFile) << "\n\n";
        }
    }
    int _numSeqs = alignPtr->getNumSeqs();
    int _lenSeq1 = alignPtr->getSeqLength(1);
    quickDistMat.reset(new clustalw::DistMatrix(_numSeqs + 1));
    int _gapPos1 = clustalw::userParameters->getGapPos1();
    int _gapPos2 = clustalw::userParameters->getGapPos2();
    int lenSeqM, lenSeqN;
    // for every pair of sequence 
    for (m = 1; m < _numSeqs; ++m)
    {
        const vector<int>* _seqM = alignPtr->getSequence(m);
        lenSeqM = alignPtr->getSeqLength(m);
        for (n = m + 1; n <= _numSeqs; ++n)
        {
            const vector<int>* _seqN = alignPtr->getSequence(n);
            lenSeqN = alignPtr->getSeqLength(n);
            p = e = 0.0;
            quickDistMat->SetAt(m, n, 0.0);
            quickDistMat->SetAt(n ,m, 0.0);
            for (i = 1; i <= _lenSeq1; ++i) // It may be this here!
            {
                j = bootPositions[i];
                if (clustalw::userParameters->getTossGaps() && (treeGaps[j] > 0))
                {
                    goto skip;
                }
                // gap position
/** *******************************************************************************
 * BUG!!!!!!!                                                                     *
 * NOTE: the following if statements were coded in so as to produce               *
 * the same distance results as the old clustal. Old clustal compares             *
 * up to the length of the first sequence. If this is longer than the             *
 * other sequences, then the -3 and 0's are compared at the end of the            *
 * array. These should not be compared, but I need to stick to this to            *
 * produce the same results as the old version for testing!                       *
 **********************************************************************************/
                if(j > lenSeqM)
                {
                    if(j == lenSeqM + 1)
                    {
                        res1 = -3;
                    }
                    else
                    {    
                        res1 = 0;
                    }
                }
                else
                { 
                    res1 = (*_seqM)[j];
                }
                if(j > lenSeqN)
                {
                    if(j == lenSeqN + 1)
                    {
                        res2 = -3;
                    }
                    else
                    {    
                        res2 = 0;
                    }
                }
                else
                {
                    res2 = (*_seqN)[j];
                }
                if ((res1 == _gapPos1) || (res1 == _gapPos2) || (res2 == _gapPos1)
                    || (res2 == _gapPos2))
                {
                    goto skip;
                }
                // gap in a seq
                e = e + 1.0;
                if (res1 != res2)
                {
                    p = p + 1.0;
                }
                skip: ;
            }

            if (p <= 0.0)
            {
                k = 0.0;
            }
            else
            {
                k = p / e;
            }

            if (clustalw::userParameters->getKimura())
            {
                if (k < 0.75)
                {
                    // use Kimura's formula
                    if (k > 0.0)
                    {
                        k =  - log(1.0 - k - (k * k / 5.0));
                    }
                }
                else
                {
                    if (k > 0.930)
                    {
                        overspill++;
                        k = 10.0; // arbitrarily set to 1000%
                    }
                    else // dayhoff_pams is from dayhoff.h file
                    {
                        tableEntry = (k * 1000.0) - 750.0;
                        k = (double)dayhoff_pams[(int)tableEntry];
                        k = k / 100.0;
                    }
                }
            }

            quickDistMat->SetAt(m, n, k);
            quickDistMat->SetAt(n ,m, k);
            if (verbose)
            {
                (*treeFile) << setw(4) << m << " vs." << setw(4) << n 
                            << "  DIST = " << fixed << setprecision(4)
                            << k << ";  length = " << setw(6)
                            << setprecision(0) << e << "\n";
            }
        }
    }
    return overspill;
}

bool ClusterTree::isAmbiguity(int c)
{
    int i;
    char codes[] = "ACGTU";

    if (clustalw::userParameters->getUseAmbiguities() == true)
    {
        return false;
    }

    for (i = 0; i < 5; i++)
        if (clustalw::userParameters->getAminoAcidCode(c) == codes[i])
        {
            return false;
        }

    return true;
}

/*
 * The function calcPercIdentity calculates the percent identity of the sequences
 * and outputs it to a the file pfile. NOTE this is not used at the moment. It was in
 * the old code, but there was no way to access it from the menu. This may change. 
 */
void ClusterTree::calcPercIdentity(ofstream* pfile, clustalw::Alignment *alignPtr)
{
    clustalw::DistMatrix percentMat;
  
    float ident;
    int nmatch;
  
    int val1, val2;
  
    int i,j,k, length_longest;
    int length_shortest;
    
    int rs = 0, rl = 0;
    // findout sequence length, longest and shortest;
    length_longest = 0;
    length_shortest = 0;

    int _numSeqs = alignPtr->getNumSeqs();
    int _seqLength;
    for (i = 1; i <= _numSeqs; i++) 
    {
        _seqLength = alignPtr->getSeqLength(i);
        if (length_longest < _seqLength)
        {
            length_longest = _seqLength;
            rs = i;
        }
        if (length_shortest > _seqLength) 
        {
            length_shortest = _seqLength;
            rl = i;
        }
    } 

    percentMat.ResizeRect(_numSeqs + 1);
    nmatch = 0;
    int _lenSeqI, _lenSeqJ;
    int _maxAA = clustalw::userParameters->getMaxAA();
    
    for (i = 1; i <= _numSeqs; i++) 
    {
        const vector<int>* _seqI = alignPtr->getSequence(i);
        _lenSeqI = alignPtr->getSeqLength(i);
        
        for (j = i; j <= _numSeqs; j++) 
        {
            const vector<int>* _seqJ = alignPtr->getSequence(j);
            _lenSeqJ = alignPtr->getSeqLength(j);
            
            cout << "\n           " << alignPtr->getName(j) << " ";
            ident = 0;
            nmatch = 0;
            for(k = 1; k <= length_longest; k++) 
            {
                if((k > _lenSeqI) || (k > _lenSeqJ))
                {
                    break;
                }
                val1 = (*_seqI)[k];
                val2 = (*_seqJ)[k];

                if (((val1 < 0) || (val1 > _maxAA)) || ((val2 < 0) || (val2 > _maxAA)))
                { 
                    continue; // residue = '-';
                }
                if (val1 == val2) 
                {
                    ident++ ;
                    nmatch++;
                }
                else 
                {
                    nmatch++ ;
                }
            }
            ident = ident/nmatch * 100.0 ;
            percentMat.SetAt(i, j, ident);
            percentMat.SetAt(j, i, ident);
        }
    }

    int _maxNameSize = alignPtr->getMaxNames();
    
    (*pfile) << "#\n#\n#  Percent Identity  Matrix - created by Clustal"
             << clustalw::userParameters->getRevisionLevel() << " \n#\n#\n";
    for(i = 1; i <= _numSeqs; i++) 
    {
        (*pfile) << "\n " << right << setw(5) << i << ": "; 
        (*pfile) << left << setw(_maxNameSize) << alignPtr->getName(i);
        
        for(j = 1; j <= _numSeqs; j++) 
        {
            (*pfile) << setw(8) << right << fixed << setprecision(0) << percentMat(i, j);
        }
    }
    (*pfile) << "\n";

}

void ClusterTree::compareTree(PhyloTree* tree1, PhyloTree* tree2, vector<int>* hits, int n)
{
    int i, j, k;
    int nhits1, nhits2;

    for (i = 1; i <= n - 3; i++)
    {
        for (j = 1; j <= n - 3; j++)
        {
            nhits1 = 0;
            nhits2 = 0;
            for (k = 1; k <= n; k++)
            {
                if (tree1->treeDesc[i][k] == tree2->treeDesc[j][k])
                {
                    nhits1++;
                }
                if (tree1->treeDesc[i][k] != tree2->treeDesc[j][k])
                {
                    nhits2++;
                }
            }
            if ((nhits1 == lastSeq - firstSeq + 1) || (nhits2 == lastSeq -
                firstSeq + 1))
            {
                (*hits)[i]++;
            }
        }
    }
}

/**
 * NOTE this will go into the OutputFile class and will not be needed here anymore.
 *
 */
/*string ClusterTree::getOutputFileName(const string prompt, string path, 
                                      const string fileExtension)
{
    string temp;
    string _fileName; // Will return this name.
    string message;
    _fileName = path + fileExtension;

    if(_fileName.compare(clustalw::userParameters->getSeqName()) == 0) 
    {
        cout << "Output file name is the same as input file.\n";
        if (clustalw::userParameters->getMenuFlag()) 
        {
            message = "\n\nEnter new name to avoid overwriting  [" + _fileName + "]: ";
            clustalw::utilityObject->getStr(message, temp);
            if(temp != "")
            {
                _fileName = temp;
            }
        }
    }
    else if (clustalw::userParameters->getMenuFlag()) 
    {

        message = prompt + " [" + _fileName + "]";
        clustalw::utilityObject->getStr(message, temp);
        if(temp != "")
        {
            _fileName = temp;
        }
    }   
    return _fileName;

}*/

bool ClusterTree::transition(int base1, int base2)
{
// assumes that the bases of DNA sequences have been translated as
// a,A = 0;   c,C = 1;   g,G = 2;   t,T,u,U = 3;  N = 4;
// a,A = 0;   c,C = 2;   g,G = 6;   t,T,u,U =17;
// A <--> G  and  T <--> C  are transitions;  all others are transversions.
    if (((base1 == 0) && (base2 == 6)) || ((base1 == 6) && (base2 == 0)))
    {
        return true;
    }
    // A <--> G 
    if (((base1 == 17) && (base2 == 2)) || ((base1 == 2) && (base2 == 17)))
    {
        return true;
    }
    // T <--> C 
    return false;
}

/**
 * This function is used to open all the bootstrap tree files. It opens them with the 
 * correct message prompt.
 */
bool ClusterTree::openFilesForBootstrap(clustalw::OutputFile* clustalFile, clustalw::OutputFile* phylipFile,
                         clustalw::OutputFile* nexusFile, clustalw::TreeNames* treeNames, string* path)
{
    if (clustalw::userParameters->getOutputTreeClustal())
    {
        if(!clustalFile || !clustalFile->openFile(&(treeNames->clustalName), 
                                  bootstrapPrompt, path, "njb", bootstrapFileTypeMsg))
        {
            return false;
        } 
    }   
        
    if (clustalw::userParameters->getOutputTreePhylip())
    {     
        if(!phylipFile || !phylipFile->openFile(&(treeNames->phylipName), 
                                bootstrapPrompt, path, "phb", bootstrapFileTypeMsg))
        {
            return false;
        }                     
    }    

    if (clustalw::userParameters->getOutputTreeNexus())
    {
        if(!nexusFile || !nexusFile->openFile(&(treeNames->nexusName), 
                                bootstrapPrompt, path, "treb", bootstrapFileTypeMsg))
        {
            return false;
        }                   
    }        
    return true;
}                         

bool ClusterTree::openFilesForTreeFromAlignment(clustalw::OutputFile* clustalFile, 
    clustalw::OutputFile* phylipFile, clustalw::OutputFile* distFile, clustalw::OutputFile* nexusFile, clustalw::OutputFile* pimFile, 
                            clustalw::TreeNames* treeNames, string* path)
{
    if (clustalw::userParameters->getOutputTreeClustal())
    {
        if(!clustalFile || !clustalFile->openFile(&(treeNames->clustalName), 
                                  "\nEnter name for CLUSTAL    tree output file  ",
                                  path, "nj", "Phylogenetic tree"))
        {
            return false;
        } 
    }
        
    if (clustalw::userParameters->getOutputTreePhylip())
    {     
        if(!phylipFile || !phylipFile->openFile(&(treeNames->phylipName), 
                             "\nEnter name for PHYLIP     tree output file  ", path, "ph",
                             "Phylogenetic tree"))
        {
            return false;
        }
    }

    if (clustalw::userParameters->getOutputTreeDistances())
    {   
        if(!distFile || !distFile->openFile(&(treeNames->distName), 
                               "\nEnter name for distance matrix output file  ",
                               path, "dst", "Distance matrix"))
        {
            return false;
        }         
    }

    if (clustalw::userParameters->getOutputTreeNexus())
    {
        if(!nexusFile || !nexusFile->openFile(&(treeNames->nexusName), 
                                "\nEnter name for NEXUS tree output file  ", path,
                                     "tre", "Nexus tree"))
        {
            return false;
        }
    }

    if (clustalw::userParameters->getOutputPim())
    {       
        if(!pimFile || !pimFile->openFile(&(treeNames->pimName), 
                            "\nEnter name for % Identity matrix output file  ", path, "pim", 
                            "perc identity"))
        {
            return false;
        }                       
    }
    return true;
}

int ClusterTree::calcQuickDistMatForAll(ofstream* clustalFile, ofstream* phylipFile, 
                                  ofstream* nexusFile, ofstream* pimFile, ofstream* distFile,
                                  clustalw::Alignment* alignPtr)
{
    int overspill = 0;
    bool _DNAFlag = clustalw::userParameters->getDNAFlag();
    
    overspill = calcQuickDistMatForSubSet(clustalFile, phylipFile, nexusFile, alignPtr);    

    if (pimFile && clustalw::userParameters->getOutputPim())
    {
        verbose = false; // Turn off file output
        if (_DNAFlag)
        {
            calcPercIdentity(pimFile, alignPtr);
        }
        else
        {
            calcPercIdentity(pimFile, alignPtr);
        }
    }


    if (distFile && clustalw::userParameters->getOutputTreeDistances())
    {
        verbose = false; // Turn off file output
        if (_DNAFlag)
        {
            overspill = dnaDistanceMatrix(distFile, alignPtr);
        }
        else
        {
            overspill = protDistanceMatrix(distFile, alignPtr);
        }
        distanceMatrixOutput(distFile, quickDistMat.get(),
                             alignPtr);
    }
    return overspill;
}                    

int ClusterTree::calcQuickDistMatForSubSet(ofstream* clustalFile, ofstream* phylipFile, 
                                          ofstream* nexusFile, clustalw::Alignment* alignPtr, 
                                          bool inBootLoop)
{
    int overspill = 0;
    bool _DNAFlag = clustalw::userParameters->getDNAFlag();
    
    if (clustalFile && clustalw::userParameters->getOutputTreeClustal())
    {
        if(!inBootLoop)
        {
            verbose = true; // Turn on file output
        }
        else
        {
            verbose = false; // Turn off when we are in the loop in bootstrap!
        }   
        if (_DNAFlag)
        {
            overspill = dnaDistanceMatrix(clustalFile, alignPtr);
        }
        else
        {
            overspill = protDistanceMatrix(clustalFile, alignPtr);
        }
    }

    if (phylipFile && clustalw::userParameters->getOutputTreePhylip())
    {
        verbose = false; // Turn off file output 
        if (_DNAFlag)
        {
            overspill = dnaDistanceMatrix(phylipFile, alignPtr);
        }
        else
        {
            overspill = protDistanceMatrix(phylipFile, alignPtr);
        }
    }

    if (nexusFile && clustalw::userParameters->getOutputTreeNexus())
    {
        verbose = false; // Turn off file output 
        if (_DNAFlag)
        {
            overspill = dnaDistanceMatrix(nexusFile, alignPtr);
        }
        else
        {
            overspill = protDistanceMatrix(nexusFile, alignPtr);
        }
    }
    return overspill;    
}

void ClusterTree::printBootstrapHeaderToClustalFile(ofstream* clustalFile)
{
    if(clustalFile)
    {
        (*clustalFile) << "\n\n\t\t\tBootstrap Confidence Limits\n\n";
        (*clustalFile) << "\n Random number generator seed = " 
                       << setw(7)
                       << clustalw::userParameters->getBootRanSeed() << "\n";
        (*clustalFile) << "\n Number of bootstrap trials   = " << setw(7)
                       << clustalw::userParameters->getBootNumTrials() << "\n";
        (*clustalFile) << "\n\n Diagrammatic representation of the above tree: \n";
        (*clustalFile) << "\n Each row represents 1 tree cycle;";
        (*clustalFile) << " defining 2 groups.\n";
        (*clustalFile) << "\n Each column is 1 sequence; ";
        (*clustalFile) << "the stars in each line show 1 group; ";
        (*clustalFile) << "\n the dots show the other\n";
        (*clustalFile) << "\n Numbers show occurences in bootstrap samples.";
    }
}

void ClusterTree::promptForBoolSeedAndNumTrials()
{
    if (clustalw::userParameters->getMenuFlag())
    {
        unsigned int tempSeed;
        tempSeed = clustalw::utilityObject->getInt(
                "\n\nEnter seed no. for random number generator ", 1, 1000,
        clustalw::userParameters->getBootRanSeed());
        clustalw::userParameters->setBootRanSeed(tempSeed);

        clustalw::userParameters->setBootNumTrials(
                     clustalw::utilityObject->getInt("\n\nEnter number of bootstrap trials ", 
                        1, 10000, clustalw::userParameters->getBootNumTrials()));    
    }
}

void ClusterTree::printErrorMessageForBootstrap(int totalOverspill, int totalDists, 
                                                int nfails)
{
    cerr << "\n";
    cerr << "\n WARNING: " << totalOverspill 
         << " of the distances out of a total of " 
         << totalDists << " times" << clustalw::userParameters->getBootNumTrials();
    cerr << "\n were out of range for the distance correction.";
    cerr << "\n This affected " << nfails << " out of " 
         << clustalw::userParameters->getBootNumTrials() << " bootstrap trials.";
    cerr << "\n This may not be fatal but you have been warned!" << "\n";
    cerr << "\n SUGGESTIONS: 1) turn off the correction";
    cerr << "\n           or 2) remove the most distant sequences";
    cerr << "\n           or 3) use the PHYLIP package.\n\n";
            
    if (clustalw::userParameters->getMenuFlag())
    {
        string dummy;
        clustalw::utilityObject->getStr(string("Press [RETURN] to continue"), dummy);
    }
}

bool ClusterTree::checkIfConditionsMet(int numSeqs, int min)
{
    if (clustalw::userParameters->getEmpty())
    {
        clustalw::utilityObject->error("You must load an alignment first");
        return false;
    }

    if (numSeqs < min)
    {
        clustalw::utilityObject->error("Alignment has only %d sequences", numSeqs);
        return false;
    }
    
    return true;
}
                                 
}
