/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This file contains the implementation of the UserParameter functions
 * Mark Larkin Dec 8 2005
 *
 * Modified: 17 January 2008 Paul McGettigan added in the O aminoacid residue to handle Pyrrolysine
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <stdio.h>
#include <string>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <climits>
#include <iomanip>
#include <fstream>
#include "UserParameters.h"
#include "clustalw_version.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
namespace clustalw
{
using namespace std;
 
UserParameters::UserParameters(bool log)
{
   // FIXME: huge parts should be merged/replaced with
   // setParamsToDefault (which is not used at all)
   
    gapPos1 = NUMRES - 2; /* code for gaps inserted by clustalw */
    gapPos2 = NUMRES - 1; /* code for gaps already in alignment */
    
    //revisionLevel = CLU_SHORT_VERSION_STRING;
    revisionLevel = CLUSTALW_VERSION;
    interactive = false;
    gui = false;

    seqName = "";
    DNAGapOpen = 15.0;
    DNAGapExtend = 6.66;
    AAGapOpen = 10.0;
    AAGapExtend = 0.2;
    gapDist = 4;
    outputOrder = ALIGNED; // Note: All macros should be replaced by const
    divergenceCutoff = 30;

    hydResidues = "GPSNDQEKR";
    noWeights = false;
    negMatrix = false;
    noHydPenalties = false;
    noVarPenalties = true;
    noPrefPenalties = false;
    useEndGaps = false;
    endGapPenalties = false;
    resetAlignmentsNew = false;
    resetAlignmentsAll = false;
    outputStructPenalties = OUTSECST;
    structPenalties1 = NONE;
    structPenalties2 = NONE;
    useSS1 = true;
    useSS2 = true;
    helixPenalty = 4;
    strandPenalty = 4;
    loopPenalty = 1;
    helixEndMinus = 3;
    helixEndPlus = 0;
    strandEndMinus = 1;
    strandEndPlus = 1;
    helixEndPenalty = 2;
    strandEndPenalty = 2;
    useAmbiguities = false;
    DNAPWGapOpen = 15.0;
    DNAPWGapExtend = 6.66;
    AAPWGapOpen = 10.0;
    AAPWGapExtend = 0.1;

    quickPairAlign = false;
    transitionWeight = 0.5;
    DNAKtup = 2;
    DNAWindowGap = 5;
    DNASignif = 4;
    DNAWindow = 4;
    AAKtup = 1;
    AAWindowGap = 3;
    AASignif = 5;
    AAWindow = 5;
    percent = true;
    tossgaps = false;
    kimura = false;
    bootNumTrials = 1000;
    bootRanSeed = 111;
    debug = 0;
    explicitDNAFlag = false;
    lowercase = true;
    clSeqNumbers = false;
    seqRange = false;
    outputClustal = true;
    outputGcg = false;
    outputPhylip = false;
    outputNbrf = false;
    outputGde = false;
    outputNexus = false;
    outputFasta = false;
    showAlign = true;
    saveParameters = false;
    outputTreeClustal = false;
    outputTreePhylip = true;
    outputTreeDistances = false;
    outputTreeNexus = false;
    outputPim = false;
    bootstrapFormat = BS_BRANCH_LABELS;
    profile1Name = ""; // Initialise to blank strings
    profile2Name = "";
    empty = true;
    profile1Empty = true;
    profile2Empty = true;
    outfileName = "";
    //profile1NumSeqs = 0; // MARK: Set to default before used.
    useTreeFile = false;
    newTreeFile = false;
    useTree1File = false;
    useTree2File = false;
    newTree1File = false;
    newTree2File = false;
    aminoAcidCodes = "ABCDEFGHIKLMNOPQRSTUVWXYZ-";
    maxAA = aminoAcidCodes.length() - 2;
    
    // Some variables need the alignment to be read in before they can be set.
    // I am putting the default as protein. Note: this should not make a difference
    // as they are not used before they have been set a value again!!
    
    gapOpen = AAGapOpen;
    gapExtend = AAGapExtend;
    PWGapOpen = AAPWGapOpen;
    PWGapExtend = AAPWGapExtend;
    
    gapPos1 = NUMRES - 2;
    gapPos2 = NUMRES - 1;
    profileNum = 0;
    menuFlag = false; // MARK: I set to default value.
    DNAFlag = false; // MARK: I set to default value.
    distanceTree = true; // MARK: I set to default value.
    ktup = AAKtup;
    window = AAWindow;
    windowGap = AAWindowGap;
    signif = AASignif;
    rangeFrom = -1;
    rangeTo = -1;
    rangeFromToSet = false;
    QTscorePlotScale = 5;
    QTresExceptionCutOff = 5;
    
    QTseqWeightCalculated = false;
    QTminLenLowScoreSegment = 1;
    QTlowScoreDNAMarkingScale = 5;
    
    // Set defaults for iteration variables.
    numIterations = 3;
    doRemoveFirstIteration = NONE; 
    maxAllowedSeqLength = INT_MAX;
    
    clusterAlgorithm = NJ;
    displayInfo = true;
    helpFlag = false;
    fullHelpFlag = false;
    quiet = false;
    
}

// FIXME:never used
void UserParameters::setParamsToDefault()
{
    DNAGapOpen = 15.0;
    DNAGapExtend = 6.66;
    AAGapOpen = 10.0;
    AAGapExtend = 0.2;
    gapDist = 4;
    outputOrder = ALIGNED; 
    divergenceCutoff = 30;
    
    hydResidues = "GPSNDQEKR";
    noWeights = false;
    negMatrix = false;
    noHydPenalties = false;
    noVarPenalties = true;
    noPrefPenalties = false;
    useEndGaps = false;
    endGapPenalties = false;
    resetAlignmentsNew = false;
    resetAlignmentsAll = false;
    outputStructPenalties = OUTSECST;
    structPenalties1 = NONE;
    structPenalties2 = NONE;
    useSS1 = true;
    useSS2 = true;
    helixPenalty = 4;
    strandPenalty = 4;
    loopPenalty = 1;
    helixEndMinus = 3;
    helixEndPlus = 0;
    strandEndMinus = 1;
    strandEndPlus = 1;
    helixEndPenalty = 2;
    strandEndPenalty = 2;
    useAmbiguities = false;
    DNAPWGapOpen = 15.0;
    DNAPWGapExtend = 6.66;
    AAPWGapOpen = 10.0;
    AAPWGapExtend = 0.1;
    quickPairAlign = false;
    transitionWeight = 0.5;
    DNAKtup = 2;
    DNAWindowGap = 5;
    DNASignif = 4;
    DNAWindow = 4;
    AAKtup = 1;
    AAWindowGap = 3;
    AASignif = 5;
    AAWindow = 5;
    percent = true;
    tossgaps = false;
    kimura = false;
    bootNumTrials = 1000;
    bootRanSeed = 111;
    debug = 0;
    lowercase = true;
    clSeqNumbers = false;
    seqRange = false;
    outputClustal = true;
    outputGcg = false;
    outputPhylip = false;
    outputNbrf = false;
    outputGde = false;
    outputNexus = false;
    outputFasta = false;
    outputTreeClustal = false;
    outputTreePhylip = true;
    outputTreeDistances = false;
    outputTreeNexus = false;
    outputPim = false;
    bootstrapFormat = BS_BRANCH_LABELS;
    useTreeFile = false;
    newTreeFile = false;
    useTree1File = false;
    useTree2File = false;
    newTree1File = false;
    newTree2File = false;
    rangeFrom = -1;
    rangeTo = -1;
    rangeFromToSet = false;
    QTscorePlotScale = 5;
    QTresExceptionCutOff = 5;
    
    QTminLenLowScoreSegment = 1;
    QTlowScoreDNAMarkingScale = 5;
    distanceTree = true; // MARK: I set to default value.
    
    numIterations = 3;
                
    if(getDNAFlag())
    {
        setDNAParams();
    }
    else
    {
        setProtParams();
    }  
    
    clusterAlgorithm = NJ;   
    displayInfo = true; 
    helpFlag = false;
    fullHelpFlag = false;
    quiet = false;
    doRemoveFirstIteration = NONE;
    maxAllowedSeqLength = INT_MAX;
}

/*
 * The function createParameterOutput is used to put all the user parameters in
 * a file. It is used for testing and for saving parameters.
 *
 *
 * FIXME: AW: Some parameters are missing here (e.g. the new ones like
 * clustering, etc)
 *
 */
void UserParameters::createParameterOutput(void)
{
    string parname, temp;
    string path;
    string message;
  
    utilityObject->getPath(seqName, &path);
    parname = path + "par";
    if(menuFlag) 
    {
        message = "\nEnter a name for the parameter output file [" + parname + "]";
        utilityObject->getStr(message, temp);
        if(temp != "")
        {
            parname = temp;
        }
    }
      
    ofstream outfile;
    outfile.open(parname.c_str(), ofstream::out);
    
    if(!outfile)
    {
        return; // Failed to open
    }
    
    outfile << "clustalw \\\n";
    if (!empty && profile1Empty) 
    {
        outfile << "-infile=" << seqName << " \\\n";
    }
    if (!profile1Empty)
    {
        outfile << "-profile1=" << profile1Name << "\\\n";
    }
    if (!profile2Empty)
    {
        outfile << "-profile2=" << profile2Name << " \\\n";
    }
    if (DNAFlag == true)
    { 
        outfile << "-type=dna \\\n";
    }
    else
    {
        outfile << "-type=protein \\\n";
    }
    if (quickPairAlign) 
    {
        outfile << "-quicktree \\\n";
        outfile << "-ktuple=" << ktup << " \\\n";
        outfile << "-window=" << window << " \\\n";
        outfile << "-pairgap=" << windowGap << " \\\n";
        outfile << "-topdiags=" << signif << " \\\n";    
        if (percent)
        {
            outfile << "-score=percent \\\n";
        }      
        else
        {
            outfile << "-score=absolute \\\n";
        }      
    }
    else 
    {
        if (!DNAFlag) 
        {
            //outfile << "-pwmatrix=" << pwMatrixName << " \\\n";
            outfile << "-pwgapopen=" << fixed << setprecision(2) << AAPWGapOpen 
                    << " \\\n";
            outfile << "-pwgapext=" << AAPWGapExtend << " \\\n";
        }
        else 
        {
            outfile << "-pwgapopen=" << fixed << setprecision(2) << PWGapOpen << " \\\n";
            outfile << "-pwgapext=" << PWGapExtend << " \\\n";
        }
    }
  
    if (!DNAFlag) 
    {
        //outfile << "-matrix=" << matrixName << " \\\n";
        outfile << "-gapopen=" << fixed << setprecision(2) << AAGapOpen << " \\\n";
        outfile << "-gapext=" << AAGapExtend << " \\\n";
    }
    else 
    {
        outfile << "-gapopen=" << fixed << setprecision(2) << DNAGapOpen << " \\\n";
        outfile << "-gapext=" << DNAGapExtend << " \\\n";
    }
  
    outfile << "-maxdiv=" << divergenceCutoff << " \\\n";
    if (!useEndGaps) 
    {
        outfile << "-endgaps \\\n";
    }    
  
    if (!DNAFlag) 
    {
        if (negMatrix) 
        {
            outfile << "-negative \\\n";
        }   
        if (noPrefPenalties)
        { 
            outfile << "-nopgap \\\n";
        }     
        if (noHydPenalties) 
        { 
            outfile << "-nohgap \\\n";
        }     
        if (noVarPenalties) 
        {
            outfile << "-novgap \\\n";
        }     
        outfile << "-hgapresidues=" << hydResidues << " \\\n";
        outfile << "-gapdist=" << gapDist << " \\\n";     
    }
    else 
    {
        outfile << "-transweight=" << transitionWeight << " \\\n";
    }
  
    if (outputGcg) 
    {
        outfile << "-output=gcg \\\n";
    }
    else if (outputGde) 
    {
        outfile << "-output=gde \\\n";
    }
    else if (outputNbrf)
    {
        outfile << "-output=pir \\\n";
    }
    else if (outputPhylip) 
    {
        outfile << "-output=phylip \\\n";
    }
    else if (outputNexus) 
    {
        outfile << "-output=nexus \\\n";
    }
    if (outfileName[0]!=EOS)
    {    
        outfile << "-outfile=" << outfileName << " \\\n";
    }
    if (outputOrder==ALIGNED)
    {
        outfile << "-outorder=aligned \\\n";
    }  
    else
    {
        outfile << "-outorder=input \\\n";
    }  
    if (outputGde)
    {
        if (lowercase)
        {
            outfile << "-case=lower \\\n";
        }
        else
        {
            outfile << "-case=upper \\\n";
        }
    }
  
  
    outfile << "-interactive\n";

    outfile.close();

}

/*
 * The function resIndex returns the index of the character c in the string t.
 *
 */
int UserParameters::resIndex(string t, char c)
{
    register int i;

    for (i = 0; t[i] && t[i] != c; i++)
        ;
    if (t[i])
    {
        return (i);
    }
    else
    {
        return  -1;
    }
}

void UserParameters::setDNAMultiGap()
{
    gapOpen = DNAGapOpen;
    gapExtend = DNAGapExtend;
}

void UserParameters::setProtMultiGap()
{
    gapOpen = AAGapOpen;
    gapExtend = AAGapExtend;
}

void UserParameters::setDNAParams()
{
    gapOpen       = DNAGapOpen;
    gapExtend     = DNAGapExtend;
    PWGapOpen  = DNAPWGapOpen;
    PWGapExtend  = DNAPWGapExtend;
    ktup           = DNAKtup;
    window         = DNAWindow;
    signif         = DNASignif;
    windowGap       = DNAWindowGap;
}

void UserParameters::setProtParams()
{
    gapOpen       = AAGapOpen;
    gapExtend     = AAGapExtend;
    PWGapOpen  = AAPWGapOpen;
    PWGapExtend  = AAPWGapExtend;
    ktup           = AAKtup;
    window         = AAWindow;
    signif         = AASignif;
    windowGap       = AAWindowGap;
}

void UserParameters::setPWParamToProtein()
{
    PWGapOpen = AAPWGapOpen;
    PWGapExtend = AAPWGapExtend;
    ktup = AAKtup;
    window = AAWindow;
    signif = AASignif;
    windowGap = AAWindowGap;
}

void UserParameters::setPWParamToDNA()
{
    PWGapOpen = DNAPWGapOpen;
    PWGapExtend = DNAPWGapExtend;
    ktup = DNAKtup;
    window = DNAWindow;
    signif = DNASignif;
    windowGap = DNAWindowGap;
}

void UserParameters::setPWProteinParam()
{
    AAPWGapOpen = PWGapOpen;
    AAPWGapExtend = PWGapExtend;
    AAKtup = ktup;
    AAWindow = window;
    AASignif = signif;
    AAWindowGap = windowGap;
}

void UserParameters::setPWDNAParam()
{
    DNAPWGapOpen = PWGapOpen;
    DNAPWGapExtend = PWGapExtend;
    DNAKtup = ktup;
    DNAWindow = window;
    DNASignif = signif;
    DNAWindowGap = windowGap;
}
/*
 * The rest of the functions are get, set and toggle functions for the variables.
 */
 
string UserParameters::getRevisionLevel()
{
    return revisionLevel;
}

void UserParameters::setRevisionLevel(string value)
{
    revisionLevel = value;
}

void UserParameters::setInteractive(bool value)
{
    interactive = value;
}

void UserParameters::setGui(bool value)
{
    gui = value;
}

void UserParameters::setGapOpen(float value)
{
    gapOpen = value;
}

void UserParameters::setGapExtend(float value)
{
    gapExtend = value;
}

void UserParameters::setPWGapOpen(float value)
{
    PWGapOpen = value;
}

void UserParameters::setPWGapExtend(float value)
{
    PWGapExtend = value;
}

void UserParameters::setMaxAA(int value)
{
    maxAA = value;
}

void UserParameters::setGapPos1(int value)
{
    gapPos1 = value;
}

void UserParameters::setGapPos2(int value)
{
    gapPos2 = value;
}

void UserParameters::setProfileNum(int value)
{
    profileNum = value;
}

void UserParameters::setMenuFlag(bool value)
{
    menuFlag = value;
}

void UserParameters::setDNAFlag(bool value)
{
    if(value == true)
    {
        setDNAParams();
    }
    else
    {
        setProtParams();
    }    
    DNAFlag = value;
}

void UserParameters::setDistanceTree(bool value)
{
    distanceTree = value;
}

void UserParameters::setSeqName(string value)
{
    seqName = value;
}

void UserParameters::setDNAGapOpen(float value)
{
    DNAGapOpen = value;
}

void UserParameters::setDNAGapExtend(float value)
{
    DNAGapExtend = value;
}

void UserParameters::setProteinGapOpen(float value)
{
    AAGapOpen = value;
}

void UserParameters::setProteinGapExtend(float value)
{
    AAGapExtend = value;
}

void UserParameters::setGapDist(int value)
{
    gapDist = value;
}

void UserParameters::setOutputOrder(int value)
{
    outputOrder = value;
}

void UserParameters::toggleOutputOrder()
{
    if (outputOrder == INPUT)
    {
        outputOrder = ALIGNED;
    }
    else
    {
        outputOrder = INPUT;
    }
}

void UserParameters::setDivergenceCutoff(int value)
{
    divergenceCutoff = value;
}

void UserParameters::setHydResidues(string value)
{
    //hydResidues = value;
    char hydResidue;
    string tempHydRes = "";
    int inputStringLength = value.length();
    if(inputStringLength > 0) 
    {
    // NOTE this was causing an error, but I fixed it. Was giving an 
    // out of range error.
        for (int i = 0; i < MAXHYDRESIDUES && i < inputStringLength; i++) 
        {
            hydResidue = toupper(value.at(i));

            if (isalpha(hydResidue))
            {
                tempHydRes += hydResidue;
            }
            else // Not Alphabetic character!
            {
                break;
            }
        }
        if(tempHydRes.size() > 0)
        {
            hydResidues = tempHydRes;
        }
    }
}

void UserParameters::setNoWeights(bool value)
{
    noWeights = value;
}

void UserParameters::setUseNegMatrix(bool value)
{
    negMatrix = value;
}

void UserParameters::toggleUseNegMatrix()
{
    negMatrix ^= true;
}

void UserParameters::setNoHydPenalties(bool value)
{
    noHydPenalties = value;
}

void UserParameters::toggleNoHydPenalties()
{
    noHydPenalties ^= true;
}

void UserParameters::setNoVarPenalties(bool value)
{
    noVarPenalties = value;
}

void UserParameters::setNoPrefPenalties(bool value)
{
    noPrefPenalties = value;
}

void UserParameters::toggleNoPrefPenalties()
{
    noPrefPenalties ^= true;
}

void UserParameters::setUseEndGaps(bool value)
{
    useEndGaps = value;
}

void UserParameters::toggleUseEndGaps()
{
    useEndGaps ^= true;
}

void UserParameters::setEndGapPenalties(bool value)
{
    endGapPenalties = value;
}

void UserParameters::toggleResetAlignmentsNew()
{
    resetAlignmentsNew ^= true;
}

void UserParameters::setResetAlignmentsNew(bool value)
{
    resetAlignmentsNew = value;
}

void UserParameters::setResetAlignmentsAll(bool value)
{
    resetAlignmentsAll = value;
}

void UserParameters::setOutputStructPenalties(int value)
{
    outputStructPenalties = value;
}

void UserParameters::setStructPenalties1(int value)
{
    structPenalties1 = value;
}

void UserParameters::setStructPenalties2(int value)
{
    structPenalties2 = value;
}

void UserParameters::setUseSS1(bool value)
{
    useSS1 = value;
}

void UserParameters::toggleUseSS1()
{
    useSS1 ^= true;
}

void UserParameters::setUseSS2(bool value)
{
    useSS2 = value;
}

void UserParameters::toggleUseSS2()
{
    useSS2 ^= true;
}

void UserParameters::setHelixPenalty(int value)
{
    helixPenalty = value;
}

void UserParameters::setStrandPenalty(int value)
{
    strandPenalty = value;
}

void UserParameters::setLoopPenalty(int value)
{
    loopPenalty = value;
}

void UserParameters::setHelixEndMinus(int value)
{
    helixEndMinus = value;
}

void UserParameters::setHelixEndPlus(int value)
{
    helixEndPlus = value;
}

void UserParameters::setStrandEndMinus(int value)
{
    strandEndMinus = value;
}

void UserParameters::setStrandEndPlus(int value)
{
    strandEndPlus = value;
}

void UserParameters::setHelixEndPenalty(int value)
{
    helixEndPenalty = value;
}

void UserParameters::setStrandEndPenalty(int value)
{
    strandEndPenalty = value;
}

void UserParameters::setUseAmbiguities(bool value)
{
    useAmbiguities = value;
}

void UserParameters::setDNAPWGapOpenPenalty(float value)
{
    DNAPWGapOpen = value;
}

void UserParameters::setDNAPWGapExtendPenalty(float value)
{
    DNAPWGapExtend = value;
}

void UserParameters::setProteinPWGapOpenPenalty(float value)
{
    AAPWGapOpen = value;
}

void UserParameters::setProteinPWGapExtendPenalty(float value)
{
    AAPWGapExtend = value;
}

void UserParameters::toggleQuickPairAlign()
{
    quickPairAlign ^= true;
}

void UserParameters::setQuickPairAlign(bool value)
{
    quickPairAlign = value;
}

void UserParameters::setTransitionWeight(float value)
{
    transitionWeight = value;
}

void UserParameters::setDNAKtup(int value)
{
    DNAKtup = value;
}

void UserParameters::setDNAWindowGap(int value)
{
    DNAWindowGap = value;
}

void UserParameters::setDNASignif(int value)
{
    DNASignif = value;
}

void UserParameters::setDNAWindow(int value)
{
    DNAWindow = value;
}

void UserParameters::setAAKtup(int value)
{
    AAKtup = value;
}

void UserParameters::setAAWindowGap(int value)
{
    AAWindowGap = value;
}

void UserParameters::setAASignif(int value)
{
    AASignif = value;
}

void UserParameters::setAAWindow(int value)
{
    AAWindow = value;
}

void UserParameters::setPercent(bool value)
{
    percent = value;
}

void UserParameters::toggleTossGaps()
{
    tossgaps ^= true;
}

void UserParameters::setTossGaps(bool value)
{
    tossgaps = value;
}

void UserParameters::setKimura(bool value)
{
    kimura = value;
}

void UserParameters::toggleKimura()
{
    kimura ^= true;
}

void UserParameters::setBootNumTrials(int value)
{
    bootNumTrials = value;
}

void UserParameters::setBootRanSeed(unsigned int value)
{
    bootRanSeed = value;
}

void UserParameters::setDebug(int value)
{
    debug = value;
}

void UserParameters::setExplicitDNAFlag(bool value)
{
    explicitDNAFlag = value;
}

void UserParameters::setLowercase(bool value)
{
    lowercase = value;
}

void UserParameters::toggleLowercase()
{
    lowercase ^= true;
}

void UserParameters::setClSeqNumbers(bool value)
{
    clSeqNumbers = value;
}

void UserParameters::toggleClSeqNumbers()
{
    clSeqNumbers ^= true;
}

void UserParameters::setSeqRange(bool value)
{
    seqRange = value;
}

void UserParameters::toggleSeqRange()
{
    seqRange ^= true;
}

void UserParameters::setOutputClustal(bool value)
{
    outputClustal = value;
}

void UserParameters::toggleOutputClustal()
{
    outputClustal ^= true;
}

void UserParameters::setOutputGCG(bool value)
{
    outputGcg = value;
}

void UserParameters::toggleOutputGCG()
{
    outputGcg ^= true;
}

void UserParameters::setOutputPhylip(bool value)
{
    outputPhylip = value;
}

void UserParameters::toggleOutputPhylip()
{
    outputPhylip ^= true;
}

void UserParameters::setOutputNbrf(bool value)
{
    outputNbrf = value;
}

void UserParameters::toggleOutputNbrf()
{
    outputNbrf ^= true;
}

void UserParameters::setOutputGde(bool value)
{
    outputGde = value;
}

void UserParameters::toggleOutputGde()
{
    outputGde ^= true;
}

void UserParameters::setOutputNexus(bool value)
{
    outputNexus = value;
}

void UserParameters::toggleOutputNexus()
{
    outputNexus ^= true;
}

void UserParameters::setOutputFasta(bool value)
{
    outputFasta = value;
}

void UserParameters::toggleOutputFasta()
{
    outputFasta ^= true;
}

void UserParameters::setShowAlign(bool value)
{
    showAlign = value;
}

void UserParameters::toggleShowAlign()
{
    showAlign ^= true;
}

void UserParameters::setSaveParameters(bool value)
{
    saveParameters = value;
}

void UserParameters::toggleSaveParameters()
{
    saveParameters ^= true;
}

void UserParameters::setOutputTreeClustal(bool value)
{
    outputTreeClustal = value;
}

void UserParameters::toggleOutputTreeClustal()
{
    outputTreeClustal ^= true;
}

void UserParameters::setOutputTreePhylip(bool value)
{
    outputTreePhylip = value;
}

void UserParameters::toggleOutputTreePhylip()
{
    outputTreePhylip ^= true;
}

void UserParameters::setOutputTreeDistances(bool value)
{
    outputTreeDistances = value;
}

void UserParameters::toggleOutputTreeDistances()
{
    outputTreeDistances ^= true;
}

void UserParameters::setOutputTreeNexus(bool value)
{
    outputTreeNexus = value;
}

void UserParameters::toggleOutputTreeNexus()
{
    outputTreeNexus ^= true;
}

void UserParameters::setOutputPim(bool value)
{
    outputPim = value;
}

void UserParameters::setBootstrapFormat(int value)
{
    bootstrapFormat = value;
}

void UserParameters::toggleBootstrapFormat()
{
    if (bootstrapFormat == BS_NODE_LABELS)
    {
        bootstrapFormat = BS_BRANCH_LABELS;
    }
    else
    {
        bootstrapFormat = BS_NODE_LABELS;
    }
}

void UserParameters::setProfile1Name(string value)
{
    profile1Name = value;
}

void UserParameters::setProfile2Name(string value)
{
    profile2Name = value;
}

void UserParameters::setEmpty(bool value)
{
    empty = value;
}

void UserParameters::setProfile1Empty(bool value)
{
    profile1Empty = value;
}

void UserParameters::setProfile2Empty(bool value)
{
    profile2Empty = value;
}

void UserParameters::setOutfileName(string value)
{
    outfileName = value;
}

/*void UserParameters::setProfile1NumSeqs(int value)
{
    profile1NumSeqs = value;
}*/ // MARK CHANGE Jan 10

void UserParameters::setUseTreeFile(bool value)
{
    useTreeFile = value;
}

void UserParameters::setNewTreeFile(bool value)
{
    newTreeFile = value;
}

void UserParameters::setUseTree1File(bool value)
{
    useTree1File = value;
}

void UserParameters::setUseTree2File(bool value)
{
    useTree2File = value;
}

void UserParameters::setNewTree1File(bool value)
{
    newTree1File = value;
}

void UserParameters::setNewTree2File(bool value)
{
    newTree2File = value;
}

void UserParameters::setAminoAcidCodes(string value)
{
    aminoAcidCodes = value;
}

void UserParameters::setKtup(int value)
{
    ktup = value;
}

void UserParameters::setWindow(int value)
{
    window = value;
}

void UserParameters::setWindowGap(int value)
{
    windowGap = value;
}

void UserParameters::setSignif(int value)
{
    signif = value;
}

void UserParameters::setRangeFrom(int from)
{
    rangeFrom = from;
    rangeFromToSet = true;
}

void UserParameters::setRangeTo(int to)
{
    rangeTo = to;
    rangeFromToSet = true;
}

int UserParameters::getQTScorePlotScale()
{
    return QTscorePlotScale;
}

void UserParameters::setQTScorePlotScale(int score)
{
    QTscorePlotScale = score;
}

int UserParameters::getQTResExceptionCutOff()
{
    return QTresExceptionCutOff;
}

void UserParameters::setQTResExceptionCutOff(int cutOff)
{
    QTresExceptionCutOff = cutOff;
}

bool UserParameters::getQTseqWeightCalculated()
{
    return QTseqWeightCalculated;
}

void UserParameters::setQTseqWeightCalculated(bool calculated)
{
    QTseqWeightCalculated = calculated;
}

int UserParameters::getQTminLenLowScoreSegment()
{
    return QTminLenLowScoreSegment;
}

void UserParameters::setQTminLenLowScoreSegment(int minLen)
{
    QTminLenLowScoreSegment = minLen;
}

int UserParameters::getQTlowScoreDNAMarkingScale()
{
    return QTlowScoreDNAMarkingScale;
}

void UserParameters::setQTlowScoreDNAMarkingScale(int dnaScale)
{
    QTlowScoreDNAMarkingScale = dnaScale;
}

bool UserParameters::IterationIsEnabled()
{
    if (doRemoveFirstIteration == clustalw::NONE)
        return false;
    else
        return true;
}
}




