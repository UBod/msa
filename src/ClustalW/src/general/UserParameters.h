/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Mark Larkin 9 Dec 2005. 
 * Note that most of the variables in this class are from param.h.
 *
 * Change log:
 * Jan 18th 2006: Changed all of the c style strings into c++ string objects!
 * Mark 30-1-2007: Added useLEScoringFunction and access functions.
 * 31-01-07,Nigel Brown(EMBL): Revision string now obtained from clustalw
 * specific versionw.h and updated from 1.83 to 2.1.
 */
 
#ifndef USERPARAMETERS_H
#define USERPARAMETERS_H

#include <string>
#include <iostream>
#include "utils.h"
#include "general.h"
#include "clustalw.h"

namespace clustalw
{

class UserParameters
{
    public:
        /* Functions */
        UserParameters(bool log = false);
        void setParamsToDefault();
        void createParameterOutput();
        int resIndex(string t,char c);
        void setDNAMultiGap();
        void setProtMultiGap();
        void setDNAParams();
        void setProtParams();
        void setPWProteinParam();
        void setPWDNAParam();
        void setPWParamToProtein();
        void setPWParamToDNA();
        string getRevisionLevel();
        void setRevisionLevel(string value);
        bool getInteractive(){return interactive;};
        void setInteractive(bool value);
        bool getGui(){return gui;};
        void setGui(bool value);
        float getGapOpen(){return gapOpen;};
        void setGapOpen(float value);
        float getGapExtend(){return gapExtend;};
        void setGapExtend(float value);
        float getPWGapOpen(){return PWGapOpen;};
        void setPWGapOpen(float value);
        float getPWGapExtend(){return PWGapExtend;};
        void setPWGapExtend(float value);
        float getAAGapOpen(){return AAGapOpen;}
        void setAAGapOpen(float gap){AAGapOpen = gap;}
        float getAAGapExtend(){return AAGapExtend;}
        void setAAGapExtend(float gap){AAGapExtend = gap;}
        float getAAPWGapOpen(){return AAPWGapOpen;}
        void setAAPWGapOpen(float gap){AAPWGapOpen = gap;}
        float getAAPWGapExtend(){return AAPWGapExtend;}
        void setAAPWGapExtend(float gap){AAPWGapExtend = gap;}
        int getMaxAA(){return maxAA;};
        void setMaxAA(int value);
        int getGapPos1(){return gapPos1;};
        void setGapPos1(int value);
        int getGapPos2(){return gapPos2;};
        void setGapPos2(int value);
        int getProfileNum(){return profileNum;};
        void setProfileNum(int value);
        bool getMenuFlag(){return menuFlag;};
        void setMenuFlag(bool value);
        bool getDNAFlag(){return DNAFlag;};
        void setDNAFlag(bool value);
        bool getDistanceTree(){return distanceTree;};
        void setDistanceTree(bool value);
        string getSeqName(){return seqName;};
        void setSeqName(string value);
        float getDNAGapOpen(){return DNAGapOpen;};
        void setDNAGapOpen(float value);
        float getDNAGapExtend(){return DNAGapExtend;};
        void setDNAGapExtend(float value);
        float getProteinGapOpen(){return AAGapOpen;};
        void setProteinGapOpen(float value);
        float getProteinGapExtend(){return AAGapExtend;};
        void setProteinGapExtend(float value);
        int getGapDist(){return gapDist;};
        void setGapDist(int value);
        int getOutputOrder(){return outputOrder;};
        void setOutputOrder(int value);
        void toggleOutputOrder();
        int getDivergenceCutoff(){return divergenceCutoff;};
        void setDivergenceCutoff(int value);


        string getHydResidues(){return hydResidues;};
        void setHydResidues(string value);
        bool getNoWeights(){return noWeights;};
        void setNoWeights(bool value);
        bool getUseNegMatrix(){return negMatrix;};
        void setUseNegMatrix(bool value);
        void toggleUseNegMatrix();
        bool getNoHydPenalties(){return noHydPenalties;};
        void setNoHydPenalties(bool value);
        void toggleNoHydPenalties();
        bool getNoVarPenalties(){return noVarPenalties;};
        void setNoVarPenalties(bool value);
        bool getNoPrefPenalties(){return noPrefPenalties;};
        void setNoPrefPenalties(bool value);
        void toggleNoPrefPenalties();
        bool getUseEndGaps(){return useEndGaps;};
        void setUseEndGaps(bool value);
        void toggleUseEndGaps();
        bool getEndGapPenalties(){return endGapPenalties;};
        void setEndGapPenalties(bool value);
        bool getResetAlignmentsNew(){return resetAlignmentsNew;};
        void setResetAlignmentsNew(bool value);
        bool getResetAlignmentsAll(){return resetAlignmentsAll;};
        void toggleResetAlignmentsNew();
        void setResetAlignmentsAll(bool value);
        int getOutputStructPenalties(){return outputStructPenalties;};
        void setOutputStructPenalties(int value);
        int getStructPenalties1(){return structPenalties1;};
        void setStructPenalties1(int value);
        int getStructPenalties2(){return structPenalties2;};
        void setStructPenalties2(int value);
        bool getUseSS1(){return useSS1;};
        void setUseSS1(bool value);
        void toggleUseSS1();
        bool getUseSS2(){return useSS2;};
        void setUseSS2(bool value);
        void toggleUseSS2();
        int getHelixPenalty(){return helixPenalty;};
        void setHelixPenalty(int value);
        int getStrandPenalty(){return strandPenalty;};
        void setStrandPenalty(int value);
        int getLoopPenalty(){return loopPenalty;};
        void setLoopPenalty(int value);
        int getHelixEndMinus(){return helixEndMinus;};
        void setHelixEndMinus(int value);
        int getHelixEndPlus(){return helixEndPlus;};
        void setHelixEndPlus(int value);
        int getStrandEndMinus(){return strandEndMinus;};
        void setStrandEndMinus(int value);
        int getStrandEndPlus(){return strandEndPlus;};
        void setStrandEndPlus(int value);
        int getHelixEndPenalty(){return helixEndPenalty;};
        void setHelixEndPenalty(int value);
        int getStrandEndPenalty(){return strandEndPenalty;};
        void setStrandEndPenalty(int value);
        bool getUseAmbiguities(){return useAmbiguities;};
        void setUseAmbiguities(bool value);
        float getDNAPWGapOpenPenalty(){return DNAPWGapOpen;};
        void setDNAPWGapOpenPenalty(float value);
        float getDNAPWGapExtendPenalty(){return DNAPWGapExtend;};
        void setDNAPWGapExtendPenalty(float value);
        float getProteinPWGapOpenPenalty(){return AAPWGapOpen;};
        void setProteinPWGapOpenPenalty(float value);
        float getProteinPWGapExtendPenalty(){return AAPWGapExtend;};
        void setProteinPWGapExtendPenalty(float value);

        bool getQuickPairAlign(){return quickPairAlign;};
        void setQuickPairAlign(bool value);
        void toggleQuickPairAlign(); // Mark new!!!
        float getTransitionWeight(){return transitionWeight;};
        void setTransitionWeight(float value);
        int getDNAKtup(){return DNAKtup;};
        void setDNAKtup(int value);
        int getDNAWindowGap(){return DNAWindowGap;};
        void setDNAWindowGap(int value);
        int getDNASignif(){return DNASignif;};
        void setDNASignif(int value);
        int getDNAWindow(){return DNAWindow;};
        void setDNAWindow(int value);
        int getAAKtup(){return AAKtup;};
        void setAAKtup(int value);
        int getAAWindowGap(){return AAWindowGap;};
        void setAAWindowGap(int value);
        int getAASignif(){return AASignif;};
        void setAASignif(int value);
        int getAAWindow(){return AAWindow;};
        void setAAWindow(int value);
        bool getPercent(){return percent;};
        void setPercent(bool value);
        bool getTossGaps(){return tossgaps;};
        void setTossGaps(bool value);
        void toggleTossGaps();
        bool getKimura(){return kimura;};
        void setKimura(bool value);
        void toggleKimura();
        int getBootNumTrials(){return bootNumTrials;};
        void setBootNumTrials(int value);
        unsigned int getBootRanSeed(){return bootRanSeed;};
        void setBootRanSeed(unsigned int value);
        int getDebug(){return debug;};
        void setDebug(int value);
        bool getExplicitDNAFlag(){return explicitDNAFlag;};
        void setExplicitDNAFlag(bool value);
        bool getLowercase(){return lowercase;};
        void setLowercase(bool value);
        void toggleLowercase();
        bool getClSeqNumbers(){return clSeqNumbers;};
        void setClSeqNumbers(bool value);
        void toggleClSeqNumbers();
        bool getSeqRange(){return seqRange;};
        void setSeqRange(bool value);
        void toggleSeqRange();
        bool getOutputClustal(){return outputClustal;};
        void setOutputClustal(bool value);
        void toggleOutputClustal(); 
        bool getOutputGCG(){return outputGcg;};
        void setOutputGCG(bool value);
        void toggleOutputGCG(); 
        bool getOutputPhylip(){return outputPhylip;};
        void setOutputPhylip(bool value);
        void toggleOutputPhylip(); 
        bool getOutputNbrf(){return outputNbrf;};
        void setOutputNbrf(bool value);
        void toggleOutputNbrf(); 
        bool getOutputGde(){return outputGde;};
        void setOutputGde(bool value);
        void toggleOutputGde(); 
        bool getOutputNexus(){return outputNexus;};
        void setOutputNexus(bool value);
        void toggleOutputNexus(); 
        bool getOutputFasta(){return outputFasta;};
        void setOutputFasta(bool value);
        void toggleOutputFasta(); 
        bool getShowAlign(){return showAlign;};
        void setShowAlign(bool value);
        void toggleShowAlign();
        bool getSaveParameters(){return saveParameters;};
        void setSaveParameters(bool value);
        void toggleSaveParameters();
        bool getOutputTreeClustal(){return outputTreeClustal;};
        void setOutputTreeClustal(bool value);
        void toggleOutputTreeClustal();
        bool getOutputTreePhylip(){return outputTreePhylip;};
        void setOutputTreePhylip(bool value);
        void toggleOutputTreePhylip();
        bool getOutputTreeDistances(){return outputTreeDistances;};
        void setOutputTreeDistances(bool value);
        void toggleOutputTreeDistances();
        bool getOutputTreeNexus(){return outputTreeNexus;};
        void setOutputTreeNexus(bool value);
        void toggleOutputTreeNexus();
        bool getOutputPim(){return outputPim;};
        void setOutputPim(bool value);
        int getBootstrapFormat(){return bootstrapFormat;};
        void setBootstrapFormat(int value);
        void toggleBootstrapFormat();
        string getProfile1Name(){return profile1Name;};
        void setProfile1Name(string value);
        string getProfile2Name(){return profile2Name;};
        void setProfile2Name(string value);
        bool getEmpty(){return empty;};
        void setEmpty(bool value);
        bool getProfile1Empty(){return profile1Empty;};
        void setProfile1Empty(bool value);
        bool getProfile2Empty(){return profile2Empty;};
        void setProfile2Empty(bool value);
        string getOutfileName(){return outfileName;};
        void setOutfileName(string value);
        bool getUseTreeFile(){return useTreeFile;};
        void setUseTreeFile(bool value);
        bool getNewTreeFile(){return newTreeFile;};
        void setNewTreeFile(bool value);
        bool getUseTree1File(){return useTree1File;};
        void setUseTree1File(bool value);
        bool getUseTree2File(){return useTree2File;};
        void setUseTree2File(bool value);
        bool getNewTree1File(){return newTree1File;};
        void setNewTree1File(bool value);
        bool getNewTree2File(){return newTree2File;};
        void setNewTree2File(bool value);
        string getAminoAcidCodes(){return aminoAcidCodes;};
        char getAminoAcidCode(int i){return aminoAcidCodes[i];};
        void setAminoAcidCodes(string value);
        int getKtup(){return ktup;};
        void setKtup(int value);
        int getWindow(){return window;};
        void setWindow(int value);
        int getWindowGap(){return windowGap;};
        void setWindowGap(int value);
        int getSignif(){return signif;};
        void setSignif(int value);
        int getRangeFrom(){return rangeFrom;};
        int getRangeTo(){return rangeTo;};
        void setRangeFrom(int from);
        void setRangeTo(int to);
        bool getRangeFromToSet(){return rangeFromToSet;};
        void setRangeFromToSet(bool set){rangeFromToSet = set;};
        int getQTScorePlotScale();
        void setQTScorePlotScale(int score);
        int getQTResExceptionCutOff();
        void setQTResExceptionCutOff(int cutOff);
        bool getQTseqWeightCalculated();
        void setQTseqWeightCalculated(bool calculated);
        int getQTminLenLowScoreSegment();
        void setQTminLenLowScoreSegment(int minLen);
        int getQTlowScoreDNAMarkingScale();
        void setQTlowScoreDNAMarkingScale(int dnaScale);

        
        // Access functions for the iteration variables.
        void setNumIterations(int num){numIterations = num;}
        int getNumIterations(){return numIterations;}
        void setDoRemoveFirstIteration(int doIter){doRemoveFirstIteration = doIter;}
        int getDoRemoveFirstIteration(){return doRemoveFirstIteration;}
        bool IterationIsEnabled();
            
        void setClusterAlgorithm(int clust){clusterAlgorithm = clust;}
        int getClusterAlgorithm(){return clusterAlgorithm;}
        
        void setDisplayInfo(bool display){displayInfo = display;}
        bool getDisplayInfo(){return displayInfo;}
        bool getHelpFlag() {return helpFlag;}
        void setHelpFlag(bool b) {helpFlag = b;}
        bool getFullHelpFlag() {return fullHelpFlag;}
        void setFullHelpFlag(bool b) {fullHelpFlag = b;}
        void setMaxAllowedSeqLength(int num){maxAllowedSeqLength = num;}
        int getMaxAllowedSeqLength(){return maxAllowedSeqLength;}

        bool ResetGapsIsEnabled() {return (resetAlignmentsNew || resetAlignmentsAll);};

        /* Attributes */

    private:
        /* Functions */


        /* Attributes */
        
        string revisionLevel;
        bool interactive;
        bool gui;
        float gapOpen;
        float gapExtend;
        float PWGapOpen;
        float PWGapExtend;
        int maxAA;
        int gapPos1;
        int gapPos2;
        int profileNum;
        bool menuFlag;
        bool DNAFlag;
        bool distanceTree;
        string seqName;
        float DNAGapOpen;
        float DNAGapExtend;
        float AAGapOpen;
        float AAGapExtend;
        int gapDist;
        int outputOrder;
        int divergenceCutoff;
        string hydResidues;
        bool noWeights;
        bool negMatrix;
        bool noHydPenalties;
        bool noVarPenalties;
        bool noPrefPenalties;
        bool useEndGaps;
        bool endGapPenalties;
        bool resetAlignmentsNew;
        bool resetAlignmentsAll;
        int outputStructPenalties;
        int structPenalties1;
        int structPenalties2;
        bool useSS1;
        bool useSS2;
        int helixPenalty;
        int strandPenalty;
        int loopPenalty;
        int helixEndMinus;
        int helixEndPlus;
        int strandEndMinus;
        int strandEndPlus;
        int helixEndPenalty;
        int strandEndPenalty;
        bool useAmbiguities;
        float DNAPWGapOpen;
        float DNAPWGapExtend;
        float AAPWGapOpen;
        float AAPWGapExtend;

        bool quickPairAlign;
        float transitionWeight;
        int DNAKtup;
        int DNAWindowGap;
        int DNASignif;
        int DNAWindow;
        int AAKtup;
        int AAWindowGap;
        int AASignif;
        int AAWindow;
        bool percent;
        bool tossgaps;
        bool kimura;
        int bootNumTrials;
        unsigned int bootRanSeed;
        int debug;
        bool explicitDNAFlag;
        bool lowercase; /* Flag for GDE output - set on comm. line*/
        bool clSeqNumbers;
        bool seqRange;
        bool outputClustal;
        bool outputGcg;
        bool outputPhylip;
        bool outputNbrf;
        bool outputGde;
        bool outputNexus;
        bool outputFasta;
        bool showAlign;
        bool saveParameters;
        bool outputTreeClustal;
        bool outputTreePhylip;
        bool outputTreeDistances;
        bool outputTreeNexus;
        bool outputPim;
        int bootstrapFormat;
        string profile1Name;
        string profile2Name;
        bool empty;
        bool profile1Empty;
        bool profile2Empty;
        string outfileName;
        //int profile1NumSeqs; MARK CHANGE Jan 10
        bool useTreeFile;
        bool newTreeFile;
        bool useTree1File;
        bool useTree2File;
        bool newTree1File;
        bool newTree2File;
        string aminoAcidCodes;
        int ktup;
        int window;
        int windowGap;
        int signif;
        
        int rangeFrom;
        int rangeTo;
        bool rangeFromToSet;
        
        int QTscorePlotScale;
        int QTresExceptionCutOff;
        bool QTseqWeightCalculated;
        int QTminLenLowScoreSegment;
        int QTlowScoreDNAMarkingScale;

        
        // New variables for iteration
        int numIterations;
        int doRemoveFirstIteration;
        
        int clusterAlgorithm;
        bool displayInfo;
        bool helpFlag;
        bool fullHelpFlag;
        bool quiet;
        
        int maxAllowedSeqLength;
};
}
#endif

