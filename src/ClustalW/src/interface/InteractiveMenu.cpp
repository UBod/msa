/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Mark Larkin Dec 12 2005.
 * This provides the implementation of the interactive menu functions.
 * Changes:
 * 15-5-07: Added changes to clustering algorithm choice in function phylogenticTreeMenu 
 *          Added iteration to the multipleAlignMenu function.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <stdio.h>
#include <string>
#include <ctype.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include "InteractiveMenu.h"

namespace clustalw
{
using namespace std;
 
InteractiveMenu::InteractiveMenu()
{
    try
    {
        clustalObj = new Clustal();
        lin1 = "";
        
        secStructOutputTxt[0] = string("Secondary Structure");
        secStructOutputTxt[1] = string("Gap Penalty Mask"); 
        secStructOutputTxt[2] = string("Structure and Penalty Mask"); 
        secStructOutputTxt[3] = string("None");
    }
    catch(bad_alloc)
    {
        cerr<<"The memory heap is exhausted. The program must terminate now\n";
        throw 1;
    }
    
    /* Initialise the menu structs */
    matrixMenu.noptions = 5;
    strcpy(matrixMenu.opt[AABLOSUM].title, "BLOSUM series");
    strcpy(matrixMenu.opt[AABLOSUM].string, "blosum");
    strcpy(matrixMenu.opt[AAPAM].title, "PAM series");
    strcpy(matrixMenu.opt[AAPAM].string, "pam");
    strcpy(matrixMenu.opt[AAGONNET].title, "Gonnet series");
    strcpy(matrixMenu.opt[AAGONNET].string, "gonnet");
    strcpy(matrixMenu.opt[AAIDENTITY].title, "Identity matrix");
    strcpy(matrixMenu.opt[AAIDENTITY].string, "id");
    strcpy(matrixMenu.opt[AAUSERDEFINED].title, "User defined");
    strcpy(matrixMenu.opt[AAUSERDEFINED].string, "");

    dnaMatrixMenu.noptions = 3;
    strcpy(dnaMatrixMenu.opt[DNAIUB].title, "IUB");
    strcpy(dnaMatrixMenu.opt[DNAIUB].string, "iub");
    strcpy(dnaMatrixMenu.opt[DNACLUSTALW].title, "CLUSTALW(1.6)");
    strcpy(dnaMatrixMenu.opt[DNACLUSTALW].string, "clustalw");
    strcpy(dnaMatrixMenu.opt[DNAUSERDEFINED].title, "User defined");
    strcpy(dnaMatrixMenu.opt[DNAUSERDEFINED].string, ""); 
    
    pwMatrixMenu.noptions = 5;
    strcpy(pwMatrixMenu.opt[PWAABLOSUM].title, "BLOSUM 30");
    strcpy(pwMatrixMenu.opt[PWAABLOSUM].string, "blosum");
    strcpy(pwMatrixMenu.opt[PWAAPAM].title, "PAM 350");
    strcpy(pwMatrixMenu.opt[PWAAPAM].string, "pam");
    strcpy(pwMatrixMenu.opt[PWAAGONNET].title, "Gonnet 250");
    strcpy(pwMatrixMenu.opt[PWAAGONNET].string, "gonnet");
    strcpy(pwMatrixMenu.opt[PWAAIDENTITY].title, "Identity matrix");
    strcpy(pwMatrixMenu.opt[PWAAIDENTITY].string, "id");
    strcpy(pwMatrixMenu.opt[PWAAUSER].title, "User defined");
    strcpy(pwMatrixMenu.opt[PWAAUSER].string, "");
    
}

InteractiveMenu::~InteractiveMenu()
{
    delete clustalObj;
}

void InteractiveMenu::mainMenu()
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<" **************************************************************\n";
        cout<<" ******** CLUSTAL "
            << userParameters->getRevisionLevel()
            <<" Multiple Sequence Alignments  ********\n";
        cout<<" **************************************************************\n";
        cout<<"\n\n";

        cout<<"     1. Sequence Input From Disc\n";
        cout<<"     2. Multiple Alignments\n";
        cout<<"     3. Profile / Structure Alignments\n";
        cout<<"     4. Phylogenetic trees\n\n";
        cout<<"     S. Execute a system command\n";
        cout<<"     H. HELP\n";
        cout<<"     X. EXIT (leave program)\n\n\n";

        choice = utilityObject->getChoice(string("Your choice"));
        
        string offendingSeq; // unused here
        switch (toupper(choice))
        {
            case '1':
                clustalObj->sequenceInput(false, &offendingSeq); 
                break;
            case '2':
                multipleAlignMenu();
                break;
            case '3':
                profileAlignMenu();
                break;
            case '4':
                phylogeneticTreeMenu();
                break;
            case 'S':
                doSystem(); 
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('1'); 
                break;
            case 'Q':
            case 'X':
            	throw 0;
                break;
            default:
                cout<<"\n\nUnrecognised Command\n\n";
                break;
        }
    }
}
void InteractiveMenu::multipleAlignMenu()
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<"****** MULTIPLE ALIGNMENT MENU ******\n\n\n";

        cout<<"    1.  Do complete multiple alignment now "
            <<(!userParameters->getQuickPairAlign() ? "Slow/Accurate" : "Fast/Approximate")
            <<"\n";
        cout<<"    2.  Produce guide tree file only\n";
        cout<<"    3.  Do alignment using old guide tree file\n\n";
        cout<<"    4.  Toggle Slow/Fast pairwise alignments = "
            <<((!userParameters->getQuickPairAlign()) ? "SLOW" : "FAST")
            <<"\n\n";
        cout<<"    5.  Pairwise alignment parameters\n";
        cout<<"    6.  Multiple alignment parameters\n\n";
        cout<<"    7.  Reset gaps before alignment?";

        if (userParameters->getResetAlignmentsNew())
        {
            cout<<" = ON\n";
        }
        else
        {
            cout<<" = OFF\n";
        }

        cout<<"    8.  Toggle screen display          = "
            <<((!userParameters->getShowAlign()) ? "OFF" : "ON")
            <<"\n";
        cout<<"    9.  Output format options\n";
        cout<<"    I. Iteration = ";
        
        if(userParameters->getDoRemoveFirstIteration() == ALIGNMENT)
        {
            cout << "ALIGNMENT\n\n";
        }
        else if(userParameters->getDoRemoveFirstIteration() == TREE)
        {
            cout << "TREE\n\n";
        }
        else
        {
            cout << "NONE\n\n";
        }
        
        cout<<"    S.  Execute a system command\n";
        cout<<"    H.  HELP\n";
        cout<<"    or press [RETURN] to go back to main menu\n\n\n";

        choice = utilityObject->getChoice(string("Your choice"));
        if (choice == '\n')
        {
            return ;
        }

        switch (toupper(choice))
        {
            case '1':
                clustalObj->align(&phylipName, output);
                break;
            case '2':
                clustalObj->doGuideTreeOnly(&phylipName); 
                break;
            case '3':
                clustalObj->doAlignUseOldTree(&phylipName);
                break;
            case '4':
                userParameters->toggleQuickPairAlign();
                break;
            case '5':
                pairwiseMenu();
                break;
            case '6':
                multiMenu();
                break;
            case '7':
                userParameters->toggleResetAlignmentsNew();
                if (userParameters->getResetAlignmentsNew() == true)
                {
                    userParameters->setResetAlignmentsAll(false);
                }
                break;
            case '8':
                userParameters->toggleShowAlign();
                break;
            case '9':
                formatOptionsMenu();
                break;
            case 'I':
                iterationMenu();
                break;    
            case 'S':
                doSystem();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('2');
                break;
            case 'Q':
            case 'X':
                return ;

            default:
                fprintf(stdout, "\n\nUnrecognised Command\n\n");
                break;
        }
    }
}

void InteractiveMenu::profileAlignMenu(void)
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<"****** PROFILE AND STRUCTURE ALIGNMENT MENU ******\n\n\n";

        cout<<"    1.  Input 1st. profile             ";

        if (!userParameters->getProfile1Empty())
        {
            cout<<"(loaded)";
        }

        cout<<"\n";
        cout<<"    2.  Input 2nd. profile/sequences   ";

        if (!userParameters->getProfile2Empty())
        {
            cout<<"(loaded)";
        }

        cout<<"\n\n";
        cout<<"    3.  Align 2nd. profile to 1st. profile\n";
        
        cout<<"    4.  Align sequences to 1st. profile "
            <<((!userParameters->getQuickPairAlign()) 
                ? "(Slow/Accurate)\n\n" : "(Fast/Approximate)\n\n");
                
        cout<<"    5.  Toggle Slow/Fast pairwise alignments = "
            <<((!userParameters->getQuickPairAlign()) 
                   ? "SLOW\n\n" : "FAST\n\n");
                   
        cout<<"    6.  Pairwise alignment parameters\n";
        cout<<"    7.  Multiple alignment parameters\n\n";
        
        cout<<"    8.  Toggle screen display                = "
            <<((!userParameters->getShowAlign()) ? "OFF\n" : "ON\n");
            
        cout<<"    9.  Output format options\n";
        cout<<"    0.  Secondary structure options\n\n";
        cout<<"    S.  Execute a system command\n";
        cout<<"    H.  HELP\n";
        cout<<"    or press [RETURN] to go back to main menu\n\n\n";

        choice = utilityObject->getChoice(string("Your choice"));
        if (choice == '\n')
        {
            return ;
        }

        switch (toupper(choice))
        {
            case '1':
                clustalObj->profile1Input();
                break;
            case '2':
                clustalObj->profile2Input();
                break;
            case '3':
                clustalObj->profileAlign(&p1TreeName, &p2TreeName); 
                break;
            case '4':
                /* align new sequences to profile 1 */
                clustalObj->sequencesAlignToProfile(&phylipName);
                break;
            case '5':
                userParameters->toggleQuickPairAlign();
                break;
            case '6':
                pairwiseMenu();
                break;
            case '7':
                multiMenu();
                break;
            case '8':
                userParameters->toggleShowAlign();
                break;
            case '9':
                formatOptionsMenu();
                break;
            case '0':
                ssOptionsMenu();
                break;
            case 'S':
                doSystem();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('6');
                break;
            case 'Q':
            case 'X':
                return ;

            default:
                cout<<"\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

void InteractiveMenu::ssOptionsMenu()
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<" ********* SECONDARY STRUCTURE OPTIONS *********\n\n\n";

        cout<<"     1. Use profile 1 secondary structure / penalty mask  ";       
        if (userParameters->getUseSS1())
        {
            cout<<"= YES\n";
        }
        else
        {
            cout<<"= NO\n";
        }
        
        cout<<"     2. Use profile 2 secondary structure / penalty mask  ";
        if (userParameters->getUseSS2())
        {
            cout<<"= YES\n\n";
        }
        else
        {
            cout<<"= NO\n\n";
        }

        cout<<"     3. Output in alignment  ";
        cout<<"= "<< secStructOutputTxt[userParameters->getOutputStructPenalties()] <<"\n\n";

        cout<<"     4. Helix gap penalty                     : "
            <<userParameters->getHelixPenalty()<<"\n";
        cout<<"     5. Strand gap penalty                    : "
            <<userParameters->getStrandPenalty()<<"\n";
        cout<<"     6. Loop gap penalty                      : "
            <<userParameters->getLoopPenalty()<<"\n";
        cout<<"     7. Secondary structure terminal penalty  : "
            <<userParameters->getHelixEndPenalty()<<"\n";
        cout<<"     8. Helix terminal positions       within : "
            <<userParameters->getHelixEndMinus()
            <<"      outside : "
            <<userParameters->getHelixEndPlus()<<"\n";
        cout<<"     9. Strand terminal positions      within : "
            <<userParameters->getStrandEndMinus()   
            <<"      outside : " 
            <<userParameters->getStrandEndPlus()<<"\n\n\n";
        cout<<"     H. HELP\n\n\n";

        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }
        
        switch (toupper(choice))
        {
            case '1':
                userParameters->toggleUseSS1();
                break;
            case '2':
                userParameters->toggleUseSS2();
                break;
            case '3':
                userParameters->setOutputStructPenalties(secStrOutputOptions());
                break;
            case '4': 
                cout<<"Helix Penalty Currently: "
                    << userParameters->getHelixPenalty()<<"\n";
                userParameters->setHelixPenalty(utilityObject->getInt("Enter number",
                                                  1,9, userParameters->getHelixPenalty()));
                break;
            case '5':
                cout<<"Strand Gap Penalty Currently: "
                    << userParameters->getStrandPenalty() <<"\n";
                userParameters->setStrandPenalty(utilityObject->getInt("Enter number",
                                                  1,9, userParameters->getStrandPenalty()));
                break;
            case '6':
                cout<<"Loop Gap Penalty Currently: "
                    << userParameters->getLoopPenalty() <<"\n";
                userParameters->setLoopPenalty(utilityObject->getInt("Enter number", 
                                                 1, 9, userParameters->getLoopPenalty()));
                break;
            case '7':
                cout<<"Secondary Structure Terminal Penalty Currently: "
                    << userParameters->getHelixEndPenalty() <<"\n";
                userParameters->setHelixEndPenalty(utilityObject->getInt("Enter number",
                                               1, 9,userParameters->getHelixEndPenalty()));
                userParameters->setStrandEndPenalty(userParameters->getHelixEndPenalty());
                break;
            case '8':
                cout<<"Helix Terminal Positions Currently: \n";
                cout<<"        within helix: "
                    << userParameters->getHelixEndMinus()
                    << "      outside helix: " 
                    << userParameters->getHelixEndPlus() <<"\n";
                    
                userParameters->setHelixEndMinus(utilityObject->getInt(
                                 "Enter number of residues within helix", 
                                 0, 3, userParameters->getHelixEndMinus()));
                userParameters->setHelixEndPlus(utilityObject->getInt(
                                 "Enter number of residues outside helix", 0, 3,
                                 userParameters->getHelixEndPlus()));
                break;
            case '9':
                cout<<"Strand Terminal Positions Currently: \n";
                cout<<"        within strand: "
                    << userParameters->getStrandEndMinus()
                    << "      outside strand: "
                    << userParameters->getStrandEndPlus() <<"\n";
                userParameters->setStrandEndMinus(utilityObject->getInt(
                                 "Enter number of residues within strand", 0, 3,
                                 userParameters->getStrandEndMinus()));
                userParameters->setStrandEndPlus(utilityObject->getInt(
                                 "Enter number of residues outside strand", 0, 3,
                                 userParameters->getStrandEndPlus()));
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('B');
                break;
            default:
                cout<<"\n\nUnrecognised Command\n\n";
                break;
        }
    }
}
int InteractiveMenu::secStrOutputOptions()
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<" ********* Secondary Structure Output Menu *********\n\n\n";

        cout<<"     1. "<< secStructOutputTxt[0] <<"\n";
        cout<<"     2. "<< secStructOutputTxt[1] <<"\n";
        cout<<"     3. "<< secStructOutputTxt[2] <<"\n";
        cout<<"     4. "<< secStructOutputTxt[3] <<"\n";
        cout<<"     H. HELP\n\n";
        cout<<"     -- Current output is "
            << secStructOutputTxt[userParameters->getOutputStructPenalties()];
        cout<<" --\n";

        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return (userParameters->getOutputStructPenalties());;
        }

        switch (toupper(choice))
        {
            case '1':
                return (OUTSECST);
            case '2':
                return (OUTGAP);
            case '3':
                return (OUTBOTH);
            case '4':
                return (OUTNONE);
            case '?':
            case 'H':
                clustalObj->getHelp('C');
            case 'Q':
            case 'X':
                return (0);

            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

void InteractiveMenu::phylogeneticTreeMenu()
{
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<"****** PHYLOGENETIC TREE MENU ******\n\n\n";

        cout<<"    1.  Input an alignment\n";
        cout<<"    2.  Exclude positions with gaps?        ";

        if (userParameters->getTossGaps())
        {
            cout<<"= ON\n";
        }
        else
        {
            cout<<"= OFF\n";
        }

        cout<<"    3.  Correct for multiple substitutions? ";

        if (userParameters->getKimura())
        {
            cout<<"= ON\n";
        }
        else
        {
            cout<<"= OFF\n";
        }

        cout<<"    4.  Draw tree now\n";
        cout<<"    5.  Bootstrap tree\n";
        cout<<"    6.  Output format options\n";
        cout<<"    7.  Clustering algorithm = "; // Mark change 15-5-2007
        if(userParameters->getClusterAlgorithm() == NJ)
        {
            cout << "NJ\n\n";
        }
        else
        {
            cout << "UPGMA\n\n";
        }        
        cout<<"    S.  Execute a system command\n";
        cout<<"    H.  HELP\n";
        cout<<"    or press [RETURN] to go back to main menu\n\n\n";

        choice = utilityObject->getChoice(string("Your choice"));
        if (choice == '\n')
        {
            return ;
        }
        
        string offendingSeq; // unused here
        switch (toupper(choice))
        {
            case '1':
                clustalObj->sequenceInput(false, &offendingSeq); 
                break;
            case '2':
                userParameters->toggleTossGaps();
                break;
            case '3':
                userParameters->toggleKimura();
                break;
            case '4':
                clustalObj->phylogeneticTree(&phylipName, &clustalName, &distName,
                    &nexusName, "amenu.pim");
                break;
            case '5':
                clustalObj->bootstrapTree(&phylipName, &clustalName, &nexusName);
                break;
            case '6':
                treeFormatOptionsMenu();
                break;
            case '7':
                clusteringAlgorithmMenu();
                break;
            case 'S':
                doSystem();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('7');
                break;
            case 'Q':
            case 'X':
                return ;

            default:
                cout<<"\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

void InteractiveMenu::treeFormatOptionsMenu()
{
    while (true)
    {
        lin1 = "";
        cout<< "\n\n\n";
        cout<<" ****** Format of Phylogenetic Tree Output ******\n\n\n";
        cout<<"     1. Toggle CLUSTAL format tree output    =  "
            << ((!userParameters->getOutputTreeClustal()) ? "OFF" : "ON")<<"\n";
        cout<<"     2. Toggle Phylip format tree output     =  "
            << ((!userParameters->getOutputTreePhylip()) ? "OFF" : "ON")<<"\n";
        cout<<"     3. Toggle Phylip distance matrix output =  "
            << ((!userParameters->getOutputTreeDistances()) ? "OFF" : "ON")<<"\n";
        cout<<"     4. Toggle Nexus format tree output      =  "
            << ((!userParameters->getOutputTreeNexus()) ? "OFF" : "ON")<<"\n\n";
        cout<<"     5. Toggle Phylip bootstrap positions    =  "
            <<((userParameters->getBootstrapFormat() == BS_NODE_LABELS) ? "NODE LABELS" :
            "BRANCH LABELS") <<"\n\n\n";
        cout<<"     H. HELP\n\n\n";
        
        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }
        switch (toupper(choice))
        {
            case '1':
                userParameters->toggleOutputTreeClustal();
                break;
            case '2':
                userParameters->toggleOutputTreePhylip();
                break;
            case '3':
                userParameters->toggleOutputTreeDistances();
                break;
            case '4':
                userParameters->toggleOutputTreeNexus();
                break;
            case '5':
                userParameters->toggleBootstrapFormat();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('0');
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }

}

void InteractiveMenu::formatOptionsMenu()
{
    
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<" ********* Format of Alignment Output *********\n\n\n";
        cout<<"     F. Toggle FASTA format output       =  "
            << ((!userParameters->getOutputFasta()) ? "OFF" : "ON") <<"\n\n";
        cout<<"     1. Toggle CLUSTAL format output     =  "
            << ((!userParameters->getOutputClustal()) ? "OFF" : "ON") <<"\n";
        cout<<"     2. Toggle NBRF/PIR format output    =  "
            << ((!userParameters->getOutputNbrf()) ? "OFF" : "ON") << "\n";
        cout<<"     3. Toggle GCG/MSF format output     =  "
            << ((!userParameters->getOutputGCG()) ? "OFF" : "ON") << "\n";
        cout<<"     4. Toggle PHYLIP format output      =  "
            << ((!userParameters->getOutputPhylip()) ? "OFF" : "ON") << "\n";
        cout<<"     5. Toggle NEXUS format output       =  "
            << ((!userParameters->getOutputNexus()) ? "OFF" : "ON") << "\n";
        cout<<"     6. Toggle GDE format output         =  "
            << ((!userParameters->getOutputGde()) ? "OFF" : "ON") << "\n\n";
        cout<<"     7. Toggle GDE output case           =  "
            << ((!userParameters->getLowercase()) ? "UPPER" : "LOWER") << "\n";

        cout<<"     8. Toggle CLUSTALW sequence numbers =  "
            << ((!userParameters->getClSeqNumbers()) ? "OFF" : "ON") << "\n";
        cout<<"     9. Toggle output order              =  "
            << ((userParameters->getOutputOrder() == 0) ? "INPUT FILE" : "ALIGNED")
            <<"\n\n";

        cout<<"     0. Create alignment output file(s) now?\n\n";
        cout<<"     T. Toggle parameter output          = "
            << ((!userParameters->getSaveParameters()) ? "OFF" : "ON") << "\n";
        cout<<"     R. Toggle sequence range numbers =  "
            <<((!userParameters->getSeqRange()) ? "OFF" : "ON") << "\n\n";

        cout<<"     H. HELP\n\n\n";

        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }
        switch (toupper(choice))
        {
            case '1':
                userParameters->toggleOutputClustal();
                break;
            case '2':
                userParameters->toggleOutputNbrf();
                break;
            case '3':
                userParameters->toggleOutputGCG();
                break;
            case '4':
                userParameters->toggleOutputPhylip();
                break;
            case '5':
                userParameters->toggleOutputNexus();
                break;
            case '6':
                userParameters->toggleOutputGde();
                break;
            case '7':
                userParameters->toggleLowercase();
                break;
            case '8':
                userParameters->toggleClSeqNumbers();
                break;
            case '9':
                userParameters->toggleOutputOrder();
                break;
            case 'F':
                userParameters->toggleOutputFasta();
                break;
            case 'R':
                userParameters->toggleSeqRange();
                break;
            case '0':
                clustalObj->outputNow();
                break;
            case 'T':
                userParameters->toggleSaveParameters();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('5');
                break;
            default:
                cout<<"\n\nUnrecognised Command\n\n";
                break;
        }
    }    
}
 
void InteractiveMenu::pairwiseMenu()
{
    if (userParameters->getDNAFlag())
    {
        userParameters->setPWParamToDNA();
    }
    else
    {
        userParameters->setPWParamToProtein();
    }

    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n";
        cout<<" ********* PAIRWISE ALIGNMENT PARAMETERS *********\n\n\n";

        cout<<"     Slow/Accurate alignments:\n\n";

        cout<<"     1. Gap Open Penalty       : "
            << fixed << setprecision(2) << userParameters->getPWGapOpen() << "\n";
        cout<<"     2. Gap Extension Penalty  : "
            << fixed << setprecision(2) << userParameters->getPWGapExtend() << "\n";
        cout<< "     3. Protein weight matrix  :"
            << matrixMenu.opt[subMatrix->getPWMatrixNum() - 1].title 
            << "\n";
        cout<<"     4. DNA weight matrix      :"
            << dnaMatrixMenu.opt[subMatrix->getPWDNAMatrixNum() - 1].title 
            << "\n\n";

        cout<<"     Fast/Approximate alignments:\n\n";

        cout<<"     5. Gap penalty            :"
            << userParameters->getWindowGap() << "\n";
        cout<<"     6. K-tuple (word) size    :"
            << userParameters->getKtup() << "\n";
        cout<<"     7. No. of top diagonals   :"
            << userParameters->getSignif() << "\n";
        cout<<"     8. Window size            :"
            << userParameters->getWindow() << "\n\n";

        cout<<"     9. Toggle Slow/Fast pairwise alignments ";
        if (userParameters->getQuickPairAlign())
        {
            cout<<"= FAST\n\n";
        }
        else
        {
            cout<<"= SLOW\n\n";
        }


        cout<<"     H. HELP\n\n\n";
        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            if (userParameters->getDNAFlag())
            {
                userParameters->setPWDNAParam();
            }
            else
            {
                userParameters->setPWProteinParam();
            }

            return ;
        }


        switch (toupper(choice))
        {
            case '1':
                cout<<"Gap Open Penalty Currently: "
                    << userParameters->getPWGapOpen() << "\n";
                userParameters->setPWGapOpen(
                              (float)utilityObject->getReal("Enter number", (double)0.0,
                              (double)100.0, (double)userParameters->getPWGapOpen()));
                break;
            case '2':
                cout<<"Gap Extension Penalty Currently: "
                    << userParameters->getPWGapExtend() << "\n";
                userParameters->setPWGapExtend(
                              (float)utilityObject->getReal("Enter number", (double)0.0,
                              (double)10.0, (double)userParameters->getPWGapExtend()));
                break;
            case '3':            
                readMatrix(Protein, Pairwise, pwMatrixMenu);
                break;
            case '4': 
                readMatrix(DNA, Pairwise, dnaMatrixMenu);
                break;
            case '5':
                cout<<"Gap Penalty Currently: "
                    << userParameters->getWindowGap() << "\n";
                userParameters->setWindowGap(
                                 utilityObject->getInt("Enter number", 1, 500, 
                                                userParameters->getWindowGap()));
                break;
            case '6':
                cout<<"K-tuple Currently: "
                    << userParameters->getKtup() << "\n";
                if (userParameters->getDNAFlag())
                {
                    int _ktup = utilityObject->getInt("Enter number", 1, 4,
                                                      userParameters->getKtup());
                    userParameters->setKtup(_ktup);
                    // see bug 185
                    userParameters->setDNAKtup(_ktup);
                    userParameters->setWindowGap(_ktup + 4);
                    userParameters->setDNAWindowGap(_ktup + 4);

                }
                else
                {
                    int _ktup = utilityObject->getInt("Enter number", 1, 2, 
                                                      userParameters->getKtup());
                    userParameters->setKtup(_ktup);
                    // see bug 185
                    userParameters->setAAKtup(_ktup);
                    userParameters->setWindowGap(_ktup + 3);
                    userParameters->setAAWindowGap(_ktup + 3);
                     
                }
                break;
            case '7':
                cout<<"Top diagonals Currently: "
                    << userParameters->getSignif() << "\n";
                userParameters->setSignif(
                                utilityObject->getInt("Enter number", 1, 50, 
                                          userParameters->getSignif()));
                break;
            case '8':
                cout<<"Window size Currently: "
                    << userParameters->getWindow() << "\n";
                userParameters->setWindow(
                                utilityObject->getInt("Enter number", 1, 50, 
                                           userParameters->getWindow()));
                break;
            case '9':
                userParameters->toggleQuickPairAlign();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('3');
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

void InteractiveMenu::multiMenu()
{
    if (userParameters->getDNAFlag())
    {
        userParameters->setDNAMultiGap();
    }
    else
    {
        userParameters->setProtMultiGap();
    }

    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n ********* MULTIPLE ALIGNMENT PARAMETERS *********\n\n\n";
        cout<<"     1. Gap Opening Penalty              :"
            << fixed << setprecision(2) << userParameters->getGapOpen() << "\n";
        cout<<"     2. Gap Extension Penalty            :"
            << fixed << setprecision(2) << userParameters->getGapExtend() << "\n";

        cout<<"     3. Delay divergent sequences        :"
            << userParameters->getDivergenceCutoff() << " %\n\n";

        cout<<"     4. DNA Transitions Weight           :"
            << fixed << setprecision(2) << userParameters->getTransitionWeight() << "\n\n";
        cout<<"     5. Protein weight matrix            :"
            << matrixMenu.opt[subMatrix->getMatrixNum() - 1].title 
            << "\n";
        cout<<"     6. DNA weight matrix                :"
            << dnaMatrixMenu.opt[subMatrix->getDNAMatrixNum() - 1].title 
            << "\n";
        cout<<"     7. Use negative matrix              :"
            << ((!userParameters->getUseNegMatrix()) ? "OFF" : "ON") << "\n\n";
        cout<<"     8. Protein Gap Parameters\n\n";
        cout<<"     H. HELP\n\n\n";

        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            if (userParameters->getDNAFlag())
            {
                //userParameters->setDNAMultiGap();
                userParameters->setDNAGapOpen(userParameters->getGapOpen());
                userParameters->setDNAGapExtend(userParameters->getGapExtend());
            }
            else
            {
                //userParameters->setProtMultiGap();
                userParameters->setAAGapOpen(userParameters->getGapOpen());
                userParameters->setAAGapExtend(userParameters->getGapExtend());
            }
            return ;
        }


        switch (toupper(choice))
        {
            case '1':
                cout<<"Gap Opening Penalty Currently: "
                    << userParameters->getGapOpen() << "\n";
                userParameters->setGapOpen(
                                  (float)utilityObject->getReal("Enter number", (double)0.0, 
                                  (double)100.0, (double)userParameters->getGapOpen()));
                break;
            case '2':
                cout<<"Gap Extension Penalty Currently: "
                    << userParameters->getGapExtend() << "\n";
                userParameters->setGapExtend(
                                (float)utilityObject->getReal("Enter number", (double)0.0,
                                (double)10.0, (double)userParameters->getGapExtend()));
                break;
            case '3':
                cout<<"Min Identity Currently: "
                    << userParameters->getDivergenceCutoff() << "\n";
                userParameters->setDivergenceCutoff(
                                     utilityObject->getInt("Enter number", 0, 100,
                                     userParameters->getDivergenceCutoff()));
                break;
            case '4':
                cout<<"Transition Weight Currently: "
                    << userParameters->getTransitionWeight() << "\n";
                userParameters->setTransitionWeight(
                                (float)utilityObject->getReal("Enter number", (double)0.0,
                                (double)1.0,
                                (double)userParameters->getTransitionWeight()));
                break;
            case '5':
                readMatrix(Protein, MultipleAlign, matrixMenu);
                break;
            case '6':
                readMatrix(DNA, MultipleAlign, dnaMatrixMenu);
                break;
            case '7':
                userParameters->toggleUseNegMatrix();
                break;
            case '8':
                gapPenaltiesMenu();
                break;
            case '?':
            case 'H':
                //clustalObj->getHelp('4');
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
}
 
void InteractiveMenu::gapPenaltiesMenu()
{
    
    while (true)
    {
        lin1 = "";
        cout<< "\n\n\n ********* PROTEIN GAP PARAMETERS *********\n\n\n\n";
        cout<<"     1. Toggle Residue-Specific Penalties :"
            << ((userParameters->getNoPrefPenalties()) ? "OFF" : "ON") << "\n\n";
        cout<<"     2. Toggle Hydrophilic Penalties      :"
            << ((userParameters->getNoHydPenalties()) ? "OFF" : "ON") << "\n";
        cout<<"     3. Hydrophilic Residues              :"
            << userParameters->getHydResidues() << "\n\n";
        cout<<"     4. Gap Separation Distance           :"
            << userParameters->getGapDist() << "\n";
        cout<<"     5. Toggle End Gap Separation         :"
            << ((!userParameters->getUseEndGaps()) ? "OFF" : "ON") << "\n\n";
        cout<<"     H. HELP\n\n\n";

        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }

        switch (toupper(choice))
        {
            case '1':
                userParameters->toggleNoPrefPenalties();
                break;
            case '2':
                userParameters->toggleNoHydPenalties();
                break;
            case '3':
                cout<<"Hydrophilic Residues Currently: "
                    << userParameters->getHydResidues() << "\n";
                lin1 = "";
                utilityObject->getStr(string("Enter residues (or [RETURN] to quit)"), lin1);
                if (lin1.size() > 0)
                {
                    userParameters->setHydResidues(lin1);
                }
                break;
            case '4':
                cout<<"Gap Separation Distance Currently: "
                    << userParameters->getGapDist() << "\n";
                userParameters->setGapDist(
                                utilityObject->getInt("Enter number", 0, 100, 
                                userParameters->getGapDist()));
                break;
            case '5':
                userParameters->toggleUseEndGaps();
                break;
            case '?':
            case 'H':
                clustalObj->getHelp('A');
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

/*
 * This function displays the menu for selecting the weight matrix to use.
 * It is used for both the protein and DNA menu for pairwise and Multiple alignment.
 * This is why it requires the MatMenu struct.
 */
int InteractiveMenu::readMatrix(int alignResidueType, int alignType, MatMenu menu)
{
    static char userFile[FILENAMELEN + 1];
    int i, option;
    char title[10];
    int matn; // Used to show which is the current matrix.
    
    if(alignResidueType == Protein)
    {
        strcpy(title, "PROTEIN");
    }
    else // DNA
    {
        strcpy(title, "DNA");
    }
    
    while (true)
    {
        lin1 = "";
        cout<<"\n\n\n ********* "<< title <<" WEIGHT MATRIX MENU *********\n\n\n";
        
        // Find out what the currently selected matrix is.
        matn = subMatrix->getMatrixNumForMenu(alignResidueType, alignType);
        
        for (i = 0; i < menu.noptions; i++)
        {
            cout<< "     " << i + 1 << ". " << menu.opt[i].title << "\n";
        }

        cout<<"     H. HELP\n\n";
        cout<<"     -- Current matrix is the "
            << menu.opt[matn -1].title << " ";
        
        if (matn == menu.noptions)
        {
            cout<<"(file = "<< userFile <<")";;
        }

        cout<<"--\n";
        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return matn;
        }

        option = toupper(choice) - '0';
        // Select the matrix series to be using 
        if (option > 0 && option < menu.noptions)
        {
            subMatrix->setCurrentNameAndNum(string(menu.opt[i - 1].string), option, 
                                            alignResidueType, alignType);
        }
        else if (option == menu.noptions) // Read in a User defined matrix.
        {
            // NOTE this will be changed to deal with matrix series.
            if (subMatrix->getUserMatFromFile(userFile, alignResidueType, alignType))
            {
                subMatrix->setCurrentNameAndNum(userFile, option, 
                                            alignResidueType, alignType);
            }
        }
        else
        switch (toupper(choice))
        {
            case '?':
            case 'H':
                clustalObj->getHelp('8');
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
} 

void InteractiveMenu::doSystem()
{
    lin1 = "";
    utilityObject->getStr(string("\n\nEnter system command"), lin1);

    if (lin1.size() > 0)
    {
        system(lin1.c_str());
    }
    cout<< "\n\n";
}

void InteractiveMenu::clusteringAlgorithmMenu()
{
    string currentAlgorithm = "";
    cout<<"****** CLUSTERING ALGORITHM MENU ******\n\n\n";
    while (true)
    {
        if(userParameters->getClusterAlgorithm() == NJ)
        {
            currentAlgorithm = "Neighbour Joining";
        }
        else
        {
            currentAlgorithm = "UPGMA";
        }
        lin1 = "";
        cout<< "\n\n\n";
        cout<<" ****** Clustering Algorithms ******\n\n\n";
        cout<<"     1. Neighbour Joining \n";
        cout<<"     2. UPGMA \n";
        cout << "-- Current algorithm is "<< currentAlgorithm << " --\n\n\n";

        
        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }
        switch (toupper(choice))
        {
            case '1':
                userParameters->setClusterAlgorithm(NJ);
                break;
            case '2':
                userParameters->setClusterAlgorithm(UPGMA);
                break;
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }    
}

void InteractiveMenu::iterationMenu()
{
    string currentIteration = "";
    cout<<"****** ITERATION MENU ******\n\n\n";
    while (true)
    {
        if(userParameters->getDoRemoveFirstIteration() == ALIGNMENT)
        {
            currentIteration = "Alignment iteration";
        }
        else if(userParameters->getDoRemoveFirstIteration() == TREE)
        {
            currentIteration = "Tree based iteration";
        }
        else
        {
            currentIteration = "None";
        }
        lin1 = "";
        cout<< "\n\n\n";
        cout<<" ****** Iteration Choices ******\n\n\n";
        cout<<"     1. Off \n";
        cout<<"     2. Tree based iteration (iterates each profile alignment step) \n";
        cout<<"     3. Alignment iteration (iterates final alignment only) \n\n";
        cout << "-- Current selection is "<< currentIteration << " --\n\n\n";

        
        choice = utilityObject->getChoice(string("Enter number (or [RETURN] to exit)"));
        if (choice == '\n')
        {
            return ;
        }
        switch (toupper(choice))
        {
            case '1':
                userParameters->setDoRemoveFirstIteration(NONE);
                break;
            case '2':
                userParameters->setDoRemoveFirstIteration(TREE);
                break;
            case '3':
                userParameters->setDoRemoveFirstIteration(ALIGNMENT);
                break;    
            default:
                cout<< "\n\nUnrecognised Command\n\n";
                break;
        }
    }
}

}

