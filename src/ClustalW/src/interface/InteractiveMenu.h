/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Mark Larkin 12 Dec 2005.
 * Interactive menu class. It uses the object userParameters for
 * all the session variables.
 */
#ifndef INTERACTIVEMENU_H
#define INTERACTIVEMENU_H

#include <string>
#include "../Clustal.h"
#include "../general/clustalw.h"
#include "../general/userparams.h"
#include "../general/utils.h"
#include "../substitutionMatrix/globalmatrix.h"

namespace clustalw
{

using namespace std;

class InteractiveMenu
{
    public:
        /* Functions */
        InteractiveMenu(); 
        ~InteractiveMenu();
        void mainMenu();

        /* Attributes */

    private:
        /* Functions */
        void doSystem();
        void multipleAlignMenu();
        void profileAlignMenu();
        void ssOptionsMenu();
        int secStrOutputOptions();
        void phylogeneticTreeMenu();
        void treeFormatOptionsMenu();
        void formatOptionsMenu(); 
        void pairwiseMenu();
        void multiMenu();
        void gapPenaltiesMenu();
        int readMatrix(int alignResidueType, int alignType, MatMenu menu);
        void clusteringAlgorithmMenu();
        void iterationMenu();
         
        /* Attributes */
        Clustal* clustalObj;
        string phylipName;
        string clustalName;
        string distName;
        string nexusName;
        //string fasta_name;
        string p1TreeName;
        string p2TreeName;
        string secStructOutputTxt[4]; // Changed to a string array
        string lin1;
        MatMenu dnaMatrixMenu;
        MatMenu matrixMenu;
        MatMenu pwMatrixMenu;
        char choice;
};
}
#endif

