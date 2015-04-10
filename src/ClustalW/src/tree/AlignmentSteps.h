/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The AlignmentSteps class is used to hold the progressive alignment steps that have
 * been calculated from the guide tree.
 *
 * Note: I have pushed an empty vector onto steps, so that the steps will match up
 * with the old sets array.
 *
 ***************************************************************************************/
 
#ifndef ALIGNMENTSTEPS_H
#define ALIGNMENTSTEPS_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;

namespace clustalw
{

class AlignmentSteps
{
    public:
        /* Functions */
        AlignmentSteps() : numSteps(0){steps.push_back(vector<int>());}; // Empty vector 
        void saveSet(int n, int *groups);
        void saveSet(vector<int>* groups);
        int getNumSteps();
        string getNextStep();
        void printAlignSteps();
        const vector<vector<int> >* getSteps(){return &steps;};
        vector<int>* getStep(int i){return &steps[i];};
        void clear();
        /* Attributes */

    private:
        /* Functions */
        
        /* Attributes */
        vector<vector<int> > steps;
        int numSteps;
};

}
#endif
