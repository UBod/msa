/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "AlignmentSteps.h"

namespace clustalw
{

void AlignmentSteps::saveSet(int n, int *groups)
{
    vector<int> tempVec;
    tempVec.resize(n + 1);
    tempVec[0] = 0;
    for(int i = 1; i < n + 1; i++)
    {
        tempVec[i] = groups[i - 1];
    } 
    steps.push_back(tempVec);
    numSteps++;
}

void AlignmentSteps::saveSet(vector<int>* groups)
{
    steps.push_back(*groups);
    numSteps++;
}

int AlignmentSteps::getNumSteps()
{
    return numSteps;
}

void AlignmentSteps::printAlignSteps()
{
    int rows = steps.size();
    for(int i = 1; i < rows; i++)
    {
        for(int j = 1; j < (int)steps[i].size(); j++)
        {
            cout << " " << steps[i][j];
        }
        cout << "\n";
    }
    cout << "\n\n";
}

void AlignmentSteps::clear()
{
    int size = steps.size();
    for(int i = 0; i < size; i++)
    {
        steps[i].clear();
    }
    steps.clear();
    steps.push_back(vector<int>());
    numSteps = 0;
}

}
