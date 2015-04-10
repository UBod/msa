/**
 * Author: Andreas Wilm
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
/**
 * This is the clustalw help class which replaces the old help file
 *
 */
#ifndef HELP_H
#define HELP_H



#include <string>
#include <iostream>
using namespace std;

typedef struct {
    string marker;
    string title;
    string content;
} section;

class Help {
    
public:
    /* Functions */
    Help();
    ~Help();
    string GetSection(string marker);
    string GetSection(char marker);
    string GetSectionTitle(string marker);
    string GetSectionTitle(char marker);
    vector<string> ListSectionMarkers();
    /* Attributes */
    
private:
    /* Functions */
    /* Attributes */
    vector<section> sections;
};

#endif
