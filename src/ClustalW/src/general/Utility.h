/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef UTILITY_H
#define UTILITY_H

#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <ctype.h>
#include "clustalw.h"
class QLabel;
using namespace std;

namespace clustalw
{
class Utility
{
    public:
        Utility();
	virtual ~Utility(){};

        char* rTrim(char *str);
        void rTrim(string *str);
        char* blankToUnderscore(char *str);
        string blankToUnderscore(string str);
        void getStr(string instr, string& outstr);
        char getChoice(string instr);
        double getReal(const char *instr, double minx, double maxx, double def);
        int getInt(const char *instr,int minx,int maxx, int def);
        unsigned long getUniqueSequenceIdentifier();
        bool lineType(char *line, const char *code);

        bool blankLine(char *line);
        bool blankLineNumericLabel(char *line);
        
        void getPath(string str, string *path);
        
        virtual char promptForYesNo(char *title, const char *prompt);
        virtual char promptForYesNo(const char *title, const char *prompt);
        virtual void error( char *msg,...);
        virtual void error( const char *msg,...);
        virtual void warning( char *msg,...);
        virtual void warning( const char *msg,...);
        virtual void info( char *msg,...);
        virtual void info( const char *msg,...);
        virtual void myname( char *myname );
        template <class T> T MIN(T x, T y){if (x <= y){return x;}return y;}
        template <class T> T MAX(T x, T y){if (x >= y){return x;}return y;}
        // Note the following function will never be used from this class. It is for the
        // QT version. I need it here so that I can add the function to QTUtility.
        virtual void setInfoLabelPtr(QLabel* ptrToLabelObj);
        bool isNumeric(char ch);
        double average(std::vector<double>& v);
        double stdDev(std::vector<double>& v);
        double median(std::vector<double> v);
        /* Attributes */
        void beQuiet(bool b) {quiet=b;};
        
    private:
        /* Functions */


        /* Attributes */
        bool quiet;

};
}
#endif

