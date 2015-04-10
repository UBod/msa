/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 *
 * @author Mark Larkin, Conway Institute, UCD. mark.larkin@ucd.ie
 *
 * Changes:
 *
 *  2007-12-03, Andreas Wilm (UCD): replaced gets with fgets, and
 *  got rid of some compiler warning 
 *
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream> 
#include <algorithm> 
#include "Utility.h"
#include "general/userparams.h"

namespace clustalw
{

/**   constructor
 */
Utility::Utility()
{
    quiet=false;
}

    
/**
 * Removes trailing blanks from a string
 */
void Utility::rTrim(string *str)
{
    string::reverse_iterator rit = str->rbegin();
     
     while(rit != str->rend() && isspace(*rit))
     {
         str->erase((++rit).base());
     }
}

/**
 * Removes trailing blanks from a string
 * @param str String to remove trailing blanks from.
 * @return Pointer to the processed string
 */
char * Utility::rTrim(char *str)
{
    register int p;

    p = strlen(str) - 1;

    while (isspace(str[p]))
    {
        p--;
    }

    str[p + 1] = EOS;

    return str;
}

/**
 * Replace blanks in a string with underscores. Also replaces , ; : ( or ) with _.
 * @param str String to process.
 * @return Pointer to the processed string 
 */
char * Utility::blankToUnderscore(char *str)
{
    int i, p;

    p = strlen(str) - 1;

    for (i = 0; i <= p; i++)
        if ((str[i] == ' ') || (str[i] == ';') || (str[i] == ',') || (str[i] ==
            '(') || (str[i] == ')') || (str[i] == ':'))
        {
            str[i] = '_';
        }

    return str;
}

/**
 * Replace blanks in a string with underscores. Also replaces , ; : ( or ) with _.
 * @param str String to process.
 * @return Pointer to the processed string 
 */
string Utility::blankToUnderscore(string str)
{
    int i, p;

    p = str.size() - 1;

    for (i = 0; i <= p; i++)
        if ((str[i] == ' ') || (str[i] == ';') || (str[i] == ',') || (str[i] ==
            '(') || (str[i] == ')') || (str[i] == ':'))
        {
            str[i] = '_';
        }

    return str;
}

/**
 * 
 * @param instr 
 * @return 
 */
char Utility::getChoice(string instr)
{
    cout << instr << ": ";
    cout.flush();
    char choice;
    cin.get(choice);
    // We only want one character, so we ignore the rest.
    if(choice != '\n')
    {
        cin.ignore(150, '\n');
    }
    cin.clear();
    if(isalpha(choice) || isNumeric(choice))
    {
        return choice;
    }
    else if(choice == '\n')
    {
        return '\n';
    }
    else
    {
        return ' ';
    }
}

bool Utility::isNumeric(char ch)
{
    if(ch >= '0' && ch <= '9')
    {
        return true;
    }
    return false;
}


/**
 * 
 * @param instr 
 * @param outstr 
 */
void Utility::getStr(string instr, string& outstr)
{
    cout << instr << ": ";
    cout.flush();    
    string temp;
    getline(cin, temp, '\n');
    outstr = temp;
    cin.clear();
}

/**
 * 
 * @param instr 
 * @param minx 
 * @param maxx 
 * @param def 
 * @return 
 */
double Utility::getReal(const char *instr, double minx, double maxx, double def)
{
    int status;
    float ret;
    char line[MAXLINE];

    while (true)
    {
        fprintf(stdout, "%s (%.1f-%.1f)   [%.1f]: ", instr, minx, maxx, def);
        //gets(line);
        fgets(line, MAXLINE, stdin);
        status = sscanf(line, "%f", &ret);
        if (status == EOF)
        {
            return def;
        }
        if (ret > maxx)
        {
            fprintf(stderr, "ERROR: Max. value=%.1f\n\n", maxx);
            continue;
        }
        if (ret < minx)
        {
            fprintf(stderr, "ERROR: Min. value=%.1f\n\n", minx);
            continue;
        }
        break;
    }
    return (double)ret;
}

/**
 * 
 * @param instr 
 * @param minx 
 * @param maxx 
 * @param def 
 * @return 
 */
int Utility::getInt(const char *instr, int minx, int maxx, int def)
{
    int ret, status;
    char line[MAXLINE];

    while (true)
    {
        fprintf(stdout, "%s (%d..%d)    [%d]: ", instr, minx, maxx, def);
        //gets(line);
        fgets(line, MAXLINE, stdin);
        status = sscanf(line, "%d", &ret);
        if (status == EOF)
        {
            return def;
        }
        if (ret > maxx)
        {
            fprintf(stderr, "ERROR: Max. value=%d\n\n", maxx);
            continue;
        }
        if (ret < minx)
        {
            fprintf(stderr, "ERROR: Min. value=%d\n\n", minx);
            continue;
        }
        break;
    }
    return ret;
}

/**
 * 
 * @param line 
 * @param code 
 * @return 
 */
bool Utility::lineType(char *line, const char *code)
{
   // AW: introduced sanity check
   int n;
   if (strlen(line)<strlen(code))
     n=strlen(line);
   else
     n=strlen(code);
   
   return (strncmp(line, code, strlen(code)) == 0);
}

/**
 * 
 * @param line 
 * @return 
 */
bool Utility::blankLine(char *line)
{
    int i;

    for (i = 0; line[i] != '\n' && line[i] != EOS; i++)
    {
        if (isdigit(line[i]) || isspace(line[i]) || (line[i] == '*') ||
            (line[i] == ':') || (line[i] == '.'))
            ;
        else
        {
            return false;
        }
    }
    return true;
}

/**
 * 
 * @param line 
 * @return 
 */
// Utility::blankLine thinks that a line like
// 2125209         .......... .......... .......... .......... .......... 
// is blank, causing problems for countSeqs when sequence labels are numeric
// this function will return false if there's a number and a bunch of dots (not ideal, i know)
bool Utility::blankLineNumericLabel(char *line)
{
    int i;
    int dots = 0;
    bool isnumeric = false;

    for (i = 0; line[i] != '\n' && line[i] != EOS; i++)
    {
        if (isdigit(line[i]) || isspace(line[i]) || (line[i] == '*') ||
            (line[i] == ':') || (line[i] == '.'))
            ;
        else
        {
            return false;
        }
        if(line[i] == '.')
            dots++;
        if(isdigit(line[i]))
            isnumeric = true;
    }
    if(isnumeric && dots > 10)
        return false;
    else
        return true;
}

/**
 * 
 * @param str 
 * @param path 
 */
void Utility::getPath(string str, string *path)
{
    int i;
    string _temp;
    _temp = str;
    
    for (i = _temp.length() - 1; i >  - 1; --i)
    {
        if (str[i] == DIRDELIM)
        {
            i =  - 1;
            break;
        }
        if (str[i] == '.')
        {
            break;
        }    
    }
    
    if (i < 0)
    {
        _temp += ".";
    }
    else
    {
        _temp = _temp.substr(0, i + 1);
    }
    *path = _temp;   
}

/**
 * 
 * @param title 
 * @param prompt 
 * @return
 *
 *
 */
char Utility::promptForYesNo(char *title, const char *prompt)
{
    cout << "\n" << title << "\n";
    string promptMessage = string(prompt) + "(y/n) ? [y]";
    
    string answer;
    getStr(promptMessage, answer);
    
    if(!answer.empty())
    {
        if ((answer[0] != 'n') && (answer[0] != 'N'))
        {
            return ('y');
        }
    }
    return ('n');
}

/**
 * 
 * @param title 
 * @param prompt 
 * @return
 *
 */
char Utility::promptForYesNo(const char *title, const char *prompt)
{
    cout << "\n" << title << "\n";
    string promptMessage = string(prompt) + "(y/n) ? [y]";
    
    string answer;
    getStr(promptMessage, answer);
    
    if(!answer.empty())
    {
        if ((answer[0] != 'n') && (answer[0] != 'N'))
        {
            return ('y');
        }
    }
    return ('n');
}

/**
 * 
 * @param msg 
 */
void Utility::error( char *msg,...)
{
    va_list ap;

    va_start(ap, msg);
    fprintf(stderr, "\n\nERROR: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
}

/**
 * 
 * @param msg 
 */
void Utility::warning( char *msg,...)
{
    va_list ap;

    va_start(ap, msg);
    fprintf(stderr, "\n\nWARNING: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
}

/**
 * 
 * @param msg 
 */
void Utility::error( const char *msg,...)
{
    va_list ap;

    va_start(ap, msg);
    fprintf(stderr, "\n\nERROR: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
}

/**
 * 
 * @param msg 
 */
void Utility::warning( const char *msg,...)
{
    va_list ap;

    va_start(ap, msg);
    fprintf(stderr, "\n\nWARNING: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
}

/**
 * 
 * @param msg 
 */
void Utility::info( char *msg,...)
{
    va_list ap;
    va_start(ap, msg);

    if(! quiet)
    {
        fprintf(stdout, "\n");
        vfprintf(stdout, msg, ap);
        va_end(ap);
    }
}

/**
 * 
 * @param msg 
 */
void Utility::info(const char *msg,...)
{
    va_list ap;
    if(! quiet)
    {
        va_start(ap, msg);
        fprintf(stdout, "\n");
        vfprintf(stdout, msg, ap);
        va_end(ap);
    }
}


    
/**
 * 
 */
void Utility::myname( char *myname)
{
    strcpy(myname, "clustalw\0");
}


    
/**
 * Change:
 * Mark 25-1-2007. I made this change to get around the problem of having to keep track
 * of output indexes in the alignment stage. This function returns the next unique id.
 */
unsigned long Utility::getUniqueSequenceIdentifier()
{
    static unsigned long nextSequenceIdentifier = MinIdentifier;
    return nextSequenceIdentifier++;
}

void Utility::setInfoLabelPtr(QLabel* ptrToLabelObj)
{
    // Dont do anything. This is here to allow the function to be added to QTUtility.
    // Polymorphism wont work with different functions in the classes.
}


double Utility::average(std::vector<double>& v)
{
    double tmp = 0.0;
    std::vector<double>::iterator i;

    if (v.size() == 0)
        return 0.0;

    for (i=v.begin(); i != v.end(); i++)
        tmp += *i;
    return (tmp / v.size());
}

double Utility::stdDev(std::vector<double>& v)
{
    std::vector<double>::iterator i;
    double tmp = 0.0;
    double avg = average(v);

    if (v.size() == 0)
        return 0.0;
    
    for(i=v.begin(); i != v.end(); ++i)
        tmp += (*i - avg) * (*i - avg);
    return sqrt(tmp / v.size());
}

double Utility::median(std::vector<double> v)
{
    // From Moo & Koenig, "Accelerated C++:
    typedef vector<double>::size_type vec_sz;
    vec_sz size = v.size();
    
    if (v.size() == 0)
        return 0.0;
    
    std::sort(v.begin(), v.end());
    vec_sz mid = size/2;
    return size % 2 == 0 ? (v[mid] + v[mid-1]) / 2 : v[mid];
}

}

