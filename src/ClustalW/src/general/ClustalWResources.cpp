/**
 * Implements a singleton that maintains program resources.
 * The single instance is (re)instantiated on demand like:
 *     Resources *res = Resources::Instance();
 *
 * 24-05-07,Nigel Brown(EMBL): created.
 * 3-7-07, Mark Larkin, modified this class for clustalw 
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ClustalWResources.h"
#include <iostream>
#include "clustalw.h"
#include <fstream>
using namespace std;
 
namespace clustalw
{

//environment variables
static const char *CLUW_INSTALL_DIR = "CLUW_INSTALL_DIR";

//return the sole instance
ClustalWResources *ClustalWResources::Instance() {
    static ClustalWResources instance;
    return &instance;
}

ClustalWResources::ClustalWResources()
{
    //defaultPath
    defaultPath = ".";

    //executablePath
    executablePath = ".";
    
    //installPath
    installPath = ".";
    char *env;
    if ((env = getenv(CLUW_INSTALL_DIR)) != 0) 
    {
        installPath = string(env);
    }

    homePath = "";
}

void ClustalWResources::setPathToExecutable(string path) 
{
    executablePath = dirname(path);
}

string ClustalWResources::dirname(string path) 
{
    string tempString;
    int size = path.size();
    tempString = path;
    for (int i = size - 1; i > 0; i--) 
    {
        if (tempString[i] == DIRDELIM) // Mark, no standard function in c++
        { 
            tempString.erase(i);
            break;
        }
    }
    return tempString;
}

void ClustalWResources::dump() 
{
    printf("%s => %s [%s]\n%s => %s\n%s => %s\n",
           "installPath", installPath.c_str(), CLUW_INSTALL_DIR,
           "executablePath", executablePath.c_str(),
           "homePath", homePath.c_str()
           );
}

string ClustalWResources::findFile(const char *file, const ClustalWResourcePathType where) const 
{
    return findFile(string(file), where);
}

string ClustalWResources::findFile(const string file, const ClustalWResourcePathType where) const 
{
    const string *path;
    ifstream ifs;

    switch (where) 
    {
        case InstallPath:
            path = &installPath;
            break;
        case ExecutablePath:
            path = &executablePath;
            break;
        case HomePath:
            path = &homePath;
            break;
        default:
            path = &defaultPath;
            break;
    }
    char delim[1];
    delim[0] = DIRDELIM;
    delim[1] = 0;
    
    string fileName = *path + string(delim) + file;

    ifs.open(fileName.c_str(), ifstream::in);
    if (ifs.fail()) {
        return string();
    }
    
    if (ifs.is_open() && ifs.good()) 
    {
        ifs.close();
        return fileName;
    }
    return string(); //not found/readable
}

// Search for a (string) file in a succession of likely locations and
// return the full path as (string).
//
string ClustalWResources::searchPathsForFile(const string fileName) const 
{
    string file;
    while (1) {
        file = findFile(fileName, InstallPath);
        if (file != "") break;
        
        file = findFile(fileName, ExecutablePath);
        if (file != "") break;
        
        file = findFile(fileName, HomePath);
        if (file != "") break;
        
        file = findFile(fileName);
        if (file != "") break;
        
        file = fileName; // give up
        break;
    }
    return file;
}

}
