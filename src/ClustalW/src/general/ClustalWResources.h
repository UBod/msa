/**
 * Implements a singleton that maintains program resources.
 * The single instance is (re)instantiated on demand like:
 *     Resources *res = Resources::Instance();
 *
 * 24-05-07,Nigel Brown(EMBL): created.
 * 3-7-07, Mark Larkin, modified this class for clustalw
 */

#ifndef RESOURCESCLUSTALW_H
#define RESOURCESCLUSTALW_H

#include <string>

using namespace std;

namespace clustalw
{


enum ClustalWResourcePathType {
    DefaultPath,
    InstallPath,
    ExecutablePath,
    HomePath
};

class ClustalWResources 
{

public:
    /* return the Resources singleton */
    static ClustalWResources *Instance();

    /* setters */
    void setPathToExecutable(std::string pathToFiles);
    
    /* getters */
    std::string getDefaultPath()    { return defaultPath; }
    std::string getInstallPath()    { return installPath; }
    std::string getExecutablePath() { return executablePath; }
    std::string getHomePath()       { return homePath; }

    std::string findFile(const char *file, const ClustalWResourcePathType where = DefaultPath) const;
    std::string findFile(const std::string file, const ClustalWResourcePathType where = DefaultPath) const;
    std::string searchPathsForFile(const std::string fileName) const;

    /* debug */
    void dump();

protected:
    /* hide the constructors */
    ClustalWResources();
    ClustalWResources(const ClustalWResources&);
    ClustalWResources& operator= (const ClustalWResources&);
    
    std::string dirname(std::string path);    

private:
    std::string defaultPath;
    std::string installPath;
    std::string executablePath;
    std::string homePath;
};

}
#endif //RESOURCES_H
