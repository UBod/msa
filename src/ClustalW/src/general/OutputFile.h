/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
#include <memory>
#include <fstream>
namespace clustalw
{

class OutputFile
{
    public:
        OutputFile();
        ~OutputFile();
        bool openFile(std::string* fileName, const std::string msg, const std::string* path, 
                      const std::string ext, const std::string fileType);
        bool isOpen();
        //void writeToFile(std::string* info);
        std::ofstream* getPtrToFile();
    private:
        std::string getOutputFileName(const std::string prompt, std::string path, 
                                      const std::string fileExtension);
        std::auto_ptr<std::ofstream> file;
        std::string typeOfFileMsg; // used for closing message!
        std::string name;
};

}
#endif

