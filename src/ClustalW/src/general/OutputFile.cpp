/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "OutputFile.h"
#include "utils.h"
#include "userparams.h"

namespace clustalw
{

OutputFile::OutputFile()
{

}

OutputFile::~OutputFile()
{
    // If it is open, close it and say that a file has been created!!!!!
    if(file.get())
    {
        file->close();
        utilityObject->info("%s file created:   [%s]\n", typeOfFileMsg.c_str(),
                                name.c_str());
    }
}

bool OutputFile::openFile(std::string* fileName, const std::string msg, const std::string* path, 
                      const std::string ext, const std::string fileType)
{
    if (fileName->empty())
    {
        *fileName = getOutputFileName(msg, *path, ext);
            
        if(fileName->empty())
        {
            return false;
        }
    }

    file.reset(new std::ofstream(fileName->c_str(), std::ofstream::trunc));
                
    if(!file->is_open()) 
    {
        utilityObject->error("Cannot open output file [%s]\n", fileName->c_str()); 
        return false;
    }
    name = *fileName; 
    typeOfFileMsg = fileType;
    
    return true;
}

bool OutputFile::isOpen()
{
    return file->is_open();
}

std::ofstream* OutputFile::getPtrToFile()
{
    return file.get();
}
                      
std::string OutputFile::getOutputFileName(const std::string prompt, std::string path, 
                                          const std::string fileExtension)
{
    std::string temp;
    std::string _fileName; // Will return this name.
    std::string message;
    _fileName = path + fileExtension;

    if(_fileName.compare(userParameters->getSeqName()) == 0) 
    {
        cerr << "WARNING: Output file name is the same as input file.\n";
        if (userParameters->getMenuFlag()) 
        {
            message = "\n\nEnter new name to avoid overwriting  [" + _fileName + "]: ";
            utilityObject->getStr(message, temp);
            if(temp != "")
            {
                _fileName = temp;
            }
        }
    }
    else if (userParameters->getMenuFlag()) 
    {

        message = prompt + " [" + _fileName + "]";
        utilityObject->getStr(message, temp);
        if(temp != "")
        {
            _fileName = temp;
        }
    }   
    return _fileName;

}

}


