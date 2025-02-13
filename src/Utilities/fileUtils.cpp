#include "../Include/fileUtils.h"
#include "../Include/HadronizationException.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <cerrno>
#include <cstring>

// TODO: add win/mac/etc support

// TODO: funzione copiata da http://www.linuxquestions.org/questions/programming-9/c-list-files-in-directory-379323/
vector<string> getDirectoryContent(const string& i_directoryName)
{
    DIR* directoryPointer;
    struct dirent* direntPointer;
    
    // Try to open directory 
    if((directoryPointer = opendir(i_directoryName.c_str())) == NULL)
    {
        // TODO: add missing exit code
        throw HadronizationException("Error during " + i_directoryName + " opening",
                                     __FUNCTION__,
                                     711);
    }
    
    // Load file names and close directory
    vector<string> o_fileList;
    while((direntPointer = readdir(directoryPointer))!= NULL)
    {
        o_fileList.push_back(string(direntPointer->d_name));
    }    
    closedir(directoryPointer);
    
    return o_fileList;
}

void createDirectory(const string& i_directoryName)
{
    // Try to create the directory 
    const int exitCode(mkdir(i_directoryName.c_str(),0777));
    if(exitCode !=0)
    {
        // TODO: add missing exit code
        throw HadronizationException("Error during " + i_directoryName + " directory creation: " + strerror(errno),
                                     __FUNCTION__,
                                     712);
    }
    
    return;
}
