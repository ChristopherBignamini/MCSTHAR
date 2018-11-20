#include "../Include/logMessage.h"
#include <ctime>
#include <iostream>

using namespace std;

void logMessage(const string& i_logMessage, const LogType i_logType)
{
    time_t rawtime;
    time(&rawtime);
    string completeLogMessage("MCSTHAR ");
    switch(i_logType)
    {
        case INFO:
            completeLogMessage += "info";
            break;
        case WARNING:
            completeLogMessage += "warning";
            break;
        case ERROR:
            completeLogMessage += "error";
            break;
        default:
            completeLogMessage += "info";
            break;            
    }
    completeLogMessage += " - " + (string(ctime(&rawtime))).substr(0,24) + ": " + i_logMessage;
    cout<<completeLogMessage<<endl;
    
    return;
}