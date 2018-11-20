#ifndef COMPUTEWAVEFUNCTIONSTRANGECOMPONENT_H
#define COMPUTEWAVEFUNCTIONSTRANGECOMPONENT_H

#include <string>
using namespace std;

/**
* Log message type
*/
enum LogType
{
    INFO,
    WARNING,
    ERROR
};

/**
* @brief Standard output message logger function
*
* @param i_logMessage Message to be logged to standard output
* @param i_isWarning True if message is a warning, false otherwise
*
* @author Christopher Bignamini
*/
void logMessage(const string& i_logMessage, LogType i_logType=INFO);

#endif