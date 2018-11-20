#include "../Include/HadronizationException.h"
#include <sstream>

const string HadronizationException::m_defaultErrorMessage("Error in MCSTHAR++ code");// TODO: update

HadronizationException::HadronizationException(void)
                                              :m_isMessageEmpty(true)
                                              ,m_returnValue(-1)
                                              ,m_isFirstMerging(true)
{
}

HadronizationException::HadronizationException(const string& i_errorMessage)
                                              :m_errorMessage(i_errorMessage)
                                              ,m_isMessageEmpty(false)
                                              ,m_returnValue(-1)
                                              ,m_isFirstMerging(true)
{
}

HadronizationException::HadronizationException(const string& i_errorMessage,
                                               const string& i_methodName)
                                              :m_errorMessage(buildMessageWithMethodName(i_errorMessage,i_methodName))
                                              ,m_isMessageEmpty(false)
                                              ,m_returnValue(-1)
                                              ,m_isFirstMerging(true)
{
}

HadronizationException::HadronizationException(const string& i_errorMessage,
                                               const string& i_methodName,
                                               const int i_returnValue)
                                              :m_errorMessage(buildMessageWithMethodName(i_errorMessage,i_methodName))
                                              ,m_isMessageEmpty(false)
                                              ,m_returnValue(i_returnValue)
                                              ,m_isFirstMerging(true)
{
}

HadronizationException::~HadronizationException(void) throw()
{
}

void HadronizationException::merge(const HadronizationException& i_hadronizationException)
{
    if(m_isFirstMerging)
    {
        // Merge error messages

        // Check if error message already exist
        if(m_isMessageEmpty)
        {
            m_errorMessage = "Merged exception list:";
            m_isMessageEmpty = false;
        }
        else
        {
            // Store original message and add it to the list
            string initialMessageCopy(m_errorMessage);
            initialMessageCopy = buildMessageWithReturnValue(initialMessageCopy,m_returnValue);
            
            // Update global error message
            m_errorMessage = "Merged exception list:\nException: " + initialMessageCopy;
            
        }
        m_isFirstMerging = false;
    }
    
    // Update global error message with provided exception date
    m_errorMessage += "\nException: " +
                        buildMessageWithReturnValue(i_hadronizationException.m_errorMessage,
                                                    i_hadronizationException.m_returnValue);

    return;
}

const string& HadronizationException::getErrorMessage(void) const
{
    if(m_isMessageEmpty)
    {
        return m_defaultErrorMessage;
    }

    return m_errorMessage;
}

string HadronizationException::buildMessageWithMethodName(const string& i_errorMessage,const string& i_methodName)
{
    return "(" + i_methodName + "): " + i_errorMessage;
}

string HadronizationException::buildMessageWithReturnValue(const string& i_errorMessage,const int i_methodReturnValue)
{
    stringstream returnValueStr;
    returnValueStr<<i_methodReturnValue;
    return i_errorMessage + "\nError code: " + returnValueStr.str();// TODO: use "error code" everywhere
}
