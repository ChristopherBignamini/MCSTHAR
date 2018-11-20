#ifndef HADRONIZATIONEXCEPTION_H
#define HADRONIZATIONEXCEPTION_H

#include <string>
#include <exception>

using namespace std;

// TODO: do I have copied this code from somewhere???
/**
* @brief Hadronization exception handling class
*
* @author Christopher Bignamini
*/
class HadronizationException: public exception
{

    public:
    
        /**
        * @brief Constructor
        *
        */
        HadronizationException(void);

        /**
        * @brief Constructor with error message
        *
        * @param i_errorMessage Error message
        */
        HadronizationException(const string& i_errorMessage);
    
        /**
        * @brief Constructor with error message and exception throwing method name
        *
        * @param i_errorMessage Error message
        * @param i_methodName Name of the exception throwing method
        */
        HadronizationException(const string& i_errorMessage,
                               const string& i_methodName);
        
        /**
        * @brief Constructor with error message, exception throwing method name and return value
        *
        * @param i_errorMessage Error message
        * @param i_methodName Name of the method associated to the error (e.g. the method that raises the exception)
        * @param i_returnValue Return value related to the raised exception
        */
        HadronizationException(const string& i_errorMessage,
                               const string& i_methodName,
                               int i_returnValue);
    
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Destructor
        *
        */
        virtual ~HadronizationException(void) throw();
    
        /**
        * @brief Exception merging method
        *
        * This method is used to merge exceptions coming from different source 
        * (e.g., different threads). The resulting exception message is a list
        * of the single error messages (plus the corresponding return values). 
        * The return value is not automatically modified during merging procedure,
        * the original one is therefore conserved. The setReturnValue method below
        * must be used to set a new value.
        *
        * @param i_hadronizationException New exception to be merged with existing one
        */
        void merge(const HadronizationException& i_hadronizationException);
    
        /**
        * @brief Error message set method
        *
        * @param i_errorMessage Exception error message
        */
        inline void setErrorMessage(const string& i_errorMessage){
            m_errorMessage=i_errorMessage;
            m_isMessageEmpty = false; }

        /**
        * @brief Error message get method
        *
        * @return Exception error message
        */
        const string& getErrorMessage(void) const;

        /**
        * @brief Return value set method
        *
        * @param i_returnValue Exception return value
        */
        inline void setReturnValue(const int i_returnValue) { m_returnValue = i_returnValue; }
    
        /**
        * @brief Return value get method
        *
        * @return Exception return value (-1 is returned if no error value has been provided to exception)
        */
        inline int getReturnValue(void) const { return m_returnValue; }
    
    private:
        
        /**
        * @brief Error message with method name builder method
        *
        * @param i_errorMessage Error message
        * @param i_methodName Exception throwing method name
        * @return Error message with exception throwing method name
        */
        static string buildMessageWithMethodName(const string& i_errorMessage,const string& i_methodName);
    
        /**
        * @brief Error message with return value builder method
        *
        * @param i_errorMessage Error message
        * @param i_methodReturnValue Exception return value to be merged with error message
        * @return Error message with thrown exception return value
        */
        static string buildMessageWithReturnValue(const string& i_errorMessage,const int i_methodReturnValue);
    
        /**
        * Error message
        */
        string m_errorMessage;

        /**
        * Exception definition status flag: true in case of undefined exception (no error message provided). False otherwise.
        */
        bool m_isMessageEmpty;
    
        /**
        * Return value
        */
        int m_returnValue;
    
        /**
        * Exception merging status: true if this no previous merge has been previously required. False otherwise.
        */
        bool m_isFirstMerging;
    
        /**
        * Default error message used if no other message is provided
        */
        static const string m_defaultErrorMessage;

};

#endif
