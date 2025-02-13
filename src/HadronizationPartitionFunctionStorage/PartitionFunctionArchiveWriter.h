#ifndef PARTITIONFUNCTIONARCHIVEWRITER_H
#define PARTITIONFUNCTIONARCHIVEWRITER_H

// TODO: switch to forward declaration
#include "../Include/PartitionFunctionArchiveFile.h"
#include "../Include/PartitionFunctionDataSummaryFile.h"
#include "../Include/PartitionFunctionDataFile.h"
#include "partitionFunctionDataConstants.h"
#include "../../HadronizationPartitionFunction/Include/PartitionFunctionGenerationData.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <string>
#include <vector>

using namespace std;

// TODO: add description of archive structure: files, folders, etc...
/**
* @brief Partition function archive writing class
*
* This class provides the methods for partition function data archive writing.
* The archive is composed as follows: 
*
* @author Christopher Bignamini
*/
// TODO: add description of archive structure: files, folders, etc... and behavior (e.g. existing charge conf, diff grid...)
class PartitionFunctionArchiveWriter
{
    public:
        
//        /**
//        * @brief Constructor
//        *
//        * After the creation of a writer instance making use of this constructor the partition
//        * function data files can be written using the writeData method.
//        *
//        * @param i_partitionFunctionArchiveFolder Partition function data archive storage folder
//        * @param i_microcanonicalParameterGridStructure Microcanonical parameter grid structure
//        * @throw HadronizationException in case of error during existing partition function archive
//        *        file parsing or new archive file object instance creation
//        */    
        /**
         * // TODO: doxy
         */
        // TODO: change name to PartitionFunctionDataWriter if this version is conserved
        PartitionFunctionArchiveWriter(const string& i_partitionFunctionArchiveFolder,
                                       const MicrocanonicalParameterGridStructure& i_microcanonicalParameterGridStructure,
                                       const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                       const vector<double>& i_massValues);

        /**
        * @brief Destructor
        */
        ~PartitionFunctionArchiveWriter(void);

        /**
        * @brief Single mass and energy density value pair partition function data writing method
        *
        * @param i_partitionFunctionData Single mass/energy density pair partition function data
        * @param i_massIndex Mass value index, with respect to the provided mass set in class constructor
        * // TODO: doxy
        * @throw HadronizationException in case of error during data, summary or archive file
        *        writing
        */
        void updatePartitionFunctionData(const SingleMassEnergyDensityPairPartitionFunctionData& i_partitionFunctionData,
                                         unsigned int i_massIndex);
    
    private:
    
        /**
        * @brief Private copy constructor
        */
        PartitionFunctionArchiveWriter(const PartitionFunctionArchiveWriter&);

        /**
        * @brief Private copy constructor
        */
        PartitionFunctionArchiveWriter& operator=(const PartitionFunctionArchiveWriter&);

        /**
        * @brief Partition function data archive file update method
        *
        * @param i_chargeConfiguration New charge configuration
        * @param i_chargeFolderName New charge configuration data folder (relative) name
        * @throw  HadronizationException in case of error during archive file updating
        */
    
        /**
        * // TODO: doxy
        */
        void updatePartitionFunctionArchiveFile(const string& i_partitionFunctionArchiveFullFileName,
                                                const ChargeConfiguration& i_chargeConfiguration,
                                                const string& i_chargeFolderName);
    
        /**
        * // TODO: doxy
        */
        void updatePartitionFunctionDataSummaryFile(const double i_massValue,
                                                    const string& i_partitionFunctionDataFileName);

        /**
        * // TODO: doxy
        */
        void updatePartitionFunctionDataFile(const SingleMassEnergyDensityPairPartitionFunctionData& i_partitionFunctionData);
    
        /**
        * @brief New charge configuration storage folder name creation method
        *
        * @param i_chargeConfiguration New charge configuration
        * @return New charge configuration data storage folder name string
        */
        // TODO: move to utils
        string createNewChargeConfigurationPartitionFunctionFolderName(const ChargeConfiguration& i_chargeConfiguration) const;
    
        /**
        * @brief Charge value to string conversion method
        *
        * @param i_chargeValue Single charge value
        * @return Charge value string
        */
        // TODO: move to utils
        string createChargeValueString(const int i_chargeValue) const;
        
        /**
        * Microcanonical parameter grid structure
        */
        const MicrocanonicalParameterGridStructure m_microcanonicalParameterGridStructure;
    
        /**
        * Mass value list
        */
        const vector<double> m_massValues;
    
        /**
        * Partition function archive file handler
        */
        PartitionFunctionArchiveFile* m_partitionFunctionArchive;

        /**
        * TODO: doxy
        */
        PartitionFunctionDataSummaryFile* m_partitionFunctionSummaryFile;

        /**
        * TODO: doxy
        */
        PartitionFunctionDataFile* m_partitionFunctionDataFile;
    
        /**
        * TODO: doxy
        */
        string m_currentChargeFullFolderName;
        
        /**
        * TODO: doxy
        */
        string m_partitionFunctionSummaryFullFileName;

        /**
        * TODO: doxy
        */
        string m_partitionFunctionDataFullFileName;
    
        /**
        * Mass index of current partition data file
        */
        unsigned int m_referenceMassIndex;
    
        /**
        * Partition function data file update flag: true if writing is required for new mass value. False otherwise.
        */
        bool m_isNewMassData;
};

#endif