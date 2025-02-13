#ifndef PARTITIONFUNCTIONHANDLING_H
#define PARTITIONFUNCTIONHANDLING_H

#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <vector>
#include <string>
#include <map>


#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionArchiveFile.h"



using namespace std;

/**
* @brief Single charge configuration partition function 1D (mass) grid structure
*
* @author Christopher Bignamini
*/
struct PartitionFunctionMassGridData
{
    /**
    * Number of mass grid points
    */
    unsigned int massValueNumber;
    
    /**
    * Ordered mass values
    */
    vector<double> mass;

    /**
    * Partition function values
    */
    vector<double> partitionFunction;
    
    
    /**
     * Partition function values
     */
    vector<double> partitionFunctionErr;

};

/**
* @brief 1D grid linear interpolation data structure
*
* @author Christopher Bignamini
*/
struct SingleParameterInterpolationData
{
    /**
    * Constructor
    */ 
    SingleParameterInterpolationData(void)
                                    :gridPoints(0,0)
                                    ,weights(0.,0.)
    {
    }

    /**
    * Interpolation point indexes (wrt grid points)
    */
    pair<unsigned int, unsigned int> gridPoints;

    /**
    * Interpolation weights
    */
    pair<double, double> weights;
};


// TODO: check all limit situation (max mass value or min and similar for gammas and density) as well as exceptions (interp outside ranges)
// TODO: define final strategy for partition function data handling
/**
* @brief Partition function handling class
*
* Partition function data handling class, responsible for the runtime loading and 
* storage of the required microcanonical partition function data set as well as
* for hadronization event weight normalization factor calculation, from the partition
* function data set itself (See ClusterHadronization class for general hadronization
* procedure description). Given the microcanonical model free parameters, namely 
* the energy density and strangeness suppression parameter (gammaS in what follows), 
* and the hadronizing cluster mass value, the corresponding microcanonical partition
* function value is obtained by means of an interpolation procedure on the provided data 
* set for the cluster specific charge configuration. The single charge configuration data set
* structure is represented by a 3D grid whose parameters are cluster mass, energy density 
* and gammaS. 
*
* // TODO: add description of 1D data grid usage and how/when it is built
* // TODO: add description of data folder name standard
*
* @author Christopher Bignamini
*/
class PartitionFunctionHandling
{
    public:
    

        PartitionFunctionHandling(const string& i_partitionFunctionDataFolder,
                                  const double i_energyDensity,
                                  const double i_gammaS,
                                  const bool isNewMethod);

    
        double getPartitionFunctionValue(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                         const double i_mass,
                                         unsigned int& io_calculationErrorStatus,
                                         const bool isNewMethod);
    
    
    double getPartitionFunctionValueErr(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                   const double i_mass,
                                                                   unsigned int& io_calculationErrorStatus,
                                        const bool isNewMethod);
    
        /**
        * @brief Constructor
        *
        * After the construction of a PartitionFunctionHandling instance the partition 
        * function value for the provided cluster can be obtained using the getPartitionFunctionValue 
        * method. All values will be calculated using the free parameters provided to the present 
        * constructor.
        * 
        * @param i_partitionFunctionDataFolder Folder (with full path) containing the partition function data set
        * @param i_energyDensity Cluster energy density parameter
        * @param i_gammaS Strangeness suppression parameter
        * @throw HadronizationException if no partition function set is found for the provided parameters, 
        *        in case of error in partition function set loading using the provided folder name or if 
        *        the provided interpolation parameters are outside the available ranges
        */
        PartitionFunctionHandling(const string& i_partitionFunctionDataFolder,
                                  double i_energyDensity,
                                  double i_gammaS);
    
        /**
         * @brief Destructor
         *
         */
        ~PartitionFunctionHandling(void);
        
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Partition function value get method
        *
        * This method returns the partition value corresponding to the provided mass, charge
        * configuration and for the hadronization free parameters previously used in class constructor.
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @param io_calculationErrorStatus Error status equal to 0 in case of interpolation correctly
        *        performed, in case of missing partition function charge configuration and 2 in case of
        *        partition function interpolation mass value outside available range (I/O parameter)
        * @throw HadronizationException if no partition function set is found for the provided parameters
        */
        double getPartitionFunctionValue(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                         double i_mass,
                                         unsigned int& io_calculationErrorStatus);
    
    double getPartitionFunctionValueErr(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                   const double i_mass,
                                        unsigned int& io_calculationErrorStatus);

    
    private:

    
        void loadPartitionFunctionArchiveStructureData(void);

        void setPartitionFunctionDataSetFileList(const unique_ptr<PartitionFunctionArchiveFile>& i_partitionFunctionArchiveFile);

        void setMicrocanonicalParameterGridStructure(const unique_ptr<PartitionFunctionArchiveFile>& i_partitionFunctionArchiveFile);

        void computeMicrocanonicalParameterGridInterpolationData(void);

        void updatePartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                           const string& i_dataFolder,
                                           const double isNewMethod);
    
        void loadPartitionFunctionDataFile(const string& i_partitionFunctionDataFullFileName,
                                           const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                           vector<double>& io_partitionFunctionGrid,
                                           vector<double>& io_partitionFunctionErrGrid,
                                           const bool isNewMethod);
    
        /**
        * @brief Partition function full data set folder name loading method
        *
        * @throw HadronizationException if no partition function data set is found 
        */
        void setPartitionFunctionDataSetFileList(void);

        /**
        * @brief Single charge configuration partition function 1D data set grid building method
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_dataFolder Data set folder corresponding to the provided charge configuraton
        * @throw HadronizationException in case of error during partition function 1D grid creation 
        *        (e.g., missing data for the provided charge configuration)
        */
        void updatePartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                           const string& i_dataFolder);
    
        /**
        * @brief Single partition function file loading method
        *
        * @param i_dataFileName Data file name (full path)
        * @param i_chargeConfiguration Charge configuration
        * @param io_partitionFunctionGrid Full partition function grid for the provided charge configuration (I/O parameter)
        * @throw HadronizationException in case of error during file data loading (e.g., corrupted data file)
        */
        void loadPartitionFunctionDataFile(const string& i_dataFileName,
                                           const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                           vector<double>& io_partitionFunctionGrid,
                                           vector<double>& io_partitionFunctionGridErr);
    
        /**
        * @brief Single charge configuration partition function 1D grid building method
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_partitionFunctionGrid Full partition function grid for the provided charge configuration
        */
        void buildPartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                          const vector<double>& i_partitionFunctionGrid,
                                          const vector<double>& i_partitionFunctionErrGrid);

        /**
        * @brief Linear interpolation method
        *
        * This method computes the linear interpolation data corresponding for the
        * provided (ordered) data vector and interpolation point: the returned data 
        * contain both interpolation point indexes (wrt input vector elements) and weights.
        *
        * @param i_interpolationValue Interpolation point
        * @param i_valueSet Ordered data vector
        * @param io_interpolationErrorStatus Interpolation error status equal to false in case of interpolation
        *        correctly performed and true in case of interpolation point outside the available range
        * @return Interpolation data: interpolation point indexes (wrt input vector elements) and weights
        * @throw HadronizationException in case of interpolation point not inside the available data range
        */
        SingleParameterInterpolationData computeInterpolationData(const vector<double>& i_valueSet,
                                                                  double i_interpolationValue,
                                                                  bool& io_interpolationErrorStatus) const;

        /**
        * @brief Charge configuration retrieval method
        *
        * This method is used to retrieve from the partition function data folder name
        * the corresponding charge configuration. 
        *
        * @param i_dataFolderName Partition function data folder name
        * @return Charge configuration
        * @throw HadronizationException in case of not identified charge configuration
        */
        MCSTHAR::Utilities::ChargeConfiguration retrieveChargeConfiguration(const string& i_dataFolderName) const;
    
        /**
        * @brief Charge value retrieval method
        *
        * This method is used by the retrieveChargeConfiguration to extract from the
        * data folder name the single charge values.
        *
        * @param i_chargeDataString String containing the single charge value information
        * @return io_chargeValue Charge value (I/O parameter)
        * @throw HadronizationException in case of not identified charge value
        */
        void retrieveChargeData(const string& i_chargeDataString,int& io_chargeValue) const;

        /**
        * Partition function data set main folder
        */
        string m_partitionFunctionDataFolder;
    
        /**
        * Partition function 3D grid ordered energy density values
        */
        vector<double> m_energyDensityValues;
    
        /**
        * Partition function 3D grid ordered strangeness suppression parameter values
        */
        vector<double> m_gammaSValues;
    
        /**
        * Partition function full data set folder list
        */
        map<MCSTHAR::Utilities::ChargeConfiguration,string> m_partitionFunctionDataSetFolders;

        /**
        * Fixed energy density and strangeness suppression parameter partition function 
        * 1D grids
        */
        map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData> m_partitionFunctionMassGrids;

        /**
        * Energy density interpolation value
        */
        double m_energyDensity;
    
        /**
        * Strangeness suppression parameter interpolation value
        */
        double m_gammaS;
        
        /**
        * Energy density grid number of points
        */    
        unsigned int m_energyDensityPointNumber;

        /**
        * Strangeness suppression parameter grid number of points
        */
        unsigned int m_strangenessSuppressionParameterPointsNumber;
    
        /**
        * 2D grid (Energy density x gammaS) number of points
        */
        unsigned int m_grid2DDimension;
    
        /**
        * 2D grid (Energy density and gammaS) linear interpolation indexes
        */
        unsigned int m_index00;
        unsigned int m_index10;
        unsigned int m_index01;
        unsigned int m_index11;
    
        /**
        * 2D grid (Energy density and gammaS) linear interpolation weights
        */
        double m_weight00;
        double m_weight10;
        double m_weight01;
        double m_weight11;

        /**
        * Interpolation data availability status
        */
        bool m_isInterpolationDataAvailable;
    
        /**
        * 1D grid interpolation data availability flag
        */
        bool m_is1DGridDataAvailable;
};

#endif
