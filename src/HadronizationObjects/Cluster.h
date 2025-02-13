#ifndef CLUSTER_H
#define CLUSTER_H

#include "HadronizationObject.h"

/**
* @brief Cluster representation class
*
* This class is derived from the HadronizationObject class and
* is used to represent the clusters involved in the hadronization
* process, provided by the external event generator or created 
* during the hadronization process itself. By default, a new cluster
* is considered as hadronizable (m_isHadronizable flag set to true).
* See HadronizationObject class description for on-shellness condition.
*
* @author Christopher Bignamini
*/
class Cluster : public HadronizationObject
{
        
    public:
        
        /**
        * @brief Constructor
        *
        */
        Cluster(void);
    
        /**
        * @brief Constructor
        *
        * @param i_isHadronizable Cluster hadronizability flag
        * @param i_parentIndexes Hadronization object first and second parent indexes in event record
        * @param i_strangeCharge Cluster strange charge
        * @param i_charmCharge Cluster charm charge
        * @param i_bottomCharge Cluster bottom charge
        * @param i_electricCharge Cluster electric charge
        * @param i_baryonicCharge Cluster baryonic charge
        * @param i_energyDensity Cluster energy density
        * @param i_p Cluster 4-momentum
        * @param i_index Cluster index in event record
        */
        Cluster(bool i_isHadronizable,
                const pair<unsigned int,unsigned int>& i_parentIndexes,
                int i_strangeCharge,
                int i_charmCharge,
                int i_bottomCharge,
                double i_electricCharge,
                double i_baryonicCharge,
                double i_energyDensity,
                const TLorentzVector& i_p,
                unsigned int i_index = 0);
        
        /**
        * @brief Constructor
        *
        * @param i_isHadronizable Cluster hadronizability flag
        * @param i_parentIndexes Hadronization object first and second parent indexes in event record
        * @param i_strangeCharge Cluster strange charge
        * @param i_charmCharge Cluster charm charge
        * @param i_bottomCharge Cluster bottom charge
        * @param i_electricCharge Cluster electric charge
        * @param i_baryonicCharge Cluster baryonic charge
        * @param i_energyDensity Cluster energy density
        * @param i_mass Cluster mass
        * @param i_index Cluster index in event record
        */
        Cluster(bool i_isHadronizable,
                const pair<unsigned int,unsigned int>& i_parentIndexes,
                int i_strangeCharge,
                int i_charmCharge,
                int i_bottomCharge,
                double i_electricCharge,
                double i_baryonicCharge,
                double i_energyDensity,
                double i_mass,
                unsigned int i_index = 0);

    
        /**
        * @brief Destructor
        *
        */
        ~Cluster(void);
        
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Object parent indexes in event record setter
        *
        * @param i_parentIndexes Hadronization object first and second parent indexes in event record
        */
        inline void setParentIndexes(const pair<unsigned int,unsigned int>& i_parentIndexes) { m_parentIndexes = i_parentIndexes; }
    
        /**
        * @brief Cluster hadronizability flag setter
        *
        * @param i_isHadronizable Cluster hadronizability flag
        */
        inline void setHadronizabilityFlag(const bool i_isHadronizable) { m_isHadronizable = i_isHadronizable; }
        
        /**
        * @brief Cluster hadronization energy density setter
        *
        * @param i_energyDensity Cluster hadronization energy density
        */
        inline void setEnergyDensity(const double i_energyDensity) { m_energyDensity = i_energyDensity; }

        /**
        * @brief Cluster data setter
        *
        * @param i_isHadronizable Cluster hadronizability flag
        * @param i_parentIndexes Hadronization object first and second parent indexes in event record
        * @param i_strangeCharge Cluster strange charge
        * @param i_charmCharge Cluster charm charge
        * @param i_bottomCharge Cluster bottom charge
        * @param i_electricCharge Cluster electric charge
        * @param i_baryonicCharge Cluster baryonic charge
        * @param i_energyDensity Cluster energy density
        * @param i_p Cluster 4-momentum
        */
        void setClusterData(bool i_isHadronizable,
                            const pair<unsigned int,unsigned int>& i_parentIndexes, 
                            int i_strangeCharge,
                            int i_charmCharge,
                            int i_bottomCharge,
                            double i_electricCharge,
                            double i_baryonicCharge,
                            double i_energyDensity,
                            const TLorentzVector& i_p);

        /**
        * @brief Cluster charm charge setter
        *
        * @param i_charmCharge Cluster charm charge
        */
        inline void setCharmCharge(const int i_charmCharge)
                                  {
                                      m_charmCharge = i_charmCharge;
                                      if(m_charmCharge)
                                      {
                                          m_isHeavyFlavored = true;
                                      }
                                      else
                                      {
                                          m_isHeavyFlavored = false;
                                      }
                                  }
    
        /**
        * @brief Cluster bottom charge setter
        *
        * @param i_bottomCharge Cluster bottom charge
        */
        inline void setBottomCharge(const int i_bottomCharge)
                                   {
                                       m_bottomCharge = i_bottomCharge;
                                       if(m_bottomCharge)
                                       {
                                           m_isHeavyFlavored = true;
                                       }
                                       else
                                       {
                                           m_isHeavyFlavored = false;
                                       }
                                   }

    
        /**
        * @brief Cluster hadronizability flag getter
        *
        * @return True if cluster is hadronizable, false otherwise
        */
        inline bool isHadronizable(void) const { return m_isHadronizable; }

        /**
        * @brief Cluster hadronization energy density getter
        *
        * @return Cluster hadronization energy density
        */
        inline double getEnergyDensity(void) const { return m_energyDensity; }

        /**
        * @brief Cluster heavy flavorness condition flag getter
        *
        * @return True if cluster is heavy flavored, false otherwise
        */
        inline bool isHeavyFlavored(void) const { return m_isHeavyFlavored; }

        /**
        * @brief Cluster parent indexes in event record getter
        *
        * @return Cluster first and second parent indexes in event record
        */
        inline const pair<unsigned int,unsigned int>& getParentIndexes(void) const { return m_parentIndexes; }

    private:
    
        /**
        * Cluster heavy flavorness condition flag (true if c,cbar,b,bbar constituent quark exist)
        */
        bool m_isHeavyFlavored;
    
        /**
        * Cluster hadronizability flag (if true the cluster is considered to be hadronizable)
        */
        bool m_isHadronizable;
    
        // TODO: the storage of this data member is useless...(unless we switch to running density)
        /**
        * Cluster hadronization energy density
        */
        double m_energyDensity;
    
        /**
        * Cluster first and second parent indexes in event record
        */
        pair<unsigned int,unsigned int> m_parentIndexes;
    
};

/**
* @brief Cluster merging operator
*
* @return Cluster obtained from input cluster merging (new cluster hadronizability flag is set to true)
*/
Cluster operator+(const Cluster& i_firstCluster,const Cluster& i_secondCluster);


#endif