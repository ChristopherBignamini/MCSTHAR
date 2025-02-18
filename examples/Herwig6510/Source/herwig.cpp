
// TODO: debug
#include <iostream>
#include <fstream>
#include <iomanip>

// C++/Fortran interface wrapper
#include "../../../Interfaces/Herwig6510/Include/hwcforMCSTHAR.h"
#include "../../../Interfaces/Herwig6510/Include/loadClusters.h"
#include "../../../Interfaces/Herwig6510/Include/updateEventRecord.h"
#include "../../../Interfaces/Herwig6510/Include/Herwig6510Wrapper.h"

// Hadronization process setup
#include "../../../Hadronization/HadronizationSetup/Include/HadronizationSetupLoader.h"

// Hadronization process
#include "../../../Hadronization/HadronizationSteps/Include/HadronizationHandler.h"
#include "../../../Utilities/Include/HadronizationException.h"

// HepMC event writing for Rivet analysis
//#include "../../../Analysis/HepMCEventStorage/Include/buildHepMCEventRecord.h"

using namespace std;

int main(){
	
    // TODO: switch to rivet+hepmc api direct usage!
    // unsigned int hepMCEVentCounter(1);
    // HepMCEventRecord::initHepMCStorage();
    // HepMCEventRecord::HepMCStorageData hepMCStorageData("hepmc.fifo");
    
    // Hadronization process error status
    unsigned int hadronizationErrorStatus;
    
    // Herwig common block initialization
    hwigin();

    // Process selection
    bool isPPBar = false;
    
    if(isPPBar)
    {
        // Herwig process setup
        for(int i=0;i<8;i++){
            hwbmch.PART1[i] = ' ';
            hwbmch.PART2[i] = ' ';
        }
        hwbmch.PART1[0] = 'P';
        hwbmch.PART1[1] = 'B';
        hwbmch.PART1[2] = 'A';
        hwbmch.PART1[3] = 'R';
        hwbmch.PART2[0] = 'P';
        hwproc.PBEAM1 = 900.;
        hwproc.PBEAM2 = hwproc.PBEAM1;
        hwproc.MAXEV = 10;
        hwproc.IPROC = 1500;
    }
    else
    {
        // Herwig process setup
        hwbmch.PART1[0] = 'E';
        hwbmch.PART1[1] = '+';
        hwbmch.PART2[0] = 'E';
        hwbmch.PART2[1] = '-';
        for(int i=2;i<8;i++){
            hwbmch.PART1[i] = ' ';
            hwbmch.PART2[i] = ' ';
        }
        hwproc.PBEAM1 = 45.6;
        hwproc.PBEAM2 = hwproc.PBEAM1;
        hwproc.MAXEV = 10;
        hwproc.IPROC = 100;    
    }
    
    // Herwig ranlux seed set
    hwevnt.NRN[0] = 10101;
    hwevnt.NRN[1] = 10101;
//    hwevnt.NRN[0] = 1000985863;
//    hwevnt.NRN[1] = 680318663;
//    hwevnt.NRN[0] = 1521516856;
//    hwevnt.NRN[1] = 419042031;
    
    // MCSTHAR++ setup file
    cout<<"MCSTHAR++ setup file"<<endl;
    string mcstharSetupFile;
    cin>>mcstharSetupFile;
    
	// Herwig number of events 
	cout<<"Number of events"<<endl;
	cin>>hwproc.MAXEV;
    
    // Herwig print options
    hwevnt.MAXPR = 10;
    hwpram.PRVTX = 0;
	   
    // Compute Herwig parameter-dependent constants
	hwuinc();

    // Herwig elementary process initialization
	hweini();
    
    try
    {
    
        // Load MCSTHAR++ setup file
        HadronizationSetupLoader hadronizationSetupLoader(mcstharSetupFile);
        
        // Create hadronization handler
        HadronizationHandler hadronizationHandler(hadronizationSetupLoader.getHadronizationSetup());

        // Retrieve cluster energy density parameter used for cluster construction
        const double clusterEnergyDensity(hadronizationHandler.getEnergyDensity());
        
        // Herwig event generation
        for(int i=1;i<=hwproc.MAXEV;i++){
                                    
            // Initialize Herwig event
            hwuine();
        
            // Generate Herwig hard subprocess
            hwepro();
            
            // Generate Herwig parton cascades
            hwbgen();
            
            // Do heavy Herwig objects decays
            hwdhob();

            // Build clusters
            hwcforMCSTHAR();

            // Check error status
            if(hwevnt.IERROR==0)
            {
                // Upload clusters for MCSTHAR++
                const vector<Cluster> herwigClusters(loadClusters(clusterEnergyDensity));

                // Run cluster hadronization with MCSTHAR++
                // Load Herwig cluster
                hadronizationHandler.loadClusters(herwigClusters);
                
                // Run cluster hadronization
                hadronizationErrorStatus = hadronizationHandler.runHadronization();
                
                // TODO: find final flow control strategy
                // TODO: add more details to exception messages
                // TODO: add class for hadronization statistics
                // Check MCSTHAR++ hadronization error status
                if(hadronizationErrorStatus==0)
                {                
                    // Update Herwig event data
                    updateEventRecord(hadronizationHandler.getEventRecord());

                    // Do Herwig unstable particle decays
                    hwdhad();
                    
                    // Do Herwig heavy flavour decays
                    hwdhvy();
                    
                    // Add Herwig soft underlying event if needed
                    hwmevt();
                    
                }
            }
            // TODO: implement event rejection for MCSTHAR and Herwig error
            // TODO: implement HadronizationHistory/Stat/Log/??? class to store event rejection statistics
            
            // Event generation completed, wrap up Herwig event
            hwufne();
            
            // Check event status (Herwig error) and run analysis
            if((hwevnt.IERROR==0) && hadronizationErrorStatus==0)
            {
//                HepMCEventRecord::buildHepMCEventRecord(hepMCStorageData,hepMCEVentCounter);
//                ++hepMCEVentCounter;
            }
                
//            else
//            {
//                if(hadronizationErrorStatus)
//                {
//                    cout<<"Hadronization error, event "<<i<<" error "<<hadronizationErrorStatus<<endl;
//                }
//                else
//                {
//                    cout<<"Herwig error, event "<<i<<endl;                
//                }
//                i = i - 1;
//            }
            
           // TODO: debug
           ofstream eventFile;
           eventFile.open("eventFile.txt",std::fstream::app);
           eventFile<<setprecision(9);
           eventFile<<i<<" "<<hwevnt.EVWGT<<endl;
           for(int objectIndex=0;objectIndex<MCSTHAR::hepevt.NHEP;++objectIndex)
           {
               if((abs(MCSTHAR::hepevt.IDHEP[objectIndex])>=22) &&
                  (MCSTHAR::hepevt.IDHEP[objectIndex]!=94))
               {
                   eventFile<<MCSTHAR::hepevt.IDHEP[objectIndex]<<
                   " "<<(MCSTHAR::hepevt.PHEP)[objectIndex][0]<<
                   " "<<(MCSTHAR::hepevt.PHEP)[objectIndex][1]<<
                   " "<<(MCSTHAR::hepevt.PHEP)[objectIndex][2]<<
                   " "<<(MCSTHAR::hepevt.PHEP)[objectIndex][3]<<
                   " "<<(MCSTHAR::hepevt.PHEP)[objectIndex][4]<<endl;
               }
           }
           eventFile.close();
        }
        
        // Close Herwig event generation
        hwefin();
                
    }
    catch(HadronizationException& ex)
    {
        cout<<ex.getErrorMessage()<<endl;
        return ex.getReturnValue();
    }
        
	return 0;
}

