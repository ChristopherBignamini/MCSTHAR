#include <iostream>
#include "../Include/runSetupInterface.h"
#include "../Include/setParameters.h"

using namespace std;

void runSetupInterface(HadronizationSetup &io_setupParameters)
{
	
	bool setupCheck1,setupCheck2; 
	int userInput;

    
	cout<<endl<<endl;
	cout<<"		"<<"///////////////////////////////////////////////"<<endl;
	cout<<"		"<<"///////////////////////////////////////////////"<<endl;
	cout<<"		"<<"//////////                           //////////"<<endl;
	cout<<"		"<<"//				             //"<<endl;
	cout<<"		"<<"//                MCSTHAR++                  //"<<endl;
	cout<<"		"<<"//				             //"<<endl;
	cout<<"		"<<"//////////           for             //////////"<<endl;
	cout<<"		"<<"//                                           //"<<endl;
	cout<<"		"<<"//               HERWIG 6.510                //"<<endl;
	cout<<"		"<<"//                                           //"<<endl;
	cout<<"		"<<"//                                           //"<<endl;
	cout<<"		"<<"///////////////////////////////////////////////"<<endl;
	cout<<"		"<<"//                                           //"<<endl;	
	cout<<"		"<<"//  Authors: F.Becattini, C.Bignamini        //"<<endl;
	cout<<"                "<<"//           and F.Piccinini                 //"<<endl;
	cout<<"		"<<"//  Code monkey: C.Bignamini                 //"<<endl;
	cout<<"		"<<"//                                           //"<<endl;	
	cout<<"		"<<"///////////////////////////////////////////////"<<endl<<endl;
	
	do{
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"//              MCSTHAR++ SETUP              //"<<endl;
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
                
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"//        PHENOMENOLOGICAL PARAMETERS:       //"<<endl;
		cout<<"		"<<"//                                           //"<<endl;
		cout<<"		"<<"// 1) Energy density: "<<io_setupParameters.energyDensity<<" GeV/fm^3          //"<<endl;
		cout<<"		"<<"// 2) Gamma_s: "<<io_setupParameters.gammaS<<"                          //"<<endl;		
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
        
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		cout<<"		"<<"//                                           //"<<endl;
		if(io_setupParameters.clusterMergingFlag==0){
			cout<<"		"<<"//            CLUSTER MERGING: Off           //"<<endl;
			cout<<"		"<<"//                                           //"<<endl;
			cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		}
		else{
			cout<<"		"<<"//            CLUSTER MERGING: On            //"<<endl;
			cout<<"		"<<"//                                           //"<<endl;
			cout<<"		"<<"// Mass cut value: "<<	io_setupParameters.clusterMergingMinMass<<" GeV                   //"<<endl;
			cout<<"		"<<"//                                           //"<<endl;
			cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		}
        
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"//          RANLUX SEED: "<<io_setupParameters.randomNumberGeneratorSeed<<"             //"<<endl;
		cout<<"		"<<"//                                           //"<<endl;	
		cout<<"		"<<"///////////////////////////////////////////////"<<endl;

		do{
			cout<<"		"<<"---> Main:"<<endl;	 
			cout<<"		"<<"0: run"<<endl;
			cout<<"		"<<"1: parameter settings"<<endl;
			cout<<"		"<<"Insert and return: ";
			cin>>userInput;
			if(userInput==0){
				setupCheck1 = 0;
				setupCheck2 = 0;
			}
			else if(userInput==1){
				setupCheck1 = 1;
				setupCheck2 = setParameters(io_setupParameters);
			}
			else{
				cout<<"		"<<"WARNING, check your selection!"<<endl<<endl;
				setupCheck2 = 1;
			}
		}while(setupCheck2);
				
		if(io_setupParameters.clusterMergingFlag)
			if(io_setupParameters.maxClusterMass<=io_setupParameters.clusterMergingMinMass){
                io_setupParameters.clusterMergingFlag = 0;
				cout<<"		"<<"WARNING, no cluster mass range"<<endl;
				setupCheck1 = 1;
			}
		
	}while(setupCheck1);
	
}
