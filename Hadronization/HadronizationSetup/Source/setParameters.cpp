#include <iostream>
#include "../Include/HadronizationSetup.h"

using namespace std;

bool setParameters(HadronizationSetup &io_setupParameters)
{
	int userSelection;
	double inputValue;
  
	cout<<"		"<<"1: Phen. parameters"<<endl;	
	cout<<"		"<<"2: Cluster merging"<<endl;	
	cout<<"		"<<"3: Ranlux seed"<<endl;	
	cout<<"		"<<"Insert and return: ";
	cin>>userSelection;
	cout<<endl;
	if(userSelection==1){
		cout<<"		"<<"---> Parameter?"<<endl;
		cout<<"		"<<"1: Energy density"<<endl;	
		cout<<"		"<<"2: Gamma_s"<<endl;	
		cout<<"		"<<"Insert and return: ";
		cin>>userSelection;
		cout<<endl;
		if(userSelection==1){
			cout<<"		"<<"---> Energy density value (>= "
			<<io_setupParameters.minEnergyDensity<<" GeV/fm^3 and <= "<<io_setupParameters.maxEnergyDensity<<" GeV/fm^3): ";
			cin>>inputValue;
			cout<<endl;
			if(inputValue>=io_setupParameters.minEnergyDensity && inputValue<=io_setupParameters.maxEnergyDensity){
				io_setupParameters.energyDensity = inputValue;
				return 0;
			}
			else{
				cout<<"		"<<"WARNING, problem with energy density limits!"<<endl<<endl;
				return 1;
			}
		}
		else if(userSelection==2){
			cout<<"		"<<"---> Gamma_s value (>= "
			<<io_setupParameters.minGammaS<<" and <= "<<io_setupParameters.maxGammaS<<"): "; 
			cin>>inputValue;
			cout<<endl;
			if(inputValue>=io_setupParameters.minGammaS && inputValue<=io_setupParameters.maxGammaS){
				io_setupParameters.gammaS = inputValue;
				return 0;
			}
			else{
				cout<<"		"<<"WARNING, problem with Gamma_s limits!"<<endl<<endl;
				return 1; 
			}			
		}
		else{
			cout<<"		"<<"WARNING, check your selection!"<<endl<<endl;
			return 1;	
		}		
	}
	else if(userSelection==2){
		cout<<"		"<<"---> Reclustering on/of?"<<endl;
		cout<<"		"<<"0: Off"<<endl;	
		cout<<"		"<<"1: On"<<endl;	
		cout<<"		"<<"Insert and return: ";
		cin>>userSelection;
		cout<<endl;
		if(userSelection==0||userSelection==1){
			if(userSelection==1){
                io_setupParameters.clusterMergingFlag = 1;
                cout<<"		"<<"---> Cut value: ";
                cin>>inputValue;
                cout<<endl;
                io_setupParameters.clusterMergingMinMass = inputValue;
                if(inputValue<=0)
                    cout<<"		"<<"WARNING, negative or zero cut!"<<endl<<endl;
                return 0;
			}
			else{
				io_setupParameters.clusterMergingFlag = 0;
				return 0;
			}
		}
		else{
			cout<<"		"<<"WARNING, check your selection!"<<endl<<endl;
			return 1;
		}
	}
	else if(userSelection==3){
		cout<<"		"<<"RANLUX seed (integer): ";
		cin>>userSelection;
		cout<<endl;
		if(userSelection>=0){
			io_setupParameters.randomNumberGeneratorSeed = userSelection;
			return 0;
		}
		else{
			cout<<"		"<<"WARNING, negative RANLUX seed"<<endl<<endl;
			return 1;
		}
	}
	else{
		cout<<"		"<<"WARNING, check your selection!"<<endl<<endl;
		return 1;	
	}
}