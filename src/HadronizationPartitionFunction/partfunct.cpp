#include <iostream>
#include <fstream>
#include <cmath>
#include "../Include/partfunct.h"
#include "../../../Utilities/Include/Constants.h"

partfunct::partfunct(bool flagBin,bool flagSin,bool flagQin,bool interin,unsigned int ordin,
					 unsigned int gammasnumin,unsigned int rhonumin,double gammasminin,
					 double gammasmaxin,double rhominin,double rhomaxin,
					 int Sin,double Bin,double Qin,int &Cin,int &Btin,double &Min)
{
	unsigned int i,j;
		
	//Inizializzazione dati funzione di partizione
	
	
	M = Min;
	flagB = flagBin;
	flagS = flagSin;
	flagQ = flagQin;
	inter = interin;
	ord = ordin;
	gammasnum = gammasnumin;
	rhonum = rhonumin; 
	S = Sin;
	B = Bin;
	Q = Qin;
	C = Cin;
	Bt = Btin;
	
	heavy = 0;
	
	gammasmin = gammasminin;
	gammasmax = gammasmaxin;
	rhomin = rhominin;
	rhomax = rhomaxin;

	T = 0.160;
	
	numps = 1000000;
	nummult = 10000;
	
	numgridp = gammasnum*rhonum;
	
	//Inizializzazione file di output e input
	
	outfile[0] ='Z';
	outfile[1] = '0';
	outfile[2] = '0';
	outfile[3] = '0';
	outfile[4] = '0';
	outfile[5] = '.';
	outfile[6] ='d';
	outfile[7] = 'a';
	outfile[8] = 't';
	outfile[9] = '\0';
		
	//Inizializzazione vettori funz di partizione
	
	for(j=0;j<maxnumgridp;j++){
		partv[j] = 0.e0;
		partv2[j] = 0.e0;
	}	
	
	//Costruzione vettore gammas, rho e mass
		
	if(gammasnum>1)
		deltagms = (gammasmax - gammasmin)/(1.*(gammasnum - 1));
	if(rhonum>1)
		deltarho = (rhomax - rhomin)/(1.*(rhonum - 1));
		
	for(i=0;i<gammasnum;i++)
		gammasv[i] = gammasmin + i*deltagms;
	
	for(i=gammasnum;i<gammasnummax;i++)
		gammasv[i] = 0.e0;

	for(i=0;i<rhonum;i++)
		rhov[i] = rhomin + i*deltarho;
	
	for(i=rhonum;i<rhonummax;i++)
		rhov[i] = 0.e0;
		
	for(j=0;j<rhonum;j++)
		Vv[j] = (M/rhov[j])*vconv;
	
	for(j=rhonum;j<rhonummax;j++)
		Vv[j] = 0.e0;
	
	//Per il momento solo Boltzmann e senza interazioni!!!
	
	ord = 0;
	inter = 0;
	
	/////////////////////////////
	
	//WARNING Inizializzazione variabili di sampling
	
//	rhosampl = (rhomax + rhomin)/2.e0;//Questa va tunata!
	rhosampl = rhomin;
//	rhosampl = rhomax;
	
	/////////////////////////////	
	
	return;
};

partfunct::~partfunct()
{
	
};

void partfunct::mccalc(unsigned int &inmind,vector<int> &strunflmes,vector<double> &mixmodcs2,hadrongr &gr,TRandom1 &r)
{
		
	
	int num,nummcr,N,tmpnum,Ssum,numpoints;
	
	//"Parameters"
	double twopi6 = pow(2.*pi,6);
	
	bool fst,stop,checkmulsts;
	double Vsampl,modcss,microw,microwerr,tmp,tmp2,tmpV; 

	
	unsigned int ik,il,im,ii;	
	unsigned int auxsize = 0;
	
	unsigned int l;
	int i,j,k;//Usato da mult_to_kin per l'informazione sullo sp.fasi a disposizione, identifica l'indice di massa di apertura del canale
	double weightv[numgridp],weightv2varphi[numgridp],varomegaphi[numgridp];//max reticolo in 100*100 step
	
	kingen_data kingen_in; 
	
	num = nummult;
	nummcr = numps;
		
	tmpnum = 0;
	heavy = 0;
	
	numpoints = 0;
	for(ii=0;ii<numgridp;ii++){
		weightv[ii] = 0.e0;	
		weightv2varphi[ii] = 0.e0;
		varomegaphi[ii] = 0.e0;
	}
	
	multsampl mult(flagB,flagS,flagQ,inter,B,S,Q,C,Bt);
	
	//Qua inizia il sampling dei charm e beauty
		
	if((C==0) && (Bt==0)){
		
		heavy = 0;
		masssampl = M;
		Vsampl = masssampl/rhosampl*vconv;//Questa va capita
		
		gr.sethadrmean(Vsampl,T);
				
		mult.setprob(gr);
		
	}
	else{
		bool checkchann;
		heavy = 1;
		checkchann = gr.setheavyprob(C,Bt,T);
		if(checkchann){
			cout<<"No channels available!"<<endl;
			print(inmind,numpoints,varomegaphi);
			return;
		}
		mult.setheavyprob(gr);
	}
	
	/////////////////////////////////////////
	
	int i2 = 0;
	int i3 = 0;
	int i4 = 0;
	
	i4 = 50;
	i2 = 100;
	fst = 1;
	
	stop = 0;
		
	for(i=1;i<=num;i++){
		
		if(stop)
			break;
		
		tmpnum++;
				
		//Qua inizia il sampling dei charm e beauty
		
		if(heavy){
	
			double tmpmass;
						
			do{
				
				mult.heavyhadrsampl(r);
			
				tmpmass = mult.getheavypart(gr).getmass();
				
				masssampl = M - tmpmass;
								
//				masssampl = M;
				
				Vsampl = masssampl/rhosampl*vconv;//Questa va capita
				
				gr.sethadrmean(Vsampl,T);
												
				mult.setresmult(gr);
												
				checkmulsts = mult.sampling(heavy,gr,r);
				
				if(checkmulsts){
					kingen_in.masscheck = 1;
					continue;
				}
				
				kingen_in = mult_to_kin(mult,gr);
				
			}while(kingen_in.masscheck);	
						
		}
		else{
			
			do{
				mult.sampling(heavy,gr,r);
				kingen_in = mult_to_kin(mult,gr);
			}while(kingen_in.masscheck);			
			
		}
		
		//////////////////////////////////////////
		
				
		N = kingen_in.fs;
			
					
		if(N==2){
		
			double m1,m2;
			double Tk,multw;
			
			multw = mult.getweight(); 
			
			m1 = (kingen_in.mfs).at(0);
			m2 = (kingen_in.mfs).at(1);
			
			////Ciclo sulle gammas, rho (controllare che nel calcolo dei pesi non siano 
			////state fatte semplificazioni illecite)
		
			numpoints = numpoints + 1;
				
			microw = 4.*pi/M;
			Tk = sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2)))/(2.*M);
			microw = microw*Tk*sqrt(m1*m1+Tk*Tk)*sqrt(m2*m2+Tk*Tk);	
			microw = microw*kingen_in.spinfact/(twopi6*kingen_in.multfact);
						
		//	for(k=0;k<7;k++)
//				for(j=0;j<(mult.getmult()).at(k);j++)
//					cout<<(mult.parttrasl(gr,(mult.getpart(k)).at(j),k)).getpid()<<" ";
//			cout<<endl;
			
			auxsize = 0;
			for(ik=0;ik<gammasnum;ik++){
				for(il=0;il<rhonum;il++){
					im = auxsize + il;
					tmpV = Vv[il];
					weightv[im] = microw*tmpV*tmpV/multw;
					weightv2varphi[im] = 0.e0;
				}
				auxsize = im + 1;
			}
			
		}
		else{
			
			double multw;
			double twopi3N = pow(2.*pi,3*N); 
			
			multw = mult.getweight();
			
			for(int f=0;f<N;f++)
				inputdata.m[f]=(kingen_in.mfs).at(f);
			for(int f=N;f<150;f++)
				inputdata.m[f]=0.;
			inputdata.npointsmc = nummcr;
			inputdata.n = N;
			
			
			numpoints = numpoints + 1;
			inputdata.mass = M;
			omega();
			
			if(outputdata.phi==-1.e0){
				cout<<"Run terminato al sampling "<<i<<endl;
				return;
			}
			
			////Ciclo sulle gammas e sulle rho
			
			microw = outputdata.phi*kingen_in.spinfact/(twopi3N*kingen_in.multfact);
			microwerr = outputdata.ephi*kingen_in.spinfact/(twopi3N*kingen_in.multfact);

//			for(k=0;k<7;k++)
//				for(j=0;j<(mult.getmult()).at(k);j++)
//					cout<<(mult.parttrasl(gr,(mult.getpart(k)).at(j),k)).getpid()<<" ";
//			cout<<endl;
							
			auxsize = 0;
			for(ik=0;ik<gammasnum;ik++){
				for(il=0;il<rhonum;il++){
					im = auxsize + il;
					tmpV = Vv[il];
					weightv[im] = microw*pow(tmpV,N)/multw;
					weightv2varphi[im] = microwerr*pow(tmpV,N)/multw;
				}
				auxsize = im + 1;
			}
		}
	
		
		//Soppressione stranezza
		Ssum = 0;
		
		k = 0;
		for(j=0;j<(mult.getmult()).at(k);j++)
			Ssum = Ssum - ((mult.parttrasl(gr,(mult.getpart(k)).at(j),k)).getstrange());
		k = 1;
		for(j=0;j<(mult.getmult()).at(k);j++)
			Ssum = Ssum + ((mult.parttrasl(gr,(mult.getpart(k)).at(j),k)).getstrange());	
		Ssum = Ssum + (mult.getmult()).at(2) + (mult.getmult()).at(3); 
		
		//Heavy flavour
		if(heavy){
			Ssum = Ssum + abs(mult.getheavypart(gr).getstrange());
		}
		///
		
		////Ciclo sulle gammas e rho
				
		auxsize = 0;
		for(ik=0;ik<gammasnum;ik++){
			tmp = pow(gammasv[ik],Ssum);
			for(il=0;il<rhonum;il++){
				im = auxsize + il;
				weightv[im] = weightv[im]*tmp;
				weightv2varphi[im] = weightv2varphi[im]*tmp;
			}
			auxsize = im + 1;
		}
						
		k = 6;
		for(j=0;j<(mult.getmult()).at(k);j++){
			for(l=0;l<strunflmes.size();l++)
				if((mult.parttrasl(gr,(mult.getpart(k)).at(j),k)).getpid()==strunflmes.at(l)){
					
					modcss = mixmodcs2.at(l);
					
					////Ciclo sulle gammas e rho
					
					auxsize = 0;
					for(ik=0;ik<gammasnum;ik++){
						tmp = gammasv[ik];
						for(il=0;il<rhonum;il++){
							im = auxsize + il;
							weightv[im] = weightv[im]*(modcss*(tmp*tmp - 1.e0) + 1.e0);
							weightv2varphi[im] = weightv2varphi[im]*(modcss*(tmp*tmp - 1.e0) + 1.e0);
						}
						auxsize = im + 1;
					}
					
				}
		}
		
		////Ciclo sulle gammas e rho
		
		for(l=0;l<numgridp;l++){
			tmp = weightv[l];
			partv[l] = partv[l] + tmp;
			partv2[l] = partv2[l] + tmp*tmp;
			tmp = weightv2varphi[l]*weightv2varphi[l];
			varomegaphi[l] = varomegaphi[l] + tmp;
		}
				
		
		//Ciclo di controllo errore relativo
		
		if(i>=100){
			if(i==i2){
				tmp = -1.;
				
				for(ii=0;ii<numgridp;ii++){
					tmp2 = ((partv2[ii] - partv[ii]*partv[ii]/(1.*i))/
							(partv[ii]*partv[ii]))*(1.*i)/(1.*i-1.) + varomegaphi[ii]/(partv[ii]*partv[ii]);
					if(abs(tmp2)<std::numeric_limits<double>::epsilon())
						tmp2 = 0.;
					if(tmp2>tmp)
						tmp = tmp2;
				}

				if(tmp<1.e-2){
					stop = 1;
				}
				else{
					if(i>=i3){
						i3 = i*10;
						i4 = (int) i/2;
					}
					i2 = i + i4;
				}
			}
			
		}
		
	}	
	
	print(inmind,numpoints,varomegaphi);
			
}

kingen_data partfunct::mult_to_kin(multsampl &mult,hadrongr &hadr){
	
/*Funzione simile a mult_to_kin di MCSTHAR, vede se c'è sp. fasi a 
 disposizione in funzione della massa totale e calcola 
 alcune quantità come molteplicità di spin totale e fattoriali.*/
		
	kingen_data tmp;
	int a,i,j,k;
	double tmpmass,tmpmult,checkmass;
	vector<int> tmppid;
	
	tmp.masscheck = 1;
	tmp.fs = 0;
	tmp.spinfact = 1.;
	tmp.multfact = 1.;
	
	checkmass = 0.;
	
	for(i=0;i<7;i++)
		for(j=0;j<(mult.getmult()).at(i);j++){
			tmpmass = (mult.getmass(i)).at(j);
			(tmp.mfs).push_back(tmpmass);
			checkmass = checkmass + tmpmass;
			tmp.fs++;
		}	
		
	//Aggiungo l'heavy 
	
	if(heavy){
		tmpmass = mult.getheavypart(hadr).getmass();
		(tmp.mfs).push_back(tmpmass);
		checkmass = checkmass + tmpmass;
		tmp.fs++;	
	}
	
	//Controllo se il decadimento è possibile e se ci sono almeno due part nello stato finale...
	//Il caso della particella singola va sistemato, comunque...
	
//	cout<<"Il controllo sul numero di part viene fatto anche altrove?? Notare che devono essere almeno due con anche l'heavy!"<<endl;
	
	if(tmp.fs>1){
		if(M>checkmass){
			
			tmp.masscheck = 0;
			
			//Calcolo fattore di spin
			
			for(j=0;j<7;j++)
				for(k=0;k<(mult.getmult()).at(j);k++){
					tmp.spinfact = tmp.spinfact*(2.*(mult.parttrasl(hadr,(mult.getpart(j)).at(k),j)).getspin() + 1.);
					tmppid.push_back((mult.parttrasl(hadr,(mult.getpart(j)).at(k),j)).getpid());
				}	
			
			if(heavy){
				tmp.spinfact = tmp.spinfact*(2.*(mult.getheavypart(hadr).getspin()) + 1.);
			}
			
			
			//Calcolo il fattore di molteplicità...
			
			for(;tmppid.size()>0;){
				tmpmult = 1.;
				a = tmppid.at(0);
				tmppid.erase(tmppid.begin());
				for(j=0;j<tmppid.size();j++){	
					if(tmppid.at(j)==a){
						tmpmult = tmpmult + 1.;
						tmp.multfact = tmp.multfact*tmpmult;
						tmppid.erase(tmppid.begin()+j);
						j = j - 1;
					}
				}
			}
			
		}   
	}
	return tmp;
}


void partfunct::partcalc(unsigned int &inmind,vector<int> &strunflmes,vector<double> &mixmodcs2,hadrongr &gr,TRandom1 &r){
	
	bool calct;
//	cin>>calct;
	
//	if(M>10.06996){
	if(C==0 && Bt==0){
		if(M>3.06996)
			calct = 1;
		else{
			calct = 0;
		}	
	}
	if(C!=0 && Bt==0){
		if(M>4.2)
			calct = 1;
		else{
			calct = 0;
		}	
	}
	if(C==0 && Bt!=0){
		if(M>7.2)
			calct = 1;	
		else{
			calct = 0;
		}	
	}
		
	if(calct){
		cout<<"mccalc"<<endl;
		mccalc(inmind,strunflmes,mixmodcs2,gr,r);
	}
	else{
		cout<<"sumcalc"<<endl;
		sumchcalc(inmind,strunflmes,mixmodcs2,gr);
	}

	cout<<"Output in "<<outfile<<endl;
	
	return;
}

void partfunct::sumchcalc(unsigned int &inmind,vector<int> &strunflmes,vector<double> &mixmodcs2,hadrongr &gr){
	
	
	//"Parameters"
	double twopi6 = pow(2.*pi,6);
	
	unsigned int j,ik,il,ii,l,numch;
	int fnd,nummcr,i,im,auxsize;
	vector<bool> check;
	vector<int> ind;
	vector< vector<unsigned int> > chs;
	vector<particle> vectpart;
	int number_of_channel,N,Ssum;
	double weightmult,weightexp,totwgtnorm,microw,microwerr,modcss,Vsampl,tmpV,tmp;
	double weightv[numgridp],weightverr[numgridp];//max reticolo in 100*100 step
	kingen_data kingen_in;
		
	nummcr = numps;
	
	for(j=0;j<numgridp;j++){
		weightv[j] = 0.e0;
		weightverr[j] = 0.e0;
	}
	
	numch = 0;
	totwgtnorm = 0.;
	weightexp = 0.;
	weightmult = 1.;
	
	masssampl = M;
	Vsampl = masssampl/rhosampl*vconv;//Questa va capita
	
	gr.sethadrmean(Vsampl,T);
		
	if((C==0) && (Bt==0)){
		
		heavy = 0;
		
	}
	else{
		gr.setheavyprob(C,Bt,T);
		heavy = 1;
	}
	
	/////////////////////////////////////////
	
	do{
		
		heavychannelbuilder(heavy,weightexp,gr,ind,fnd,chs,check,vectpart);
		
		number_of_channel = chs.size();
//		cout<<"Numero di canali totali "<<chs.size()<<endl;

		for(i=0;i<number_of_channel;i++){
			if(check.at(i)){
				
				channel_to_kingen_data(i,gr,ind,chs,vectpart,kingen_in,weightmult);				
				totwgtnorm = totwgtnorm + weightmult;
				
				N = kingen_in.fs;
				
				numch++;
//				cout<<endl;
//				cout<<"canale numero "<<numch<<endl;
//				cout<<"molteplicità "<<N<<endl;
				if(N==2){
					double m1,m2,Tk;
					m1 = (kingen_in.mfs).at(0);
					m2 = (kingen_in.mfs).at(1);
					
					////Ciclo sulle gammas e sulle rho (controllare che nel calcolo dei pesi non siano 
					////state fatte semplificazioni illecite}
										
					microw = 4.*pi/M;
					Tk = sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2)))/(2.*M);
					microw = microw*Tk*sqrt(m1*m1+Tk*Tk)*sqrt(m2*m2+Tk*Tk);	
					microw = microw*kingen_in.spinfact/(twopi6*kingen_in.multfact);
					
					auxsize = 0;
					for(ik=0;ik<gammasnum;ik++){
						for(il=0;il<rhonum;il++){
							im = auxsize + il;
							tmpV = Vv[il];
							weightv[im] = microw*tmpV*tmpV;
							weightverr[im] = 0.e0;
						}
						auxsize = im + 1;
					}
					
				}
				else{
					double twopi3N = pow(2.*pi,3*N); 
					
					for(int f=0;f<N;f++)
							inputdata.m[f]=(kingen_in.mfs).at(f);
					for(int f=N;f<150;f++)
						inputdata.m[f]=0.;
					inputdata.npointsmc = nummcr;
					inputdata.n = N;
					
					
					inputdata.mass = M;
					omega();
					
					if(outputdata.phi==-1.e0){
						cout<<"Run terminato sul canale "<<i<<" con N particelle"<<endl;
						return;
					}
					
					////Ciclo sulle gammas e sulle rho
					
					microw = outputdata.phi*kingen_in.spinfact/(twopi3N*kingen_in.multfact);
					microwerr = outputdata.ephi*kingen_in.spinfact/(twopi3N*kingen_in.multfact); 
					
					auxsize = 0;
					for(ik=0;ik<gammasnum;ik++){
						for(il=0;il<rhonum;il++){
							im = auxsize + il;
							tmpV = Vv[il];
							weightv[im] = microw*pow(tmpV,N);
							weightverr[im] = microwerr*pow(tmpV,N);
						}
						auxsize = im + 1;
					}
					
				}
				
				Ssum = 0;
				
				double tmpstrange;
				
				for(j=0;j<(chs.at(i)).size();j++){
					tmpstrange = vectpart.at(ind.at((chs.at(i)).at(j))).getstrange();
					if(tmpstrange>0.)
						Ssum = Ssum + tmpstrange;
					else if(tmpstrange<0.)
						Ssum = Ssum - tmpstrange;
				}
				
				////Ciclo sulle gammas e rho 
				
				auxsize = 0;
				for(ik=0;ik<gammasnum;ik++){
					tmp = pow(gammasv[ik],Ssum);
					for(il=0;il<rhonum;il++){
						im = auxsize + il;
						weightv[im] = weightv[im]*tmp;
						weightverr[im] = weightverr[im]*tmp; 
					}
					auxsize = im + 1;
				}
				
				int tmppid;
				
				for(j=0;j<(chs.at(i)).size();j++){
					tmppid = vectpart.at(ind.at((chs.at(i)).at(j))).getpid();
					for(l=0;l<strunflmes.size();l++)
						if(tmppid==strunflmes.at(l)){
							
							modcss = mixmodcs2.at(l);
							
							////Ciclo sulle gammas, rho e m
							
							auxsize = 0;
							for(ik=0;ik<gammasnum;ik++){
								tmp = gammasv[ik];
								for(il=0;il<rhonum;il++){
									im = auxsize + il;
									weightv[im] = weightv[im]*(modcss*(tmp*tmp - 1.e0) + 1.e0);
									weightverr[im] = weightverr[im]*(modcss*(tmp*tmp - 1.e0) + 1.e0);
								}
								auxsize = im + 1;
							}
						}				
				}
				
				////Ciclo sulle gammas e sulle rho
				
				for(l=0;l<numgridp;l++){
					tmp = weightverr[l];
					partv[l] = partv[l] + weightv[l];
					partv2[l] = partv2[l] + tmp*tmp;
				}
				
			}
		}
				
	}while(fnd);
		
	totwgtnorm = totwgtnorm*weightexp;
	
	print(inmind,numch,totwgtnorm);
	
	return;
	
}

void partfunct::heavychannelbuilder(bool &heavy,double &weightexp,hadrongr &gr,vector<int> &ind,int &fnd,vector< vector<unsigned int> > &chs,
									vector<bool> &check,vector<particle> &vectpart){
	
	
	static bool first = 1;
	static unsigned int N;
	
	int sumS;
	double sumQ,sumB,tmpmass;
	
	//Heavy flavour
	int sumC,sumBt;
	//
	
	vector<unsigned int> tmpv;
	static vector<double> mass;
	
	if(first){
		
		//Carico i canali di particella singola (i pesanti solo se la flag è attiva) solo al primo accesso
		//includendo solo le particelle con massa inferiore alla massa del cluster
		
		weightexp = 1.;
		
		if(heavy){
			if(C>0){
				for(unsigned int j=0;j<gr.getgr(7).size();j++){
					vectpart.push_back(gr.getgr(7).at(j));
				}
				for(unsigned int j=0;j<gr.getgr(9).size();j++){
					vectpart.push_back(gr.getgr(9).at(j));
				}
			}
			if(C<0){
				for(unsigned int j=0;j<gr.getgr(8).size();j++){
					vectpart.push_back(gr.getgr(8).at(j));
				}
				for(unsigned int j=0;j<gr.getgr(10).size();j++){
					vectpart.push_back(gr.getgr(10).at(j));
				}
			}
			if(Bt<0){
				for(unsigned int j=0;j<gr.getgr(11).size();j++){
					vectpart.push_back(gr.getgr(11).at(j));
				}
				for(unsigned int j=0;j<gr.getgr(13).size();j++){
					vectpart.push_back(gr.getgr(13).at(j));
				}
			}
			if(Bt>0){
				for(unsigned int j=0;j<gr.getgr(12).size();j++){
					vectpart.push_back(gr.getgr(12).at(j));
				}
				for(unsigned int j=0;j<gr.getgr(14).size();j++){
					vectpart.push_back(gr.getgr(14).at(j));
				}
			}
			
		}
		
		for(unsigned int j=0;j<7;j++){
			N = gr.getgr(j).size();
			for(unsigned int i=0;i<N;i++){
				vectpart.push_back(gr.getgr(j).at(i));
				if(!heavy)
					weightexp = weightexp*exp(-gr.getpartmean(j).at(i));
			}
		}
		
		for(unsigned int i=0;i<vectpart.size();i++)
			ind.push_back(i);
		
		//Riordino gli adroni 
		
		for(unsigned int n=1;n<=ind.size();n++){
			for(unsigned int i=0;i<ind.size()-1;i++){
				if(vectpart.at(ind.at(i)).getmass()<vectpart.at(ind.at(i+1)).getmass())
					swap(ind.at(i),ind.at(i+1));
			}
		}
		
		N = ind.size();
		
		cout<<"Warning, questo calcolo vale solo se gli adroni pesanti hanno massa maggiore di tutti i leggeri"<<endl;
		for(unsigned int i=0;i<N;i++){
			tmpv.clear();
			if(vectpart.at(ind.at(i)).getmass()<M){
				if(heavy){
					if(vectpart.at(ind.at(i)).getcharm()!=0 || vectpart.at(ind.at(i)).getbeauty()!=0){
						tmpv.push_back(i);
						chs.push_back(tmpv);
						check.push_back(0);
						mass.push_back(vectpart.at(ind.at(i)).getmass());
					}
				}
				else{
					tmpv.push_back(i);
					chs.push_back(tmpv);
					check.push_back(0);
					mass.push_back(vectpart.at(ind.at(i)).getmass());					
				}
			}
		}		
		
		first = 0;
	}
	
	
	fnd = 0;		
	unsigned int lim = chs.size();	
	
	for(unsigned int i=0;i<lim;i++){
				
		for(unsigned int j=chs.at(0).at(chs.at(0).size()-1);j<N;j++){
			
			tmpmass = mass.at(0) + vectpart.at(ind.at(j)).getmass();
			if(tmpmass<M){				

				fnd++;
				tmpv.clear();
				tmpv = chs.at(0);
				tmpv.push_back(j);
				mass.push_back(tmpmass);
				chs.push_back(tmpv);
				
				
				sumS = 0;
				sumQ = 0.;
				sumB = 0.;
				
				//Heavy flavour
				int numheavy = 0;
				sumC = 0;
				sumBt = 0;
				//
				for(unsigned int k=0;k<chs.at(0).size();k++){
					sumS = sumS + vectpart.at(ind.at(chs.at(0).at(k))).getstrange();
					sumQ = sumQ + vectpart.at(ind.at(chs.at(0).at(k))).getcharge();
					sumB = sumB + vectpart.at(ind.at(chs.at(0).at(k))).getbar();	
					
					//Heavy flavour
					if(abs(vectpart.at(ind.at(chs.at(0).at(k))).getcharm()) || abs(vectpart.at(ind.at(chs.at(0).at(k))).getbeauty()))
						numheavy = numheavy + 1;
					
					sumC = sumC + vectpart.at(ind.at(chs.at(0).at(k))).getcharm();
					sumBt = sumBt + vectpart.at(ind.at(chs.at(0).at(k))).getbeauty();
					//					
				}
				
				sumS = sumS + vectpart.at(ind.at(j)).getstrange();
				sumQ = sumQ + vectpart.at(ind.at(j)).getcharge();
				sumB = sumB + vectpart.at(ind.at(j)).getbar();
				
				//Heavy flavour
				if(abs(vectpart.at(ind.at(j)).getcharm()) || abs(vectpart.at(ind.at(j)).getbeauty()))
					numheavy = numheavy + 1;
								
				sumC = sumC + vectpart.at(ind.at(j)).getcharm();
				sumBt = sumBt + vectpart.at(ind.at(j)).getbeauty();
				//
				
				if(numheavy <= 1){
					if(sumS==S && sumB==B && sumQ==Q && sumC==C && sumBt==Bt)
						check.push_back(1);
					else
						check.push_back(0);
					
				}	
				else
					check.push_back(0);
								
			}
			
		}
		mass.erase(mass.begin());
		chs.erase(chs.begin());
		check.erase(check.begin());		
		
	}
	
	return;

}

	
void partfunct::channel_to_kingen_data(int &k,hadrongr &gr,vector<int> &ind,vector< vector<unsigned int> > &chs,vector<particle> &vectpart,
									   kingen_data &kingen_in,double &weightmult){
	
/*Questa funzione prepara l'input per il calcolo dello spazio fasi,
 calcola il fattore di spin e il fattore di molteplicità e il contributo
 di normalizzazione per il particolare canale (ricordare che i pesi dei canali
 calcolati con il metodo std non sono normalizzati).Inoltre prepara il vettore
 delle masse da passare a kingen. 
 Per il momento il tutto funziona solo senza interazioni, la massa che uso è quella
 memorizzata nei dati delle particelle!!!
 */
	
	
	//Dopo i check sul codice togliere checkmass e le altre cose inutili
	
	int a;
	double tmpmass,tmpmult;
	vector<int> tmppid;
	vector<double> tmpmean;
	
	weightmult = 1.;
	kingen_in.masscheck = 1;
	kingen_in.fs = 0;
	kingen_in.spinfact = 1.;
	kingen_in.multfact = 1.;
	kingen_in.mfs.clear();
	
	for(unsigned int j=0;j<(chs.at(k)).size();j++){
		tmpmass = vectpart.at(ind.at((chs.at(k)).at(j))).getmass();
		(kingen_in.mfs).push_back(tmpmass);
		kingen_in.fs++;
		//Controllo se c'è l'heavy e in tal caso riaggiorno le medie
		if(heavy)
			if(vectpart.at(ind.at((chs.at(k)).at(j))).getcharm() || vectpart.at(ind.at((chs.at(k)).at(j))).getbeauty()){
				masssampl = M - tmpmass;
				double Vsampl = masssampl/rhosampl*vconv;//Questa va capita
				gr.sethadrmean(Vsampl,T);
			}
	}	
		
	//Aggiorno il contributo esponenziale alle poisson del peso
	if(heavy)
		for(unsigned int j=0;j<7;j++){
			unsigned int N = gr.getgr(j).size();
			for(unsigned int i=0;i<N;i++){
				weightmult = weightmult*exp(-gr.getpartmean(j).at(i));
			}
		}

	
	if(kingen_in.fs>1){
		kingen_in.masscheck = 0;
		for(unsigned int j=0;j<(chs.at(k)).size();j++){
			kingen_in.spinfact = kingen_in.spinfact*(2.*vectpart.at(ind.at((chs.at(k)).at(j))).getspin() + 1.);
			tmppid.push_back(vectpart.at(ind.at((chs.at(k)).at(j))).getpid());
			//Hadrone leggero
			if((vectpart.at(ind.at((chs.at(k)).at(j))).getcharm()==0) && (vectpart.at(ind.at((chs.at(k)).at(j))).getbeauty()==0))
				weightmult = weightmult*gr.partmean(vectpart.at(ind.at((chs.at(k)).at(j))));
			else{
				//Hadrone pesante
				weightmult = weightmult*exp(-vectpart.at(ind.at((chs.at(k)).at(j))).getmass()/T)/gr.getheavynorm();
			}
		}	
				
		//Calcolo il fattore di molteplicità...
		
		for(;tmppid.size()>0;){
			tmpmult = 1.;
			a = tmppid.at(0);
			tmppid.erase(tmppid.begin());
			for(unsigned int j=0;j<tmppid.size();j++){	
				if(tmppid.at(j)==a){
					tmpmult = tmpmult + 1.;
					kingen_in.multfact = kingen_in.multfact*tmpmult;
					tmppid.erase(tmppid.begin()+j);
					j = j - 1;
				}
			}
		}
		weightmult = weightmult/kingen_in.multfact;
	}
	
	return;
	
}

void partfunct::print(unsigned int &mind,unsigned int &numch,double multnorm){

	unsigned int i,l;
	int auxsize;
	int m;
	
	
	auxsize = 0;
	m = 0;
	
	if(mind<10){	
		outfile[4] = '0' + mind;
	}
	else if(mind>=10 && mind<100){
		outfile[3] = '0' + (mind)/10;
		outfile[4] = '0' + (mind) - ((mind)/10)*10;
	}
	else if(mind>=100 && mind<1000){
		outfile[2] = '0' + (mind)/100;
		outfile[3] = '0' + (mind - ((mind)/100)*100)/10;
		outfile[4] = '0' + mind - ((mind)/10)*10; 
	}
	else{
		outfile[1] = '0' + (mind)/1000;
		outfile[2] = '0' + (mind-((mind)/1000)*1000)/100;
		outfile[3] = '0' + (mind - ((mind)/100)*100)/10;
		outfile[4] = '0' + mind - ((mind)/10)*10; 
	}
	
	ofstream out;	
	out.open(outfile);
	
	out<<numch<<endl;
	out<<rhosampl<<" "<<masssampl<<endl;	
	if(multnorm>0.)
		for(i=0;i<gammasnum;i++){
			for(l=0;l<rhonum;l++){
				m = auxsize + l;
				out<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<partv[m]/multnorm<<" "<<sqrt(partv2[m])/multnorm<<endl;
				cout<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<partv[m]/multnorm<<" "<<sqrt(partv2[m])/multnorm<<endl;
			}
			auxsize = m + 1;
		}
	else
		for(i=0;i<gammasnum;i++){
			for(l=0;l<rhonum;l++){
				m = auxsize + l;
				out<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<0.<<" "<<0.<<endl;
				cout<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<0.<<" "<<0.<<endl;
			}
			auxsize = m + 1;
		}
	
	out.close();
	out.clear();

	return;
	
};


void partfunct::print(unsigned int &mind,int numpoints,double varomegaphi[maxnumgridp]){
	
	
	unsigned int i,l;
	int auxsize;
	int m;
	double tmperr,tmpprec;
	
	auxsize = 0;
	m = 0;
	
	if(mind<10){	
		outfile[4] = '0' + mind;
	}
	else if(mind>=10 && mind<100){
		outfile[3] = '0' + (mind)/10;
		outfile[4] = '0' + (mind) - ((mind)/10)*10;
	}
	else if(mind>=100 && mind<1000){
		outfile[2] = '0' + (mind)/100;
		outfile[3] = '0' + (mind - ((mind)/100)*100)/10;
		outfile[4] = '0' + mind - ((mind)/10)*10; 
	}
	else{
		outfile[1] = '0' + (mind)/1000;
		outfile[2] = '0' + (mind-((mind)/1000)*1000)/100;
		outfile[3] = '0' + (mind - ((mind)/100)*100)/10;
		outfile[4] = '0' + mind - ((mind)/10)*10; 
	}
	
	ofstream out;	
	out.open(outfile);
	
	out<<numpoints<<endl;
	out<<rhosampl<<" "<<masssampl<<endl;	
	if(numpoints>0)
		for(i=0;i<gammasnum;i++){
			for(l=0;l<rhonum;l++){
				m = auxsize + l;
				tmperr = (partv2[m] - partv[m]*partv[m]/(1.*numpoints))/(1.*numpoints)/(1.*numpoints-1.)+varomegaphi[m]/(1.*numpoints*numpoints);
				tmpprec = tmperr/(partv[m]*partv[m])*(1.*numpoints*numpoints);
				if(abs(tmpprec)<std::numeric_limits<double>::epsilon())
					tmperr = 0.;
				out<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<partv[m]/numpoints<<" "<<sqrt(tmperr)<<endl;
				cout<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" "<<partv[m]/numpoints<<" "<<sqrt(tmperr)<<endl;
			}
			auxsize = m + 1;
		}
	else
		for(i=0;i<gammasnum;i++){
			for(l=0;l<rhonum;l++){
				m = auxsize + l;
				out<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" 0. 0."<<endl;
				cout<<M<<" "<<gammasv[i]<<" "<<rhov[l]<<" 0. 0."<<endl;
			}
			auxsize = m + 1;
		}
	out.close();
	out.clear();
	
	return;
};

//void partfunct::recintcalc(hadrongr &gr,vector<int> &strunflmes,vector<bool> &types,vector<double> &mixangl,
//						   vector<Z_interpolator> &Zmem){
//	
//	//Ciclo di integrazione sulla x
//	bool found;
//	bool twopch = 1;
//	bool checkdiff = 0;
//	unsigned int Nint,i,j,aux,count,checkint;
//	double x,eps,kernel_omega[10000],Int[10000],Intold[10000];	
//	
//	eps = 1e-2;
//	Nint = 4;
//	count = 0;
//	
//	x = 1.e0;
//	
//	bool wait;
//	
//	for(j=0;j<numgridp;j++)
//		Int[j] = 0.e0; 
//	
//	sumchcalc(twopch,strunflmes,types,mixangl,gr);
//
//	//Qua si può mettere un check. Oppure all'interno di integrand, in modo che quando
//	//kernel_omega è zero evita di fare l'intero ciclo di integrazione in quanto sarà sempre nullo!
//	
//	found = 0;
//	
//	wait = 1;
//	
//	for(;checkdiff==0;){
//	
//		checkdiff = 1;
//		
////		cout<<"Nel ciclo di integrazione"<<endl;
////		cout<<"Numero di intervalli "<<Nint<<endl;
//			
//		//L'integrazione è "backward" perchè se la omega è zero in x allora
//		//lo è anche in x'< x. Il controllo è nella variabile checkint, quando
//		//è 1 bisogna proseguire nell'integrazione all'indietro!
//		
//		for(i=Nint-1;i>0;i--){
//							
//			x = (1.e0*i)/(1.e0*Nint);
//						
//			checkint = integrand(x,gr,strunflmes,types,mixangl,kernel_omega,Zmem);
//									
//			if(found == 0){
//				if(checkint == 0){
////					cout<<"Esco dalla processo di integrazione, omega sempre nulla"<<endl;
//					break;
//				}
//				else if(checkint == 1){
////					cout<<"All'estremo superiore c'è sp.fasi disponibile ma bisogna aumentare gli int"<<endl;
//					break;
//				}
//			}
//					
//			found = 1;
//			
//			if(i==Nint-1)
//				for(j=0;j<numgridp;j++){
//					Int[j] = kernel_omega[j]*x/(1.e0*Nint); 
//				}			
//			else
//				for(j=0;j<numgridp;j++){
//					Int[j] = Int[j] + kernel_omega[j]*x/(1.e0*Nint);
//				}
//								
//			
//		}
//		
//		if(found == 0){
//			if(checkint==1){
////				cout<<"C'è spazio per integrare, aumento il numero di intervalli"<<endl;
//				Nint = Nint*2;
//				checkdiff = 0;
//			}
//			else{
////				cout<<"Non c'è nulla da integrare, la omega fa zero anche a max massa invariante, esco dalla funzione di integrazione"<<endl;
//				break;
//			}
//		}
//		else{
////			cout<<"Sto integrando nel supporto, inizio a contare i cicli di integrazione"<<endl;
//			if(count<2)
//				checkdiff = 0;
//			else if(count==2){
//				checkdiff = 0;
//				for(j=0;j<numgridp;j++)
//					Intold[j] = Int[j];
//			}
//			else if(count>2){
////				cout<<"Controllo convergenza integrazione"<<endl;
//				for(j=0;j<numgridp;j++){
////					cout<<j<<" "<<Int[j]<<" "<<Intold[j]<<" "<<abs(Int[j]-Intold[j])<<" "<<eps*Intold[j]<<endl;
//					if(abs(Int[j]-Intold[j])<eps*abs(Intold[j])){
//						checkdiff = checkdiff*1;
//					}
//					else{
//						if(Int[j] == 0.e0 && Intold[j] == 0.e0)
//							checkdiff = checkdiff*1;
//						else 
//							checkdiff = 0;
//					}
////					cout<<"valore checkdiff "<<checkdiff<<endl;
//					Intold[j] = Int[j];
//				}
//			}
//			
//			//		cout<<"checkdiff "<<checkdiff<<endl;
//			//		cout<<"count "<<count<<endl;
//			Nint = Nint*2;
//			count++;
//		}
//		
//		
//	}
//
//	
//	//cout<<"risultato integrazione"<<endl;
//	for(i=0;i<gammasnum;i++){
//		aux = i*rhonum;
//		for(j=0;j<rhonum;j++){
//			partv[aux+j] = partv[aux+j] + Int[aux+j]*Vv[j]/(16.e0*pi2*pow(M,5));
//			//			cout<<gammasv[i]<<" "<<rhov[j]<<" "<<partv[aux+j]<<endl;
//		}
//	}
//	
//	return;
//}
//double partfunct::kernel(double &m,double &spin,double &x){
//	
//	double k;
//	
//	k = 2.*spin + 1.;
//	k = k*(M - m)*(M - m);
//	k = k*(M*M + m*m - (M - m)*(M - m)*x*x)*(M*M + m*m - (M - m)*(M - m)*x*x);
//	k = k*sqrt((M*M - (x*(M - m) + m)*(x*(M - m) + m))*(M*M - (x*(M - m) - m)*(x*(M - m) - m)));
//	
//	return k;
//	
//}
//
//unsigned int partfunct::integrand(double &x,hadrongr &gr,vector<int> &strunflmes,vector<bool> &types,
//						 vector<double> &mixangl,double *kernel_omega,vector<Z_interpolator> &Zmem){
//
//	bool found;
//	unsigned int i,j,checkint;	
//	int Sr,k;
//	double mass,spin,Qr,Br,maxmassinv,omega[10000];
//	
////	cout<<"Massa cluster "<<M<<endl;
//	
////	cout<<"Calcolo l'integrando per x "<<x<<endl;
////	cout<<"La conf di carica è S"<<S<<" Q "<<Q<<" B "<<B<<endl;
//	
//	bool wait;
//	
//	for(i=0;i<numgridp;i++)
//		kernel_omega[i] = 0.;
//	
//	//Controllo che per almeno uno dei set di particelle ci sia un contributo 
//	//non nullo. In caso contrario quando esco avviso la funzione chiamante che
//	//è inutile che cerchi di integrare su valori ancora più piccoli di massa invariante
//	checkint = 0;
//	
//	for(i=0;i<7;i++){
//	
////		cout<<"Gruppo particelle "<<i<<endl;
//		
//		//Massima massa invariante in fase di integrazione
////		cout<<"lightest "<<gr.getlightind(i)<<" "<<(gr.getgr(i)).size()<<endl;
//		maxmassinv = M - (gr.getgr(i)).at(gr.getlightind(i)).getmass();
//		
////		cout<<"maxmassinv "<<maxmassinv<<endl;
//		
//		if(maxmassinv<=twompi0){
//			//Non c'è nulla da integrare
////			cout<<"Dal controllo sulla massa invariante massima non c'è nulla da integrare in questo set di particelle "<<maxmassinv<<" "<<twompi0<<endl;
//			continue;
//		}
//		
//		//C'è spazio di integrazione per almeno un set di particelle, controllo per il 
//		//particolare valore di x passato alla funzione.
//		checkint = 1;
//		
//		if(x*maxmassinv<=twompi0){
//			//Non c'è nulla da integrare
////			cout<<"Dal controllo sulla massa invariante non c'è nulla da integrare in questo set di particelle per questa x "<<x*maxmassinv<<" "<<twompi0<<endl;
////			cout<<"x "<<x<<" maxmassinv "<<maxmassinv<<endl;
//			continue;
//		}
//		
//		//Almeno per una particella c'è sp.fasi disponibile con questa x
//		checkint = 2;
//		
//		//Inizio a ciclare sulle particelle
////		cout<<"numero part nel gruppo "<<gr.getgr(i).size()<<endl;
//		for(j=0;j<gr.getgr(i).size();j++){
//			
////			cout<<"spin "<<(gr.getgr(i)).at(j).getspin()<<" massa "<<(gr.getgr(i)).at(j).getmass()<<endl;
//			
//			Zmemsize = Zmem.size();
//			
//			mass = (gr.getgr(i)).at(j).getmass(); 
//			maxmassinv = M - mass;	
//			
////			cout<<"Massa invariante locale "<<maxmassinv<<endl;
//			
//			if(x*maxmassinv<=twompi0){
//				//Nulla da integrare su questa particella 
////				cout<<"Dal controllo sulla massa invariante non c'è nulla da integrare per questa particella "<<x*maxmassinv<<" "<<twompi0<<endl;
//				continue;
//			}
//			
//			Sr = S - (gr.getgr(i)).at(j).getstrange();
//			Qr = Q - (gr.getgr(i)).at(j).getcharge();
//			Br = B - (gr.getgr(i)).at(j).getbar();
//			
//			cout<<"Confiurazione di carica Sr "<<Sr<<" Qr "<<Qr<<" Br "<<Br<<endl;
////			cout<<"Dimensione memoria "<<Zmemsize<<endl;
////			cin>>wait;
//			
//			found = 0;
//			
//			//Cerco in memoria se ho già la funzione di partizione che mi serve 
//			for(k=0;k<Zmemsize;k++){
//				if(Sr==Zmem.at(k).getS() && Qr==Zmem.at(k).getQ() && Br==Zmem.at(k).getB()){
//					cout<<"configurazione trovata k "<<k<<" SQB "<<Zmem.at(k).getS()<<" "<<Zmem.at(k).getQ()<<" "<<Zmem.at(k).getB()<<endl;
//					
//					//Controllo fino a dove posso integrare con le info da file			
//					
//					//cout<<"configurazione trovata, massa massima a disposizione "<<Zmem.at(k).getmaxmass()<<" massa max da int "<<maxmassinv<<endl;
//					//cin>>wait;
//					
//					if(Zmem.at(k).getmaxmass()>=maxmassinv){
//						//Ok, ho dati a sufficienza (non guardo la rho per ora)
//						//						cout<<"ok, interpolazione"<<endl;
//						//						cin>>wait;
//						maxmassinv = x*maxmassinv;
//						Zmem.at(k).getZ(maxmassinv,numgridp,omega);
//						if(Zmem.at(k).getcheck())
//							cout<<"Warning, something went wrong during the interpolation!"<<endl;
//					}
//					else{
//						//Una parte dell'informazione esiste, calcolo la rimanente
//						//e aggiorno Zmem ciclando su una partfunct temporanea
//						//						cout<<"no, devo calcolare quello che manca"<<endl;
//						//						cin>>wait;
//						respartcalc(Sr,Qr,Br,strunflmes,types,mixangl,gr,maxmassinv,k,Zmem);
//						maxmassinv = x*maxmassinv;
//						Zmem.at(k).getZ(maxmassinv,numgridp,omega);
//						if(Zmem.at(k).getcheck())
//							cout<<"Warning, something went wrong during the interpolation!"<<endl;
//						
//					}
//					
//					found = 1;
//					break;
//				}
//				else if(Sr==-Zmem.at(k).getS() && Qr==-Zmem.at(k).getQ() && Br==-Zmem.at(k).getB()){
//					
//					cout<<"configurazione trovata k "<<k<<" SQB "<<Zmem.at(k).getS()<<" "<<Zmem.at(k).getQ()<<" "<<Zmem.at(k).getB()<<endl;
//					
//					if(Sr==-Zmem.at(k).getS()){
//						cout<<"Cambio la conf originaria"<<endl;
//						Sr = -Sr;
//						Qr = -Qr;
//						Br = -Br;
//					}
//					
//					//Controllo fino a dove posso integrare con le info da file			
//					
//					//cout<<"configurazione trovata, massa massima a disposizione "<<Zmem.at(k).getmaxmass()<<" massa max da int "<<maxmassinv<<endl;
//					//cin>>wait;
//					
//					if(Zmem.at(k).getmaxmass()>=maxmassinv){
//						//Ok, ho dati a sufficienza (non guardo la rho per ora)
////						cout<<"ok, interpolazione"<<endl;
////						cin>>wait;
//						maxmassinv = x*maxmassinv;
//						Zmem.at(k).getZ(maxmassinv,numgridp,omega);
//						if(Zmem.at(k).getcheck())
//							cout<<"Warning, something went wrong during the interpolation!"<<endl;
//					}
//					else{
//						//Una parte dell'informazione esiste, calcolo la rimanente
//						//e aggiorno Zmem ciclando su una partfunct temporanea
////						cout<<"no, devo calcolare quello che manca"<<endl;
////						cin>>wait;
//						respartcalc(Sr,Qr,Br,strunflmes,types,mixangl,gr,maxmassinv,k,Zmem);
//						maxmassinv = x*maxmassinv;
//						Zmem.at(k).getZ(maxmassinv,numgridp,omega);
//						if(Zmem.at(k).getcheck())
//							cout<<"Warning, something went wrong during the interpolation!"<<endl;
//						
//					}
//					
//					found = 1;
//					break;
//				}
//			}
//			
//			if(found == 0){
//				
//				//Non ho trovato in memoria l'informazione necessaria, la calcolo e carico in Zmem 
//				//ciclando su una partfunct temporanea
//				
//				cout<<"configurazione non trovata SQB "<<Sr<<" "<<Qr<<" "<<Br<<endl;
//				
//				k = -1;
////				cout<<"non esiste la confiurazione in memoria, calcolo tutto "<<endl;
////				cin>>wait;
//				respartcalc(Sr,Qr,Br,strunflmes,types,mixangl,gr,maxmassinv,k,Zmem);
//				maxmassinv = x*maxmassinv;
//				Zmem.at(Zmem.size()-1).getZ(maxmassinv,numgridp,omega);
////				cout<<"valori ottenuti in interpolazione"<<endl;
////				for(k=0;k<numgridp;k++)
////					cout<<omega[k]<<endl;
////				cin>>wait;
////				cout<<"dimensione memoria "<<Zmem.size()<<endl;
////				if(Zmem.at(Zmem.size()-1).getcheck())
////					cout<<"Warning, something went wrong during the interpolation!"<<endl;
//////				cin>>wait;
//			}
//			
//			spin = (gr.getgr(i)).at(j).getspin(); 
//
////			cout<<"Check costruzione kernel_omega"<<endl;
//			for(k=0;k<numgridp;k++){
////				cout<<kernel(mass,spin,x)<<" "<<omega[k]<<endl;
//				kernel_omega[k] = kernel_omega[k] + kernel(mass,spin,x)*omega[k];
//			}
//		}
//	}
//	
//	return checkint;
//
//}
//
//void partfunct::respartcalc(int &Sr,double &Qr,double &Br,vector<int> &strunflmes,vector<bool> &types,vector<double> &mixangl,
//							hadrongr &gr,double &maxmassinv,int &k,vector<Z_interpolator> &Zmem)
//{
//	bool buildir,last;
//	unsigned int origsize;
//	double mininvm,tmpmass;	
//	
//	
//	if(k==-1){
//		Z_interpolator tmpZ_int(Sr,Qr,Br); 
//		origsize = 0;
//		mininvm = 0;
//		tmpmass = twompi0;
//		buildir = 1;
//		Zmem.push_back(tmpZ_int);
//		k = Zmem.size() - 1;
//	}
//	else{
//		origsize = Zmem.at(k).getnummass();
//		mininvm = Zmem.at(k).getmaxmass();
//		if(mininvm<2.){
//			tmpmass = mininvm + 0.01; 
//			buildir = 0;
//		}
//		else{
//			tmpmass = mininvm + mpi0;
//			buildir = 0;
//		}		
//	}
//	
//	last = 0;
//	
//	if(buildir){
//		string createdir("mkdir Z/SQB");
//				
//		if(Sr<0)
//			createdir += '0' + 1;
//		else 	
//			createdir += '0';
//		
//		if(abs(Sr)<10){
//			createdir += '0';
//			createdir += '0' + abs(Sr);
//		}
//		else{
//			createdir += '0' + abs(Sr)/10;
//			createdir += '0' + abs(Sr) - (abs(Sr)/10)*10;
//		}
//		
//		
//		
//		if(Qr<0)
//			createdir += '0' + 1;
//		else
//			createdir += '0';	
//
//		if(abs(Qr)<10){
//			createdir += '0';
//			createdir += '0' + (int) abs(Qr);
//		}
//		else{
//			createdir += '0' + ((int) abs(Qr))/10;
//			createdir += '0' + (int) abs(Qr) - (((int) abs(Qr))/10)*10;
//		}
//		
//				
//		if(Br<0)
//			createdir += '0' + 1;
//		else 
//			createdir += '0';	
//
//		if(abs(Br)<10){
//			createdir += '0';
//			createdir += '0' + (int) abs(Br);
//		}
//		else{
//			createdir += '0' + ((int) abs(Br))/10;
//			createdir += '0' + (int) abs(Br) - (((int) abs(Br))/10)*10;
//		}
//				
//		system(createdir.c_str());
//	}
//	
//	
//	for(unsigned int i=origsize;;i++){ 
//	
//		//Calcolo
//		
////		cout<<"Massa invariante nel ciclo "<<tmpmass<<endl;
//		
//		partfunct tmppart(flagB,flagS,flagQ,inter,ord,gammasnum,rhonum,gammasmin,gammasmax,rhomin,rhomax,Sr,Br,Qr,tmpmass);
//		tmppart.partcalc(i,strunflmes,types,mixangl,gr,k,Zmem);
//						
//		if(tmpmass<2.){
//			tmpmass = tmpmass + 0.01; 
//		}
//		else{
//			tmpmass = tmpmass + mpi0;
//		}
//		
//		if(last)
//			break;
//		
//		if(tmpmass>=maxmassinv)
//			last = 1;
//	}
//	
//	return;
//}

