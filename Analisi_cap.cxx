#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <dirent.h>

using namespace std;

struct parametri {
	double a=0;
	double b=0;
	double S_a=0;
	double S_b=0;
	
	double Sxy=0;
};



parametri interpolazione_lin(vector<double>, vector<double>, double, double, int);


int main(){
	
	vector<double> intensita; vector<double> theta;
	double fondo = 2540.5; // +- 0.2 
	double	S_y=2.3, S_x=1; //pm 
	
	
	double t_i, I_i;
	while(cin >> t_i >> I_i){
		intensita.push_back(I_i);
		theta.push_back(t_i);
	}
	int M = theta.size();
	
	
	double Max_y=0, Max_xi=0, Max_xf, index;
	for(int i=0; i<M; i++){
		if(intensita[i] > 4*1000*1000){
			Max_xi = theta[i];
			break;
		}
		if (Max_y < intensita[i]){
			Max_y = intensita[i+1];
		}
	}
	
	for(int i=0; i<M; i++){
		if(intensita[i] > 4*1000*1000){
			if(intensita[i+1] <4*1000*1000){
				Max_xf = theta[i];
				break;
			}	
		}
	}
	
	
	cout << Max_xi << " " << Max_xf  << " " << Max_y << endl;
	
	double Max_x = (Max_xf + Max_xi)/2;
	double SMax_x = (1/sqrt(2)) * sqrt(2*pow(S_x/sqrt(12), 2));
	
	cout << Max_x << " +- " << SMax_x << endl;
	
	
//MINIMI SPERIMENTALI

	vector<double> m; 
	vector<double> ind;
	for(int i=0; i<M; i++){
		if (intensita[i-3]>intensita[i] && intensita[i+3]>intensita[i]){
			if (intensita[i-2]>intensita[i] && intensita[i+2]>intensita[i]){
					if (intensita[i] > 1.5*fondo){
						m.push_back(theta[i]);
						ind.push_back(i);
					}
			}
		}
	}

	
	
	int n=m.size();

	
	vector<double> minimi; vector<double> minn;vector<int> indi;
	for(int i=0; i<n; i++){
		if (abs(m[i])<500 && abs(m[i+1] - m[i])>20){
			minn.push_back(m[i]);
			indi.push_back(ind[i]);
		}
	}
	
	vector<int> indici;
	for(int i=0; i<minn.size(); i++){
		if (minn[i]!=296){
			if(minn[i]<0){
				if(abs(minn[i+1] - minn[i])>50){
					minimi.push_back(minn[i]);
					indici.push_back(indi[i]);
				}
			}
			else if(minn[i]>0){
				if(abs(minn[i-1] - minn[i])>50){
					minimi.push_back(minn[i]);
					indici.push_back(indi[i]);
				}
			}
			
		}
	}
	
	int K=minimi.size();
	for(auto i:minimi){
		cout << i << endl;
	}
	cout << endl;
	for(auto i:indici){
		cout << intensita[i] << endl;
	}
	
		
//CORREZIONE SCOSTAMENTO
	
	
	
	cout << endl << "Minimi corretti: " << endl;

	vector<double> corr_msper; vector<double> Scorr_msper;
	
	for(int i=0; i<K; i++){
		corr_msper.push_back(minimi[i] - Max_x);
		Scorr_msper.push_back(sqrt(pow(S_x, 2) + pow(SMax_x, 2)));
	}

	for(int i=0; i<K; i++){
		cout << corr_msper[i] <</* " " << Scorr_msper[i] <<*/ endl;
	}
	
	
//CONVERSIONE RADIANTI

	double p=0.5/1000, S_p=0.005/1000; //m, DA STIMARE
	double d=65.5/1000, S_d=0.5/1000; //m
	
	vector<double> rad_msper; vector<double> Srad_msper;
	
	for(int i=0; i<K; i++){
		rad_msper.push_back(p*corr_msper[i]/(400*d));
		Srad_msper.push_back(sqrt(pow(Scorr_msper[i]*p/(400*d), 2) + pow(corr_msper[i], 2)*(pow(S_p/400, 2) + pow(p*S_d/(400*pow(d, 2)), 2))));
	}
	
	cout << endl << "Conversione in radianti:" << endl;
	for(int i=0; i<K; i++){
		cout << rad_msper[i] <</* " " << Srad_msper[i] <<*/ endl;
	}
	cout << endl;
	for(auto i:Srad_msper) cout << i << endl;
	
	
	
//CALCOLO SPESSORE

	vector<double> ascisse; vector<double> ordinate_sper;
	
	for(int i=-3; i<=3; i++){
		if(i!=0){
			ascisse.push_back(i);
		}
	}
	int N = ascisse.size();
	
	for (int i=0; i<rad_msper.size(); i++){
			ordinate_sper.push_back(sin(rad_msper[i]));
	}
	cout << ordinate_sper.size() << endl;
	
	cout << "I valori dei seni con gli ordini: " << endl;
	
	
	for (int i=0; i<N; i++){
		cout << ascisse[i] /*<< " " << ordinate_sper[i]*/ << endl;
	}
	cout << endl;
	for (int i=0; i<N; i++){
		cout << ordinate_sper[i] << endl;
	}

//Interpolazione Caso 2 SPER


			double q1=0, q2=0, q3=0, q4=0, q5=0;
			for (int i=0; i<N; i++){
				q1 += pow(1/Srad_msper[i],2); // 1 / S_yi^2
				q2 += pow(ascisse[i]/Srad_msper[i], 2);// x_i^2 / s_yi^2
				q3 += ascisse[i]/pow(Srad_msper[i], 2), 2; // x_i / S_yi^2
				q4 += (ascisse[i] * ordinate_sper[i]) / pow(Srad_msper[i], 2); // (x_i * y_i) / S_yi^2
				q5 += ordinate_sper[i]/pow(Srad_msper[i], 2); // y_i / S_yi^2
			}
			
			double delta_sper = q1*q2 - q3*q3;
			double b_sper = 0.5*(q1*q4 - q3*q5)/delta_sper;
			double a_sper = (q2*q5 - q3*q4)/delta_sper;
			
			double Sa_sper = sqrt((1/delta_sper) * q2);
			double Sb_sper = sqrt((1/delta_sper) * q1);
			
			cout << "Parametri per i dati sper sono: " << endl << "Delta: " << delta_sper << " b= " << b_sper << " +- " << Sb_sper << " a= " << a_sper << " +- " << Sa_sper << endl;
			
			
			
			//Errore a posteriori e Coeff di correlazione lineare e Chi-Quadro

		
			double r1=0, media_y=0, media_x=0; double Sxy_sper = 0;
			for (int i=0; i<N; i++){
				r1 += pow(ordinate_sper[i] - (a_sper + b_sper*ascisse[i]), 2);
				media_y += ordinate_sper[i];
				media_x += ascisse[i];
			}
		
			
			double post = sqrt((r1)/(N-2));
			cout << "Errore a posteriori: " << post << endl;
			
			
			media_y = media_y/N;
			media_x = media_x/N;

			double w1=0, w2=0;
			for(int i=0; i<N; i++){
				Sxy_sper += (ascisse[i] - media_x) * (ordinate_sper[i] - media_y);
				w1 += pow(ascisse[i] - media_x, 2);
				w2 += pow(ordinate_sper[i] - media_y, 2);
			}
			Sxy_sper = Sxy_sper / (N-1);
			w1 = sqrt(w1/(N-1));
			w2 = sqrt(w2/(N-1));
			
			cout << "Medie: " << media_x << " +- " << w1 << " / " << media_y << " +- " << w2 << endl;
			cout << "Covarianza: " << Sxy_sper << endl; 


			double rxy=0;
			rxy = Sxy_sper / (w1*w2);
			double S_rxy = sqrt((1-pow(rxy, 2)) / (N-2));
			cout << "Coefficente di correlazione lineare: " << rxy << " +- " << S_rxy << endl;


			double X2=0;
			for (int i=0; i<N; i++){
				double Xi = pow((ordinate_sper[i] - (ascisse[i]*b_sper + a_sper))/ (Srad_msper[i]), 2) ;
				X2 += Xi;
			}   
			
			double NDOF = N - 2;
			cout << "Il chi-quadro e': " << X2 << " e il NDOF: " << NDOF << endl;
			cout << "-------------------------------------------------------" << endl << endl;
	
	
	
	double l = 670/*/(1000*1000*1000)*/, S_l= 1/*/(1000*1000*1000)*/;
	
	
	double d_sper = l/b_sper;
	double Sd_sper = (d_sper) * sqrt(pow(S_l/l, 2) + pow(Sb_sper/b_sper, 2));
	
	
	cout << "Spessore sperimentale: " << d_sper << " +- " << Sd_sper << endl;
	
	return 0;
}




parametri interpolazione_lin(vector<double> ascisse, vector<double> ordinate, double S_x, double S_y, int N){
	//Interpolazione Caso 2


			double q1=0, q2=0, q3=0, q4=0, q5=0;
			for (int i=0; i<N; i++){
				q1 += pow(1/S_y,2); // 1 / S_yi^2
				q2 += pow(ascisse[i]/S_y, 2);// x_i^2 / s_yi^2
				q3 += ascisse[i]/pow(S_y, 2), 2; // x_i / S_yi^2
				q4 += (ascisse[i] * ordinate[i]) / pow(S_y, 2); // (x_i * y_i) / S_yi^2
				q5 += ordinate[i]/pow(S_y, 2); // y_i / S_yi^2
			}
			
			double delta0 = q1*q2 - q3*q3;
			double b0 = (q1*q4 - q3*q5)/delta0;
			double a0 = (q2*q5 - q3*q4)/delta0;
			
			cout << "Parametri secondo il secondo caso:" << endl << "Delta: " << delta0 << " b= " << b0 << " a= " << a0 << endl;
			
	//Interpolazione caso 3

			double S_i = 0;
			
			double p1=0, p2=0, p3=0, p4=0, p5=0;
			for (int i=0; i<N; i++){
				S_i = sqrt(pow(S_y, 2) + pow(b0, 2) * pow(S_x, 2));
				p1 += pow(1/S_i,2); // 1 / S_i^2
				p2 += pow(ascisse[i]/S_i, 2);// x_i^2 / s_i^2
				p3 += ascisse[i]/pow(S_i, 2), 2; // x_i / S_i^2
				p4 += (ascisse[i] * ordinate[i]) / pow(S_i, 2); // (x_i * y_i) / S_yi^2
				p5 += ordinate[i]/pow(S_i, 2); // y_i / S_yi^2
			}
			
			double delta = p1*p2 - p3*p3;
			double b = (p1*p4 - p3*p5)/delta;
			double a = (p2*p5 - p3*p4)/delta;
			
			double S_a = sqrt(p2/delta);
			double S_b = sqrt(p1/delta);
			
			cout << "Parametri secondo il terzo caso: " << endl << " b= " << b << " +- " << S_b << " a = " << a << " +- " << S_a << endl;

			
			
			
		
		
	//Errore a posteriori e Coeff di correlazione lineare e Chi-Quadro

		
			double r1=0, media_y=0, media_x=0, Sxy = 0;
			for (int i=0; i<N; i++){
				r1 += pow(ordinate[i] - (a + b*ascisse[i]), 2);
				media_y += ordinate[i];
				media_x += ascisse[i];
			}
		
			
			double post = sqrt((r1)/(N-2));
			cout << "Errore a posteriori: " << post << endl;
			
			
			media_y = media_y/N;
			media_x = media_x/N;

			double w1=0, w2;
			for(int i=0; i<N; i++){
				Sxy += (ascisse[i] - media_x) * (ordinate[i] - media_y);
				w1 += pow(ascisse[i] - media_x, 2);
				w2 += pow(ordinate[i] - media_y, 2);
			}
			Sxy = Sxy / (N-1);
			w1 = sqrt(w1/(N-1));
			w2 = sqrt(w2/(N-1));
			
			cout << "Medie: " << media_x << " +- " << w1 << " / " << media_y << " +- " << w2 << endl;
			cout << "Covarianza: " << Sxy << endl; 

			double rxy = Sxy / (w1*w2);
			double S_rxy = sqrt((1-pow(rxy, 2)) / (N-2));
			cout << "Coefficente di correlazione lineare: " << rxy << " +- " << S_rxy << endl;


			double X2=0;
			for (int i=0; i<N; i++){
				double Xi = pow((ordinate[i] - (ascisse[i]*b + a))/ (S_y), 2) ;
				X2 += Xi;
			}   
			
			double NDOF = N - 2;
			cout << "Il chi-quadro e': " << X2 << " e il NDOF: " << NDOF << endl;
			
			parametri par;
			
			par.a = a;
			par.b = b;
			
			par.S_a = S_a;
			par.S_b = S_b;
			
			par.Sxy = Sxy;

			
			cout << "-------------------------------------------------------" << endl << endl;
			
			return par;
}