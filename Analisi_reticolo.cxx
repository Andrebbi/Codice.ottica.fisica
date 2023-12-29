#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <dirent.h>

//Per girare il programma basta lanciarlo da terminale e scrivere il nome della cartella contenente i file di dati corrispondente all'analisi

using namespace std;

struct parametri {
	double a=0;
	double b=0;
	double S_a=0;
	double S_b=0;
	
	double Sxy=0;
};


parametri interpolazione_lin(vector<double>, vector<double>, double, int);


int main(){
	
	
	string folder_path;
    cout << "Inserisci il percorso della cartella: ";
    getline(cin, folder_path);
    
    // Controllo se l'ultimo carattere della stringa folder_path è una "/".
    // Se non lo è, aggiungo la "/" alla fine della stringa.
    if (folder_path.back() != '/') {
        folder_path += '/';
    }
    
    DIR* directory = opendir(folder_path.c_str());
    if (directory == NULL) {
        cout << "Errore: Impossibile aprire la cartella." << endl;
        return 1;
    }
	
	double pi = 3.14159265359;
	double S_t = 2/sqrt(3); //primi
	S_t = S_t/60; //gradi		
	//S_t = S_t *(pi/180); //rad
	
	double theta0 = 107.1167, S_theta0 = S_t/sqrt(2), S_theta = sqrt(2* pow(S_theta0, 2)); //gradi
	cout << S_theta0 << " " << S_theta << endl;
	double S_rad = S_theta*(pi/180); 
	cout << S_rad << endl;
	double d = 3.33, S_d = 0.03; // um
	
	vector<double> lambda11; vector<double> S_lambda11; 
	vector<double> lambda22; vector<double> S_lambda22; 
	
	struct dirent* file;
    while ((file = readdir(directory)) != NULL) {
        //evito le directory mamme e figlie
        if (strcmp(file->d_name, ".") == 0 || strcmp(file->d_name, "..") == 0) {
            continue;
        }
        string file_path = folder_path + file->d_name; //crea percorso completo file
        ifstream input_file(file_path.c_str()); //apre il file
        if (!input_file.is_open()) {
            cout << "Errore: Impossibile aprire il file " << file_path << endl;
            continue;
        }
		cout << endl;
		
		
	//Presa dei dati
		vector<double> gradi_A; vector<double> primi_A;
		vector<double> gradi_B; vector<double> primi_B;
		vector<double> ordini;
		
		
		double gA, gB, pA, pB, o;
		while (input_file >> gA >> pA >> gB >> pB >> o){
			gradi_A.push_back(gA);
			primi_A.push_back(pA);
			gradi_B.push_back(gB);
			primi_B.push_back(pB);
			ordini.push_back(o);
		}
		
		int N = ordini.size();
		
		cout << "Dati corretti: " << endl;
		vector<double> gradi; vector<double> radianti;
		for (int i=0; i<N; i++){
			gradi.push_back((gradi_A[i] + primi_A[i]/60 + gradi_B[i] + primi_B[i]/60 - 180)/2 - theta0);
			radianti.push_back(gradi[i]*(pi/180));
			cout << gradi[i] << endl;
		}
		
		cout << endl << "Dati in radianti: " << endl;
		for(auto i:radianti) cout << i << endl;


	//Primo metodo
		vector<double> lambda1; vector<double> S_lambda1;
		
		for (int i=0; i<N; i++){
			lambda1.push_back(-d*sin(radianti[i])/ordini[i]);
			S_lambda1.push_back(sqrt(pow(sin(radianti[i])*S_d/ordini[i], 2) + pow(d*cos(radianti[i])*S_rad/ordini[i], 2)));
			
			cout << lambda1[i] << endl;
		}
		cout << endl;
		for(auto i:S_lambda1) cout << i << endl;
		
		double max=0, S_max, min=1, S_min;
		for(int i=0; i<N; i++){
			if (lambda1[i] > max){
				max = lambda1[i];
				S_max = S_lambda1[i];
			}
			if (lambda1[i] < min){
				min = lambda1[i];
				S_min = S_lambda1[i];
			}
		}
		cout << endl << max << " +- " << S_max << " " << min << " +- " << S_min << endl;
		
		double compatibilita = (max - min)/sqrt(pow(S_max, 2) + pow(S_min, 2));
		cout << "La compatibilita' e' " << compatibilita << endl;
		
		
		double mediap=0, k=0;
		for (int i=0; i<N; i++){
			k += 1 / pow(S_lambda1[i], 2);
			mediap += lambda1[i] / pow(S_lambda1[i], 2);
		}
		mediap = mediap / k;
		double S_mediap = sqrt(1 / k);
		
		cout << "Lambda secondo il primo metodo = " << mediap  << " +- " << S_mediap << endl << endl;
		
		lambda11.push_back(mediap);
		S_lambda11.push_back(S_mediap);
		
		
	//Secondo metodo
	
		vector<double> seni;
		
		for(int i=0; i<N; i++){
			seni.push_back(sin(radianti[i]));
		}
		cout << "I seni: " << endl;
		for(int i=0; i<N; i++){
			cout << seni[i] << endl;
		}
		
		
		parametri par;
		
		par = interpolazione_lin(ordini, seni, S_rad, N);
		
		double lambda2 = -d * par.b;
		double S_lambda2 = lambda2 * sqrt(pow(par.S_b*d, 2) + pow(S_d * par.b, 2));
		cout << "Lambda secondo il secondo metodo e' = " << lambda2  << " +- " << S_lambda2 << endl;
		
		lambda22.push_back(lambda2);
		S_lambda22.push_back(S_lambda2);
		
		double comp = abs(par.a)/par.S_a;
		cout << "La compatibilita' tra a e 0 e': " << comp << endl << endl;
	}
	
	
	cout << "Le compatibilita' tra le lambda del primo e secondo metodo sono: " << endl;
	for (int i=0; i<lambda11.size(); i++){
		double comp = abs(lambda11[i] - lambda22[i])/sqrt(pow(S_lambda11[i], 2) + pow(S_lambda22[i], 2));
		cout << comp << endl;
	}
	
	
	return 0;
}


parametri interpolazione_lin(vector<double> ascisse, vector<double> ordinate, double S_y, int N){
	//Interpolazione Caso 2


			double q1=0, q2=0, q3=0, q4=0, q5=0;
			for (int i=0; i<N; i++){
				q1 += pow(1/S_y,2); // 1 / S_yi^2
				q2 += pow(ascisse[i]/S_y, 2);// x_i^2 / s_yi^2
				q3 += ascisse[i]/pow(S_y, 2), 2; // x_i / S_yi^2
				q4 += (ascisse[i] * ordinate[i]) / pow(S_y, 2); // (x_i * y_i) / S_yi^2
				q5 += ordinate[i]/pow(S_y, 2); // y_i / S_yi^2
			}
			
			double delta = q1*q2 - q3*q3;
			double b = (q1*q4 - q3*q5)/delta;
			double a = (q2*q5 - q3*q4)/delta;
			
			double S_a = sqrt((1/delta) * q2);
			double S_b = sqrt((1/delta) * q1);
			
			cout << "Parametri secondo il secondo caso:" << endl << "Delta: " << delta << " b= " << b << " +- " << S_b << " a= " << a << " +- " << S_a << endl;
	
	
	
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
			
			double t = abs(rxy)/S_rxy;
			cout << "Lo student e': " << t << endl;
			
			parametri par;
			
			par.a = a;
			par.b = b;
			
			par.S_a = S_a;
			par.S_b = S_b;
			
			par.Sxy = Sxy;

			
			cout << "-------------------------------------------------------" << endl << endl;
			
			return par;
}