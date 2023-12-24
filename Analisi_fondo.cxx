#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <dirent.h>

using namespace std;

int main(){
	
	vector<double> intensita; vector<double> theta;
	double S_y=2.3, S_x=1; //pm DA STIMARE
	
	
	double t_i, I_i;
	while(cin >> t_i >> I_i){
		intensita.push_back(I_i);
		theta.push_back(t_i);
	}
	int M = theta.size();

	vector<double> massimi;

	for(int i=0; i<M; i++){
		if (intensita[i-3]<intensita[i] && intensita[i+3]<intensita[i]){
			if (intensita[i-2]<intensita[i] && intensita[i+2]<intensita[i]){
					massimi.push_back(intensita[i]);
			}
		}
	}
	
	int N = massimi.size();
	cout << N << endl << endl;
	
	double m=0, Sm=S_y/sqrt(N);
	for (auto i: massimi){
		m += i;
	}
	
	cout << "La media dei massimi e' " << m/N << " +- " << Sm << endl;
	
	//for(auto i:massimi) cout << i << endl;
	
	double media=0, Smedia;
	for (auto i: intensita){
		media += i;
	}
	media = media/M; Smedia = S_y / sqrt(M);
	
	cout << "La media totale e' " << media << " +- " << Smedia << endl;

	return 0;
}