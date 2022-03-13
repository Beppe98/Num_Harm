
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

float  k(float, float);
void numL(float *, float *, float, float, int);
void numR(float *, float *, float, float, int);

int main(){
	float E=1., h;
	int N = 1400;
	float x_min = -3.5;
    float x_max = 3.5;
    float log_der_L = 0.;
    float log_der_R = 0.;
    h = (x_max - x_min)/N;
	float *x, *y;
    x = new float[N];
    y = new float[N];
    y[0] = 0;
    y[1] = 0.01;
    x[0] = -h*N/2.;
    x[1] = -h*N/2 + h;
    numL(x, y, h, E, N);
    x = new float[N];
    y = new float[N];
    y[N-1] = 0;
    y[N-2] = 0.01;
    x[N-1] = h*N/2.;
    x[N-2] = h*N/2 - h;
    numR(x, y, h, E, N);
}

//Functions 

float k(float x, float E){
	return 2.*E - x*x;;
}

void numL(float *x, float *y, float h, float E, int N){
	int i;
	ofstream file1;
    file1.open("Num_L.txt", ios :: out | ios :: trunc);
    file1.precision(6);
    file1 << x[0] << "       " << y[0] << endl;
    file1 << x[1] << "       " << y[1] << endl;
	for(i=1; i<N/2-1; i++){
		x[i+1] = x[i] + h;
		y[i+1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i-1], E)/12.)*y[i-1];
		y[i+1] = y[i+1]/(1 + h*h*k(x[i+1], E)/12.);
		file1 << x[i+1] << "       " << y[i+1] << endl;
	}
	for(i=N/2-1; i<N-1; i++){
		x[i+1] = x[i] + h;
		y[i+1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i-1], E)/12.)*y[i-1];
		y[i+1] = y[i+1]/(1 + h*h*k(x[i+1], E)/12.);
		file1 << x[i+1] << "       " << y[i+1] << endl;
		if(k(x[i], E)*k(x[i+1], E)<0){
			break;
		}
	}
	file1.close();
}
 
void numR(float *x, float *y, float h, float E, int N){
	int i;
	ofstream file2;
    file2.open("Num_R.txt", ios :: out | ios :: trunc);
    file2.precision(6);
    file2 << x[N-1] << "       " << y[N-1] << endl;
    file2 << x[N-2] << "       " << y[N-2] << endl;
	for(i=N-2; i>=1; i--){
		x[i-1] = x[i] - h;
		y[i-1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i+1], E)/12.)*y[i+1];
		y[i-1] = y[i-1]/(1 + h*h*k(x[i-1], E)/12.);
		file2 << x[i-1] << "       " << y[i-1] << endl;
		if(k(x[i], E)*k(x[i-1], E)<0){
			break;
		}
	}
	file2.close();
}
 
