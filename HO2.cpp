//Harmonic oscillator: pay atttention: you need to impose that the wavefunction 
//goes to zero when x goes to +/- \infty. 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

float  k(float, float);
float numL(float *, float *, float, float, int);
float numR(float *, float *, float, float, int);

int main(){
	float E=0.5, h;
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
    log_der_L = numL(x, y, h, E, N);
    x = new float[N];
    y = new float[N];
    y[N-1] = 0;
    y[N-2] = 0.01;
    x[N-1] = h*N/2.;
    x[N-2] = h*N/2 - h;
    log_der_R = numR(x, y, h, E, N);
    cout << log_der_L - log_der_R << endl;
}

//Functions 

float k(float x, float E){
	return 2.*E - x*x;;
}

float numL(float *x, float *y, float h, float E, int N){
	int i;
	for(i=1; i<N/2-1; i++){
		x[i+1] = x[i] + h;
		y[i+1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i-1], E)/12.)*y[i-1];
		y[i+1] = y[i+1]/(1 + h*h*k(x[i+1], E)/12.);
	}
	for(i=N/2-1; i<N-1; i++){
		x[i+1] = x[i] + h;
		y[i+1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i-1], E)/12.)*y[i-1];
		y[i+1] = y[i+1]/(1 + h*h*k(x[i+1], E)/12.);
		if(k(x[i], E)*k(x[i+1], E)<0){
			break;
		}
	}
	return (y[i+1] - y[i-1])/(2.*h*y[i]);
}
 
float numR(float *x, float *y, float h, float E, int N){
	int i;
	for(i=N-2; i>=1; i--){
		x[i-1] = x[i] - h;
		y[i-1] = 2.*(1 - h*h*5*k(x[i], E)/12.)*y[i] - (1 + h*h*k(x[i+1], E)/12.)*y[i+1];
		y[i-1] = y[i-1]/(1 + h*h*k(x[i-1], E)/12.);
		if(k(x[i], E)*k(x[i-1], E)<0){
			break;
		}
	}
	return (y[i+1] - y[i-1])/(2.*h*y[i]);
}
