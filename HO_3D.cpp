
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

float  k(float, float, int);
float numL(float *, float *, float, float, int, int);
float numR(float *, float *, float, float, int, int);
void start_L(float *, float *, float, int, int);
void start_R(float *, float *, float, int);
void print_Psi(int, int, float *, float *, int);

int main(){
	float h, eps=0.1;
	int N = 1400, M = 1000, i, k, l;
	float r_max = 7.;
    h = r_max/N;
	float *r, *y;
    float E_max = 10., E_min;
	float dE;
	float *F;
	F = new float[M];
	for(l=0; l<7; l++){
		k = 0;
	    E_min = sqrt(l*(l+1));
	    dE = (float)(E_max - E_min)/M;
	    ofstream myfile;
	    string title;
        title = "F(E)_l=" + to_string(l) + ".txt";
        myfile.open(title, ios :: out | ios :: trunc);
        myfile.precision(6);
    
    for(i=1; i<M-1; i++){
    	r = new float[N];
        y = new float[N];
        start_L(r, y, h, N, l);
        F[i] = numL(r, y, h, E_min + i*dE, N, l);
        start_R(r, y, h, N);
        F[i] = F[i] - numR(r, y, h, E_min + i*dE, N, l); 
        myfile << E_min + i*dE << "       " << F[i] << endl;
        if(F[i]*F[i-1]<0 && fabs(F[i])<eps){
             	print_Psi(k, N, r, y, l);
             	k++;
           }
        }
	myfile.close();
	}
}

//Functions 

void print_Psi(int k, int N, float *r, float *y, int l){
	int j;
	ofstream os;
    string Psi;
    Psi = "n=" + to_string(2*k) + "_" + "l=" + to_string(l) + ".txt";
    os.open(Psi, ios :: out | ios :: trunc);
    os.precision(6);
    float c = 1.0;
    for(j=0; j<N; j++){
        if(j != 0 && fabs(y[j]-y[j-1])>0.1){
           c = y[j-1]/y[j];
        }
     os << r[j] << "       " << c*y[j] << endl;
     }
    os.close();
}

float k(float r, float E, int l){
	return 2.*E - r*r - l*(l+1)/(float)(r*r);
}

void start_L(float *r, float *y, float h, int N, int l){
	r[0] = h/4.; 
    r[1] = h;
    y[0] = pow(r[0], l+1);
    y[1] = pow(r[1], l+1);
}

void start_R(float *r, float *y, float h, int N){
    y[N-1] = 0;
    y[N-2] = pow(10, -10);
    r[N-1] = h*N;
    r[N-2] = h*N - h;
}


float numL(float *r, float *y, float h, float E, int N, int l){
	int i;
	float g1=0.0, g2=0.0, g3=0.0;
	for(i=1; i<N-1; i++){
		r[i+1] = r[i] + h;
		g1 = h*h*5*k(r[i], E, l)/12.;
		g2 = h*h*k(r[i-1], E, l)/12.;
		g3 = h*h*k(r[i+1], E, l)/12.;
		y[i+1] = 2.*(1 - g1)*y[i] - (1 + g2)*y[i-1];
		y[i+1] = y[i+1]/(1 + g3);
		if(k(r[i], E, l)*k(r[i+1], E, l)<0 && r[i]> sqrt(E - sqrt(E*E - l*(l+1)))){
			break;
		}
	}
	return (y[i+1] - y[i-1])/(2.*h*y[i]);
}

 
float numR(float *r, float *y, float h, float E, int N, int l){
	int i;
	float g1=0.0, g2=0.0, g3=0.0;
	for(i=N-2; i>=1; i--){
		r[i-1] = r[i] - h;
		g1 = h*h*5*k(r[i], E, l)/12.;
		g2 = h*h*k(r[i-1], E, l)/12.;
		g3 = h*h*k(r[i+1], E, l)/12.;
		y[i-1] = 2.*(1 - g1)*y[i] - (1 + g3)*y[i+1];
		y[i-1] = y[i-1]/(1 + g2);
		if(k(r[i], E, l)*k(r[i-1], E, l)<0){
			break;
		}
	}
	return (y[i+1] - y[i-1])/(2.*h*y[i]);
}
