//Non posso far partire la simulazione con tutti i fermioni nelle stesse posizioni!

#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

using namespace std;

const int N = 4;
const int N_part = 8;
const int M = 100;
const int M_eq = 20;
double d = 0.4;
double L = 2.;

struct pos{     	
	double x;
	double y;
	double z;
	double x_new;
	double y_new;
	double z_new;
};

struct mat_A{
	double s;
	double px;
	double py;
	double pz;
	double inv[N];
	double dA_dx;
	double dA_dy;
	double dA_dz;
	double Bx[N];
	double By[N];
	double Bz[N];
};

pos *R;
mat_A *A;


void start();
double mc_step(double x);
void mc_update();
void fill_A(double n, int set);
void fill_A_new(double n, int set);
double compute_det();
void acc_rej(double P, double *p_acc);
double p(double det, double det_new);
void compute_inverse();
double log_der(int I, int J, double n);


int main(){
	double det=0., det_new=0., P=0., p_acc=0., nabla=0., n=1.; 
	srand((unsigned int)time(NULL));
	R = new pos[N_part];
	A = new mat_A[N];
	
	//********* EQUILIBRATION: randomizes the positions ***************
	start();
	for(int i=0; i<M_eq; i++){
		mc_update();
		acc_rej(1, &p_acc);
	}
	p_acc = 0.;
	
	//************ MONTE CARLO ***********************
	for(int i=0; i<M; i++){
		mc_update();
		
		//SPIN-UP DETERMINANT
	    fill_A(n, 0);
	    det = compute_det();
	    fill_A_new(n, 0);
	    det_new = compute_det();
	    P = p(det, det_new);
	    
	    //SPIN-DOWN DETERMINANT
	    fill_A(n, 4);
	    det = compute_det();
	    fill_A_new(n, 4);
	    det_new = compute_det();
	    
	    //ACCEPTANCE/REJECTION TEST	    
	    P = P*p(det, det_new);
	    acc_rej(P, &p_acc);
	    
	    //ENERGY CALCULATION
	    fill_A(n, 0);
	}
	cout << p_acc/M << endl;
}


//********* RANDOM STARTING POSITIONS ***************
void start(){
	for(int j=0; j<N_part; j++){
		R[j].x = (double)rand()/RAND_MAX;
		R[j].y = (double)rand()/RAND_MAX;
		R[j].z = (double)rand()/RAND_MAX;
	}
}


//********* MC FUNCTIONS ***************
double mc_step(double x){
	double eta=0., x_new=0.;
    eta = (double)rand()/RAND_MAX;
	x_new = x + d*(eta - 0.5);
	return x_new;
}


void mc_update(){
	for(int j=0; j<N_part; j++){
		R[j].x_new = mc_step(R[j].x);
		R[j].y_new = mc_step(R[j].y);
		R[j].z_new = mc_step(R[j].z);
	}
}


//*************** ACCEPTANCE TEST ******************
double p(double det, double det_new){
   double P;
   P = det_new*det_new/(det*det);
   return P;
}


void acc_rej(double P, double *p_acc){
	double w = (double)rand()/RAND_MAX;
	if(P > w){
		for(int j=0; j<N_part; j++){
			R[j].x = R[j].x_new;
			R[j].y = R[j].y_new;
			R[j].z = R[j].z_new;
		}
		*p_acc += 1.; 
	}
}


//*************** MATRIX FILLING ******************
//"set" is either 0 or 4 and selects the "up" set and the "down" set of coordinates
void fill_A(double n, int set){
	double r;
	for(int i=set; i<N+set; i++){
		r = sqrt(R[i].x*R[i].x + R[i].y*R[i].y + R[i].z*R[i].z);
	    A[i-set].s = exp(-r*r/(2.*n*n));
	    A[i-set].px = R[i].x*exp(-r*r/(2.*n*n));
	    A[i-set].py = R[i].y*exp(-r*r/(2.*n*n));
	    A[i-set].pz = R[i].z*exp(-r*r/(2.*n*n));
	}
}

void fill_A_new(double n, int set){
	double r;
	for(int i=set; i<N+set; i++){
		r = sqrt(R[i].x_new*R[i].x_new + R[i].y_new*R[i].y_new + R[i].z_new*R[i].z_new);
	    A[i-set].s = exp(-r*r/(2.*n*n));
	    A[i-set].px = R[i].x_new*exp(-r*r/(2.*n*n));
	    A[i-set].py = R[i].y_new*exp(-r*r/(2.*n*n));
	    A[i-set].pz = R[i].z_new*exp(-r*r/(2.*n*n));
	}
}


//*************** DETERMINANT ******************
double compute_det(){
   int s;
   double det_A;
   gsl_matrix *a = gsl_matrix_calloc(N, N);
   for(int i=0; i<N; i++){
   	  gsl_matrix_set(a, i, 0, A[i].s);
   	  gsl_matrix_set(a, i, 1, A[i].px);
   	  gsl_matrix_set(a, i, 2, A[i].py);
   	  gsl_matrix_set(a, i, 3, A[i].pz);
   }
   gsl_permutation * p = gsl_permutation_alloc(N);
   gsl_linalg_LU_decomp(a, p, &s);
   det_A = gsl_linalg_LU_det(a, s);
   return det_A;
   gsl_matrix_free(a);
   gsl_permutation_free(p);
}


//*************** INVERSE MATRIX ******************
void compute_inverse(){
   int s;
   gsl_matrix *a = gsl_matrix_alloc(N, N);
   for(int i=0; i<N; i++){
   	  gsl_matrix_set(a, i, 0, A[i].s);
   	  gsl_matrix_set(a, i, 1, A[i].px);
   	  gsl_matrix_set(a, i, 2, A[i].py);
   	  gsl_matrix_set(a, i, 3, A[i].pz);
   }
   gsl_permutation * p = gsl_permutation_alloc(N);
   gsl_linalg_LU_decomp(a, p, &s);
   gsl_matrix *inv = gsl_matrix_alloc(N, N);
   gsl_linalg_LU_invert(a, p, inv);
   gsl_matrix_free(a);
   gsl_permutation_free(p);
    for(int i=0; i<N; i++){
    	for(int j=0; j<N; j++){
    		A[i].inv[j] = gsl_matrix_get(inv, i, j);
    	}
   }
   gsl_matrix_free(inv);
}

//*********** LOGARITHMIC DERIVATIVE OF THE DETERMINANT *************
double log_der(int I, int J, double n){
	double r = sqrt(R[j].x*R[j].x + R[j].y*R[j].y + R[j].z*R[j].z);
	double R_vec = {R[j].x, R[j].y, R[j].z};
	double Tr=0.;
	for(k=0; k<N; k++){
		Tr -= A[k].inv[j]*R_vec[j]*R_vec[k];
	}
	Tr = Tr + A[i].inv[j]*n*n;
	Tr = Tr*exp(-r*r/(2.*n*n))/(n*n);
	return Tr;
} //metti a posto qua


void fill_nablaA(int J){
	dA_dx[] = -exp(-r*r/(2.*n*n))/(n*n)*{-R[J].x, R[J].x*R[J].x-1., R[J].x*R[J].y, R[J].x*R[J].z};
	dA_dy[] = -exp(-r*r/(2.*n*n))/(n*n)*{-R[J].y, R[J].x*R[J].y, R[J].y*R[J].y-1., R[J].y*R[J].z};
	dA_dx[] = -exp(-r*r/(2.*n*n))/(n*n)*{-R[J].z, R[J].x*R[J].z, R[J].z*R[J].y, R[J].z*R[J].z-1.};
}


void compute_B(int J, int n){
	fill_nablaA(J);
	double r = sqrt(R[J].x*R[J].x + R[J].y*R[J].y + R[J].z*R[J].z);	
	for(int l=0; l<N; l++){
		for(int k=0; k<N; k++){
			A[l].Bx[k] = A[l].inv[J]*dA_dx[k];
			A[l].By[k] = A[l].inv[J]*dA_dy[k];
			A[l].Bz[k] = A[l].inv[J]*dA_dz[k];
		}
	}
}




