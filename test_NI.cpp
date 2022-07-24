//Aggiorno una posizione alla volta 
//Posso verficare se il rapporto di determinanti calcolati con S è lo stesso che si trova con la funzione per il calcolo del det. 
//attento perché se la mossa viene riegettata devi annullare l'aggiornameno della matrice inversa: questo lo fai tramite l'if che hai messo

#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

using namespace std;

const int N = 4;
const int N_part = 8;
const int M = 100;
double d = 0.5;

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
};

pos *R;
mat_A *A_up;
mat_A *A_down;

void start();
void fill_A(double n, mat_A *A, int I);
void mc_update(int I);
double mc_step(double x);
void update_A(double n, mat_A *A, int I);
void compute_inverse(mat_A* A);
double compute_det(mat_A *A);
double compute_S(mat_A *A, int I);
double compute_Z_j(mat_A *A, int I, int j);
void update_inv(mat_A *A, int I, double S);
void acc_rej(double P, double *p_acc, int I);


int main(){
	R = new pos[N_part];
	A_up = new mat_A[N];
	A_down = new mat_A[N];
	double n = 1., p_acc, tmp, S_up=0., S_down=0.; 
	srand((unsigned int)time(NULL));
	start();
	for(int j=0; j<N; j++){
		fill_A(n, A_up, j);
	  fill_A(n, A_down, j + N);
	}
	compute_inverse(A_up);
	compute_inverse(A_down);
	
   for(int k=0; k<M; k++){
	//SPIN UP PART
	S_down = compute_S(A_down, 0); //qua non conta quale I prendo perché tanto la A_down non è aggiornata.  
	for(int I=0; I<N; I++){
		mc_update(I);//aggiorna la posizione di una particella
		update_A(n, A_up, I);//aggiorna la riga i-esima di A
		S_up = compute_S(A_up, I); //è il rapporto det(new)/det(old)
		tmp = p_acc;
	  acc_rej(S_up*S_up*S_down*S_down, &p_acc, I); //accetta o rigetta l'aggiornamento
		if(tmp < p_acc){
			update_inv(A_up, I, S_up); //aggiorna l'inversa in caso
		}else{
			fill_A(n, A_up, I); //torna indietro alla configurazione precedente se la mossa è rigettata
		}
	}
	
	//SPIN DOWN PART
	S_up = compute_S(A_up, 0); 
	for(int I=0; I<N; I++){
		mc_update(I + N);
		update_A(n, A_down, I);
		S_up = compute_S(A_up, I); 
		tmp = p_acc;
	  acc_rej(S_up*S_up*S_down*S_down, &p_acc, I);
		if(tmp < p_acc){
			update_inv(A_down, I, S_down);
		}else{
			fill_A(n, A_down, I + N);
		}
	}
	}//end cycle over M
	
   cout << "p_acc : "<< p_acc/(M*8.) << endl;
	cout << compute_det(A_up) << "     " << compute_det(A_down) << endl;
}


//********* RANDOM STARTING POSITIONS ***************
void start(){
	for(int j=0; j<N_part; j++){
		R[j].x = (double)rand()/RAND_MAX;
		R[j].y = (double)rand()/RAND_MAX;
		R[j].z = (double)rand()/RAND_MAX;
	}
}

void fill_A(double n, mat_A *A, int I){
	double r;
	int q = 0;
	r = sqrt(R[I].x*R[I].x + R[I].y*R[I].y + R[I].z*R[I].z);
	if(I > 3){
		q = 4; 
	}
	A[I - q].s = exp(-r*r/(2.*n*n));
	A[I - q].px = R[I].x*exp(-r*r/(2.*n*n));
	A[I - q].py = R[I].y*exp(-r*r/(2.*n*n));
	A[I - q].pz = R[I].z*exp(-r*r/(2.*n*n));
}


//***************** INVERSE MATRIX ********************
void compute_inverse(mat_A* A){
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

//***************** DETERMINANT **********************
double compute_det(mat_A *A){
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


//************** MC FUNCTIONS ***********************
void mc_update(int I){
	R[I].x_new = mc_step(R[I].x);
	R[I].y_new = mc_step(R[I].y);
	R[I].z_new = mc_step(R[I].z);
}

double mc_step(double x){
	double eta=0., x_new=0.;
    eta = (double)rand()/RAND_MAX;
	x_new = x + d*(eta - 0.5);
	return x_new;
}

void update_A(double n, mat_A *A, int I){
	double r = sqrt(R[I].x_new*R[I].x_new + R[I].y_new*R[I].y_new + R[I].z_new*R[I].z_new);
	A[I].s = exp(-r*r/(2.*n*n));
	A[I].px = R[I].x_new*exp(-r*r/(2.*n*n));
	A[I].py = R[I].y_new*exp(-r*r/(2.*n*n));
	A[I].pz = R[I].z_new*exp(-r*r/(2.*n*n));
}

//*************** UPDATE INVERSE MATRIX ******************
double compute_S(mat_A *A, int I){
	double S = 0.;
	S = A[I].s*A[0].inv[I] + A[I].px*A[1].inv[I] + A[I].py*A[2].inv[I] + A[I].pz*A[3].inv[I];
	return S;
}

double compute_Z_j(mat_A *A, int I, int j){
	double Z_j = 0.;
	Z_j = A[I].s*A[0].inv[j] + A[I].px*A[1].inv[j] + A[I].py*A[2].inv[j] + A[I].pz*A[3].inv[j];
	return Z_j;
}

void update_inv(mat_A *A, int I, double S){
	for(int k=0; k<N; k++){
		for(int j=0; j<N; j++){
			if(j != I){
				A[k].inv[j] = A[k].inv[j] - A[k].inv[I]*compute_Z_j(A, I, j)/S;
			}
		}
		A[k].inv[I] = A[k].inv[I]/S;
	}
}

//***************** ACCEPTANCE TEST **********************
void acc_rej(double P, double *p_acc, int I){
	double w = (double)rand()/RAND_MAX;
	if(P > w){
		R[I].x = R[I].x_new;
	    R[I].y = R[I].y_new;
	    R[I].z = R[I].z_new;
		*p_acc += 1.; 
	}
}

