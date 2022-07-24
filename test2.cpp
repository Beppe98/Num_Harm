//Tentativo numero 2: I può essere solo 0, 1, 2, 3
//Forse il problema era nell'aggiornamento dell'inversa: gli elementi aggiornati di A^-1 venivano usati per gli aggiornamenti successivi 
//Lezione per la vita: meglio introdurre tante variabili
#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

using namespace std;

const int N = 4;
const int M = 1000;
double d = 1.2;

using namespace std;

struct pos{     	
	double x;
	double y;
	double z;
};

struct mat_A{
	double s;
	double px;
	double py;
	double pz;
    double inv[N];
};

struct der{
   double dA_dx;
   double dA_dy;
   double dA_dz;
   double d2A_dx2;
   double d2A_dy2;
   double d2A_dz2;
};

struct determinant{
   double nabla_x;
   double nabla_y;
   double nabla_z;
   double lap_x;
   double lap_y;
   double lap_z; 
};

pos *R_up;
pos *R_up_new;
pos *R_down;
pos *R_down_new;
mat_A *A_up;
mat_A *A_up_new;
mat_A *A_down;
mat_A *A_down_new;
der *D;
determinant *det_up;
determinant *det_down;

void start(pos *R);
void compute_inverse(mat_A* A);
void mc_update(pos *R_new, pos *R, int I);
double mc_step(double x);
void fill_A(double n, mat_A *A, pos *R, int I);
double compute_S(mat_A *A, mat_A *A_new, int I, int J);
void update_inv(mat_A *A, mat_A *A_new, int I, double S);
void acc_rej(double P, double *p_acc, pos *R, pos *R_new, int I);
void d_det(double n, int I, pos *R, mat_A *A, determinant *det);
void d2_det(double n, int I, pos *R, mat_A *A, determinant *det);
double energy(pos *R_a, pos *R_b, determinant *det_a, determinant *det_b);


int main(){
	R_up = new pos[N];
	R_up_new = new pos[N];
	R_down = new pos[N];
	R_down_new = new pos[N];
	A_up = new mat_A[N];
	A_up_new = new mat_A[N];
	A_down = new mat_A[N];
	A_down_new = new mat_A[N];
	double n, S_up, S_down, tmp, p_acc, E=0., E_tot=0., E2=0.;
	srand((unsigned int)time(NULL));
	start(R_up);
	start(R_down);
	for(int i=0; i<100; i++){
		p_acc=0., E=0., E_tot=0., E2=0.;
		n = 0.5 + i*0.02;
	for(int I=0; I<N; I++){
		fill_A(1., A_up, R_up, I);
		fill_A(1., A_down, R_down, I);
	}
	compute_inverse(A_up);
	compute_inverse(A_down);
	for(int j=0; j<M; j++){
		D = new der[N];
	    det_up = new determinant[N];
	    det_down = new determinant[N];
	//SPIN UP PART
	for(int I=0; I<N; I++){
		mc_update(R_up_new, R_up, I);
		fill_A(n, A_up_new, R_up_new, I);
		S_up = compute_S(A_up, A_up_new, I, I);
		tmp = p_acc;
		acc_rej(S_up*S_up, &p_acc, R_up, R_up_new, I);
		if(tmp < p_acc){
			update_inv(A_up, A_up_new, I, S_up);
			fill_A(n, A_up, R_up_new, I);
		}
	    d_det(n, I, R_up, A_up, det_up); 
	    d2_det(n, I, R_up, A_up, det_up);	
	}
		//SPIN DOWN PART
	for(int I=0; I<N; I++){
		mc_update(R_down_new, R_down, I);
		fill_A(n, A_down_new, R_down_new, I);
		S_down = compute_S(A_down, A_down_new, I, I);
		tmp = p_acc;
		acc_rej(S_down*S_down, &p_acc, R_down, R_down_new, I);
		if(tmp < p_acc){
			update_inv(A_down, A_down_new, I, S_down);
			fill_A(n, A_down, R_down_new, I);
		}	
		d_det(n, I, R_down, A_down, det_down); 
	    d2_det(n, I, R_down, A_down, det_down);	
	}
	E = energy(R_up, R_down, det_up, det_down)/M;
	E_tot += E;
	E2 += E*E*M;
	}
	cout << E_tot << " +/- " << sqrt((E2 - E_tot*E_tot)/M) << "      " << p_acc/(M*8.) << endl;	
	}
}

//********* RANDOM STARTING POSITIONS ***************
void start(pos *R){
	for(int j=0; j<N; j++){
		R[j].x = (double)rand()/RAND_MAX;
		R[j].y = (double)rand()/RAND_MAX;
		R[j].z = (double)rand()/RAND_MAX;
	}
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

//************** MC FUNCTIONS ***********************
void mc_update(pos *R_new, pos *R, int I){
	R_new[I].x = mc_step(R[I].x);
	R_new[I].y = mc_step(R[I].y);
	R_new[I].z = mc_step(R[I].z);
}

double mc_step(double x){
	double eta=0., x_new=0.;
    eta = (double)rand()/RAND_MAX;
	x_new = x + d*(eta - 0.5);
	return x_new;
}

void fill_A(double n, mat_A *A, pos *R, int I){
	double r;
	r = sqrt(R[I].x*R[I].x + R[I].y*R[I].y + R[I].z*R[I].z);
	A[I].s = exp(-r*r/(2.*n*n));
	A[I].px = R[I].x*exp(-r*r/(2.*n*n));
	A[I].py = R[I].y*exp(-r*r/(2.*n*n));
	A[I].pz = R[I].z*exp(-r*r/(2.*n*n));
}

//*************** UPDATE INVERSE MATRIX ******************
double compute_S(mat_A *A, mat_A *A_new, int I, int J){
	double S = 0.;
	S = A_new[I].s*A[0].inv[J] + A_new[I].px*A[1].inv[J] + A_new[I].py*A[2].inv[J] + A_new[I].pz*A[3].inv[J];
	return S;
}

void update_inv(mat_A *A, mat_A *A_new, int I, double S){
	for(int k=0; k<N; k++){
		for(int j=0; j<N; j++){
			if(j != I){
				A_new[k].inv[j] = A[k].inv[j] - A[k].inv[I]*compute_S(A, A_new, I, j)/S;
			}
		}
		A_new[k].inv[I] = A[k].inv[I]/S;
	}
	for(int k=0; k<N; k++){
		for(int j=0; j<N; j++){
			A[k].inv[j] = A_new[k].inv[j];
		}
	}
}

//***************** ACCEPTANCE TEST **********************
void acc_rej(double P, double *p_acc, pos *R, pos *R_new, int I){
	double w = (double)rand()/RAND_MAX;
	if(P > w){
		R[I].x = R_new[I].x;
	    R[I].y = R_new[I].y;
	    R[I].z = R_new[I].z;
		*p_acc += 1.; 
	}
}

//***************** COMPUTATION OF ENERGY **********************
void d_det(double n, int I, pos *R, mat_A *A, determinant *det){
   double n2 = n*n;
   double x = R[I].x, y = R[I].y, z = R[I].z;
   double r = sqrt(x*x + y*y + z*z);	
   D[0].dA_dx = x, D[1].dA_dx = x*x-1., D[2].dA_dx = x*y, D[3].dA_dx = x*z;
   D[0].dA_dy = y, D[1].dA_dy = x*y, D[2].dA_dy = y*y-1., D[3].dA_dy = y*z;
   D[0].dA_dz = z, D[1].dA_dz = x*z, D[2].dA_dz = z*y, D[3].dA_dz = z*z-1.;
   for(int l=0; l<N; l++){
   	  D[l].dA_dx *= -exp(-r*r/(2.*n2))/n2;
   	  D[l].dA_dy *= -exp(-r*r/(2.*n2))/n2;
   	  D[l].dA_dz *= -exp(-r*r/(2.*n2))/n2;
   }
   for(int j=0; j<N; j++){
   	det[I].nabla_x += D[j].dA_dx*A[j].inv[I];
   	det[I].nabla_y += D[j].dA_dy*A[j].inv[I];
   	det[I].nabla_z += D[j].dA_dz*A[j].inv[I];
   }
}

void d2_det(double n, int I, pos *R, mat_A *A, determinant *det){
	double n2 = n*n;
	double x = R[I].x, y = R[I].y, z = R[I].z;
    double r = sqrt(x*x + y*y + z*z);	
	D[0].d2A_dx2 = x*x-n2, D[1].d2A_dx2 = x*x*x-x-2.*x*n2, D[2].d2A_dx2 = x*x*y-y*n2, D[3].d2A_dx2 = x*x*z-z*n2;
   D[0].d2A_dy2 = y*y-n2, D[1].d2A_dy2 = x*y*y-x*n2, D[2].d2A_dy2 = y*y*y-y-2.*y*n2, D[3].d2A_dy2 = y*y*z-z*n2;
   D[0].d2A_dz2 = z*z-n2, D[1].d2A_dz2 = x*z*z-x*n2, D[2].d2A_dz2 = y*z*z-y*n2, D[3].d2A_dz2 = z*z*z-z-2.*z*n2;
   for(int l=0; l<N; l++){
   	  D[l].d2A_dx2 *= exp(-r*r/(2.*n2))/(n2*n2);
   	  D[l].d2A_dy2 *= exp(-r*r/(2.*n2))/(n2*n2);
   	  D[l].d2A_dz2 *= exp(-r*r/(2.*n2))/(n2*n2);
   }
    for(int j=0; j<N; j++){
   	det[I].lap_x += D[j].d2A_dx2*A[j].inv[I];
   	det[I].lap_y += D[j].d2A_dy2*A[j].inv[I];
   	det[I].lap_z += D[j].d2A_dz2*A[j].inv[I];
   }
}

double energy(pos *R_a, pos *R_b, determinant *det_a, determinant *det_b){
	double E=0.;
	for(int I=0; I<N; I++){
		E -= 0.5*(det_a[I].lap_x + det_a[I].lap_y + det_a[I].lap_z);
		E -= 0.5*(det_b[I].lap_x + det_b[I].lap_y + det_b[I].lap_z);
		E += 0.5*(R_a[I].x*R_a[I].x + R_a[I].y*R_a[I].y + R_a[I].z*R_a[I].z);
		E += 0.5*(R_b[I].x*R_b[I].x + R_b[I].y*R_b[I].y + R_b[I].z*R_b[I].z);
	}
	for(int I=0; I<N; I++){
		E -= det_a[I].nabla_x*det_b[I].nabla_x + det_a[I].nabla_y*det_b[I].nabla_y + det_a[I].nabla_z*det_b[I].nabla_z;
	}
	return E;
}
