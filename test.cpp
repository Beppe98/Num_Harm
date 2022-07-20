//Non posso far partire la simulazione con tutti i fermioni nelle stesse posizioni!

#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

using namespace std;

const int N = 4;
const int N_part = 8;
const int M = 1000;
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

struct determinant{
   double nabla_x;
   double nabla_y;
   double nabla_z;
   double lap_x;
   double lap_y;
   double lap_z; 
};

struct mat_A{
	double s;
	double px;
	double py;
	double pz;
   double inv[N];
};

struct mat_B{
	double x[N];
	double y[N];
    double z[N];
    double x2[N];
    double y2[N];
    double z2[N];
};

struct der{
   double dA_dx;
   double dA_dy;
   double dA_dz;
   double d2A_dx2;
   double d2A_dy2;
   double d2A_dz2;
};

pos *R;
mat_A *A;
der *D;
determinant *det;
mat_B *B;


void start();
double mc_step(double x);
void mc_update();
void fill_A(double n, int set);
void fill_A_new(double n, int set);
double compute_det();
void acc_rej(double P, double *p_acc);
double p(double det, double det_new);
void compute_inverse();
void fill_nabla_d2_A(int J, int n);
void compute_B(int J, int n);
void laplacian_det(int J);
void nabla_det(int J);
void get_B2();
double local_energy();


int main(){
	double det_old, det_new, P, p_acc, n, E, E2, E_tot; 
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
	//CYCLE ON THE PARAMETER n
	for(int l=0; l<50; l++){
		n = 0.5 + l*0.05;
		E=0., E2=0., E_tot=0.;
		det_old=0., det_new=0., P=0., p_acc=0.;
		//LONG CYCLE ON THE MC STEPS
	for(int i=0; i<M; i++){
		mc_update();
		
		//SPIN-UP DETERMINANT
	    fill_A(n, 0);
	    det_old = compute_det();
	    fill_A_new(n, 0);
	    det_new = compute_det();
	    P = p(det_old, det_new);
	    
	    //SPIN-DOWN DETERMINANT
	    fill_A(n, 4);
	    det_old = compute_det();
	    fill_A_new(n, 4);
	    det_new = compute_det();
	    
	    //ACCEPTANCE/REJECTION TEST	    
	    P = P*p(det_old, det_new);
	    acc_rej(P, &p_acc);
	    
	    //ENERGY CALCULATION
	    for(int set=0; set<2; set++){
	    	fill_A(n, set*4);
	        compute_inverse();
	        for(int J=0; J<N_part/2; J++){
	        	D = new der[N];
	        	det = new determinant[N_part];
	            B = new mat_B[N];
	    	    fill_nabla_d2_A(J + set*4, n);
	    	    compute_B(J + set*4, n);
	    	    get_B2();
	    	    nabla_det(J + set*4);
	    	    laplacian_det(J + set*4);
	        }
	    }
	    E = local_energy()/M;
	    E_tot += E;
	    E2 += E*E*M;
	}     //END OF THE MC STEPS CYCLE   
	cout << E_tot << " +/- " << sqrt((E2 - E_tot*E_tot)/M) << "        " <<  p_acc/M  << endl;
	}
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


//**************** MATRIX FILLING ********************
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
//If set == 4 then J must be 4..7
void fill_nabla_d2_A(int J, int n){
   double n2 = n*n;
   double x = R[J].x, y = R[J].y, z = R[J].z;
   double r = sqrt(x*x + y*y + z*z);	
   D[0].dA_dx = x, D[1].dA_dx = x*x-1., D[2].dA_dx = x*y, D[3].dA_dx = x*z;
   D[0].dA_dy = y, D[1].dA_dy = x*y, D[2].dA_dy = y*y-1., D[3].dA_dy = y*z;
   D[0].dA_dz = z, D[1].dA_dz = x*z, D[2].dA_dz = z*y, D[3].dA_dz = z*z-1.;
   for(int l=0; l<N; l++){
   	  D[l].dA_dx *= -exp(-r*r/(2.*n2))/n2;
   	  D[l].dA_dy *= -exp(-r*r/(2.*n2))/n2;
   	  D[l].dA_dz *= -exp(-r*r/(2.*n2))/n2;
   }
   D[0].d2A_dx2 = x*x-n2, D[1].d2A_dx2 = x*x*x-x-2.*x*n2, D[2].d2A_dx2 = x*x*y-y*n2, D[3].d2A_dx2 = x*x*z-z*n2;
   D[0].d2A_dy2 = y*y-n2, D[1].d2A_dy2 = x*y*y-x*n2, D[2].d2A_dy2 = y*y*y-y-2.*y*n2, D[3].d2A_dy2 = y*y*z-z*n2;
   D[0].d2A_dz2 = z*z-n2, D[1].d2A_dz2 = x*z*z-x*n2, D[2].d2A_dz2 = y*z*z-y*n2, D[3].d2A_dz2 = z*z*z-z-2.*z*n2;
   for(int l=0; l<N; l++){
   	  D[l].d2A_dx2 *= exp(-r*r/(2.*n2))/(n2*n2);
   	  D[l].d2A_dy2 *= exp(-r*r/(2.*n2))/(n2*n2);
   	  D[l].d2A_dz2 *= exp(-r*r/(2.*n2))/(n2*n2);
   }
}


void compute_B(int J, int n){
	double r = sqrt(R[J].x*R[J].x + R[J].y*R[J].y + R[J].z*R[J].z);	
	for(int l=0; l<N; l++){
		for(int k=0; k<N; k++){
			B[l].x[k] = A[l].inv[J]*D[k].dA_dx;
			B[l].y[k] = A[l].inv[J]*D[k].dA_dy;
			B[l].z[k] = A[l].inv[J]*D[k].dA_dz;
		}
	}
}


void nabla_det(int J){
	for(int k=0; k<N; k++){
		det[J].nabla_x += B[k].x[k];
		det[J].nabla_y += B[k].y[k];
		det[J].nabla_z += B[k].z[k];
	}
	
}


void get_B2(){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			for(int k=0; k<N; k++){
				B[i].x2[j] += B[i].x[k]*B[k].x[j];
				B[i].y2[j] += B[i].y[k]*B[k].y[j];
				B[i].z2[j] += B[i].z[k]*B[k].z[j];
			}
		}
	}
}


void laplacian_det(int J){
	for(int k=0; k<N; k++){
		det[J].lap_x -= B[k].x2[k];
		det[J].lap_y -= B[k].y2[k];
		det[J].lap_z -= B[k].z2[k];
		det[J].lap_x += A[k].inv[J]*D[k].d2A_dx2;
		det[J].lap_y += A[k].inv[J]*D[k].d2A_dy2;
		det[J].lap_z += A[k].inv[J]*D[k].d2A_dz2;
	}
	det[J].lap_x += det[J].nabla_x*det[J].nabla_x;
	det[J].lap_y += det[J].nabla_y*det[J].nabla_y;
	det[J].lap_z += det[J].nabla_z*det[J].nabla_z;
}


double local_energy(){
	double E=0.;
	for(int J=0; J<N_part/2; J++){
		E -= 0.5*(det[J].lap_x + det[J].lap_y + det[J].lap_z);
		E -= 0.5*(det[J + 4].lap_x + det[J + 4].lap_y + det[J + 4].lap_z);
		E -= 0.5*(det[J].nabla_x*det[J].nabla_x + det[J].nabla_y*det[J].nabla_y + det[J].nabla_z*det[J].nabla_z);
		E -= (det[J].nabla_x*det[J + 4].nabla_x + det[J].nabla_y*det[J + 4].nabla_y + det[J].nabla_z*det[J + 4].nabla_z);
		E -= 0.5*(det[J + 4].nabla_x*det[J + 4].nabla_x + det[J + 4].nabla_y*det[J + 4].nabla_y + det[J + 4].nabla_z*det[J + 4].nabla_z);
		E += 0.5*(R[J].x*R[J].x + R[J].y*R[J].y + R[J].z*R[J].z);
		E += 0.5*(R[J + 4].x*R[J + 4].x + R[J + 4].y*R[J + 4].y + R[J + 4].z*R[J + 4].z);
	} 
	return E;
}

//al posto di 4 potevo mettere anche N_part/2







