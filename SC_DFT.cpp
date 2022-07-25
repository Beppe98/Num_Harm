//Ci sono fluttuazioni. Penso che siano dovute più che altro all'orbitale 2s. Sembra quasi che a volte mi prenda l'autovalore dell'orbitale 1s. 
//Per completezza potresti anche farti stampare le energie (vere) finchè ce le hai 
//Potresti anche vedere il plot degli autovalori in funzione degli step
//Prova a guardare l'orbitale 2s
//Conviene farsi stampare gli orbitali solo alla fine per velocizzare
//La soluzione stava tutta in F(E)


//Prova ad aumentare r_max (vai fino a 30 aumentando di conseguenza il numero di punti). Se aumenti r_max il Numerov diventa più preciso perché la funzione va a zero bene. Prova poi a giocare anche con il numero di punti nella mesh in modo da avere una stima degli integrali migliore!
//set logscale y per il plot della densità 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

//PARAMETERS OF THE SIMULATION
double rs = 3.93;
int Ne = 40;
double Rc = pow(Ne, 1./3.)*rs;
double r_max = 30;
double rhob = 3./(4.*M_PI*pow(rs, 3.));

//NUMEROV MESH
const int N = 6000; 
double h=r_max/N;

//POINTS IN THE NUMEROV ENERGY MESH
const int M = 1500; 

//NUMBER OF SELF-CONSISTENT ITERATIONS
const int N_iter = 100;

//MESH OF POINTS
double r[N] = {};

//ARRAYS THAT STORE THE ENERGIES COMPUTED WITH THE TWO EQUIVALENT FORMULAS
double eps[6]; 
double E_kin[6]; 

//DENSITY ARRAYS
double p[N] = {};
double p_new[N] = {}; 

//POTENTIAL ARRAYS
double v_eff[N] = {};
double v_H[N] = {};
double v_x[N] = {};
double v_c[N] = {};
double v_KS[N] = {};


//INITIAL DENSITY GUESS
void p_guess();

//PRINCIPAL CYCLE WITH NUMEROV
void start(double y[], int l);
double numL(double *y, double E, int l, int *count);
double numR(double *y, double E, int l, int *count);
void print_Psi(double *y, int n, int l);
double norm(double *y);
void num_cycle(int n, int l);
double k(int i, double E);
double find_dE();

//FUNCTIONS FOR FILLING THE POTENTIAL ARRAYS
void fill_v_Hxc();
void fill_v_eff(int l);

//DENSITY IS UPDATED ACCORDING TO THE SELF-CONSISTENT PROCEDURE
double norm_p();
void update_p();
void self_con();

//ENERGY CHECK
double kin_energy(double y[], int l);	
double energy_1(double sum_E_kin);
double energy_2(double sum_eps);
void print_energies();

//POL. CALCULATION
double polarizability();


//******************** CORE OF THE ALGORITHM *****************************
int main(){
   int i;
   double sum_eps=0., sum_E_kin=0., E1=0., E2=0.;
   
   //INITIALIZE THE MESH OF POINTS
	r[0] = h/2;
	for(i=0; i<N-1; i++){
		r[i+1] = r[i] + h;
	}
   //STARTING ELECTRON DENSITY (FERMI DISTRIBUTION)
	p_guess();
	
   //SELF CONSISTENT STEP & EIGENFUNCTIONS PLOT
	self_con();
	
   //DENSITY PLOT 
	ofstream os;
    os.open("p.txt", ios :: out | ios :: trunc);
    os.precision(6);
     for(int j=0; j<N; j++){
        os << r[j] << "       " << p[j] << endl; 
     }
    os.close();
    
    //POLARIZABILITY CALCULATION  
    double pol=0.; 
    pol = polarizability();
    cout << "Polarizability" << "       " << pol << endl; 
}

//******************** PRINT THE EIGENVALUES AT EACH STEP *********************
void self_con(){
	int i;
	ofstream myfile;
	myfile.open("eigenv.txt", ios :: out | ios :: trunc);
    myfile.precision(6);
	for(i=0; i<N_iter; i++){	
	   fill_v_Hxc();
	   num_cycle(0, 0); //1s
	   num_cycle(0, 1); //1p
	   num_cycle(0, 2); //1d
	   num_cycle(0, 3); //1f
	   num_cycle(1, 1); //2p
	   num_cycle(1, 0); //2s
	   update_p();
	   myfile << eps[0] << "     "<< eps[1]<< "       "<< eps[2] << "       " << eps[3] << "         " << eps[4]<< "       " << eps[5] << "         " <<  norm_p() << endl;	 
	   print_energies();
	}
	myfile.close();
}


//******************** CYCLE ON THE ENERGY MESH *********************
//******************** CALCULATION OF THE NEW DENSITY ***************
void num_cycle(int n, int l){
   double dE=0., F=0., F_old=0.;	
	fill_v_eff(l);
	dE = find_dE();
	int count = 0, tmp = 0;
	//NUMEROV ALGORITHM
	int i;
	double y[N] = {};
	for(i=M-1; i>=0; i--){
	    start(y, l);
		F_old = F;
	   F = numL(y, i*dE, l, &count);
	   F -= numR(y, i*dE, l, &count);
		if(F*F_old < 0 && abs(F_old) < 0.2 && abs(F) < 0.2){
		  print_Psi(y, n, l);
		 eps[4*n + l] = 2.*(2.*l + 1)*i*dE;
		  if(tmp == n){
			  break;	
		   }
		  tmp = tmp + 1;		   
		}
	} 
	//KINETIC ENERGIES 
	E_kin[4*n + l] = kin_energy(y, l);
	for(int k=0; k<N; k++){
		p_new[k] += y[k]*y[k]*(2*l + 1)/(2.*M_PI);
   }
}


double k(int i, double E){
	return 2.*(E - v_KS[i]);       
}


double find_dE(){
	double E_min = 0.;
	for(int i=0; i<N; i++){
		if(v_KS[i] < E_min){
			E_min = v_KS[i];
		} 
	}
	return E_min/M;
}


//***************** COMPUTE AND PRINT THE POTENTIALS *********************
void fill_v_eff(int l){
	ofstream file;
	file.open("v_Hxc.txt", ios :: out | ios :: trunc);
	for(int i=0; i<N; i++){
	   if(r[i] <= Rc){
		  v_eff[i] = 2.*M_PI*rhob*(r[i]*r[i]/3. - Rc*Rc) + l*(l + 1)/(2.*r[i]*r[i]);
	    }
	   else{
		  v_eff[i] = -4.*M_PI*rhob*Rc*Rc*Rc/(3.*r[i]) + l*(l + 1)/(2.*r[i]*r[i]);
	    }
	    v_KS[i] = v_eff[i] + v_H[i] + v_x[i] + v_c[i];
		file << r[i] << "      " << v_eff[i] << "      " << v_H[i] << "      " << v_x[i] << "      " << v_c[i]  << "      " << v_KS[i] << endl; 
	}
	file.close();
}


void fill_v_Hxc(){
	int i;
	double I1=0., I2=0., dr=r_max/N;
   double f[N] = {}, g[N] = {}, h[N] = {}, l[N] = {}, r_s[N] = {};
	double A = 0.031091, a1 = 0.21370;
	double b1 = 7.5957, b2 = 3.5876, b3 = 1.6382, b4 = 0.49294;
	for(i=0; i<N; i++){
		I2 += p[i]*r[i]*dr; //for v_H
		v_x[i] = -pow(3.*p[i]/M_PI, 1./3.); //for v_x
		r_s[i] = pow(3./(4.*M_PI*p[i]), 1./3.); //for v_c
		f[i] = -2.*A*(1. + a1*r_s[i]);
		h[i] = 2.*A*(b1*pow(r_s[i], 0.5) + b2*r_s[i] + b3*pow(r_s[i], 1.5) + b4*pow(r_s[i], 2.));
		g[i] = 1. + 1./h[i];
	   l[i] = (2.*A/3.)*r_s[i]*(0.5*b1*pow(r_s[i], -0.5) + b2 + 1.5*b3*pow(r_s[i], 0.5) + 2.*b4*r_s[i]);
	   v_c[i] = f[i]*log(g[i]) + (2.*A/3.)*a1*r_s[i]*log(g[i]) + f[i]*l[i]/(h[i]*(h[i] + 1.));
	}
	for(i=0; i<N; i++){
		I1 += p[i]*r[i]*r[i]*dr;
		I2 -= p[i]*r[i]*dr;
	   v_H[i] = 4*M_PI*(I1/r[i] + I2);
	}		
}


double polarizability(){
	int i, R;
	double dr = r_max/N, I=0., alpha=0.;
	R = round(Rc/dr);
	for(i=R; i<N; i++){
		I += 4.*M_PI*p[i]*r[i]*r[i]*dr;
	}
    alpha = (1. + I/Ne);
	return alpha;
}


//**************** SECTION FOR THE TWO ENERGIES ***************
double energy_1(double sum_E_kin){
	int i;
	double I=0.;
	double dr = r_max/N;
	for(i=0; i<N; i++){
		I += 2.*M_PI*p[i]*v_H[i]*r[i]*r[i]*dr;
		I += 4.*M_PI*p[i]*v_eff[i]*r[i]*r[i]*dr;
		I -= 3.*M_PI*pow(3./M_PI, 1./3.)*pow(p[i], 4./3.)*r[i]*r[i]*dr;
	}
	I = I + sum_E_kin;
	return I;
}


double energy_2(double sum_eps){
	int i;
	double I=0.;
	double dr = r_max/N;
	for(i=0; i<N; i++){
		I -= 2.*M_PI*p[i]*v_H[i]*r[i]*r[i]*dr;
	   I -= 4.*M_PI*p[i]*(v_x[i] + v_c[i])*r[i]*r[i]*dr;
		I -= 3.*M_PI*pow(3./M_PI, 1./3.)*pow(p[i], 4./3.)*r[i]*r[i]*dr;
	}
	I = I + sum_eps;
	return I;
}


void print_energies(){
	double sum_eps = 0.;
    double sum_E_kin = 0.;
    for(int i=0; i<6; i++){
    	sum_eps += eps[i];
       sum_E_kin += E_kin[i];
    }
    double E1 = energy_1(sum_E_kin);
    double E2 = energy_2(sum_eps);
    cout << "E1: " << E1 << "      " << "E2: " << E2 << endl;
}


double kin_energy(double y[], int l){
	int i;
	double I=0.;
	double der[N] = {};
	double der2[N] = {};
	double dr = r_max/N;
	for(i=0; i<N-1; i++){
		der[i] = (y[i+1] - y[i])/dr;
	}
	der[N-1] = der[N-2];
	for(i=0; i<N; i++){
	   I += r[i]*r[i]*der[i]*der[i]*dr;
	}
	return I*2.*(2.*l + 1.);
}


void p_guess(){
	double p0 = 3.*Ne/(4.*M_PI*Rc*Rc*Rc);
	double mu = 6.;
	for(int i=0; i<N; i++){
		p[i] = p0/(1. + exp(mu*(r[i] - Rc)));
	}
}


void update_p(){
    double alpha = 0.1;
	for(int i=0; i<N; i++){
		p[i] = p_new[i]*alpha + p[i]*(1. - alpha);
	   p_new[i] = 0;
	}
}


void start(double y[], int l){ 
    y[0] = pow(r[0], l+1);
    y[1] = pow(r[1], l+1);
    y[N-1] = 0;
    y[N-2] = pow(10, -10);
}

 
double numL(double *y, double E, int l, int *count){
	double A = 0, B = 0, C = 0;
	for(int i = 1; i < N - 1; i++){
		A = h*h*5*k(i, E)/12.;
		B = h*h*k(i-1, E)/12.;
		C = h*h*k(i+1, E)/12.;
		y[i+1] = 2.*(1 - A)*y[i] - (1 + B)*y[i-1];
		y[i+1] = y[i+1]/(1 + C);
		if(k(i, E)*k(i+1, E)<0){
			*count = i;
			break;
		}
	}
	return (y[*count+1] - y[*count-1])/(2.*h*y[*count]);
}


double numR(double *y, double E, int l, int *count){
	double A = 0, B = 0, C = 0;
	for(int i = N-2; i >= 1; i--){
		A = h*h*5*k(i, E)/12.;
		B = h*h*k(i-1, E)/12.;
		C = h*h*k(i+1, E)/12.;
		y[i-1] = 2.*(1 - A)*y[i] - (1 + C)*y[i+1];
		y[i-1] = y[i-1]/(1 + B);
		if(i == *count){
			break;
		}
	}
	return (y[*count+1] - y[*count-1])/(2.*h*y[*count]);
}


void print_Psi(double *y, int n, int l){
	int j, k;
	double sum=0., node=0.;
    double c = 1.0;
    for(j=0; j<N-1; j++){
    	if(y[j]*y[j+1]<0){
    		node = r[j];
    		break;
    	}
    }
    for(j=1; j<N; j++){
    	 if(j != 0 && abs(y[j-1]/y[j])>1E3 && r[j]>node){
           k = j;
           break;
        }
    }
    c = y[k-1]/y[k];
    for(j=k; j<N; j++){
        y[j] = c*y[j];
     }
    sum = norm(y); 
    ofstream os;
    string Psi;
    Psi = "n=" + to_string(n+1) + "_" + "l=" + to_string(l) + ".txt"; 
    os.open(Psi, ios :: out | ios :: trunc);
    os.precision(6);
     for(j=0; j<N; j++){
        y[j] = y[j]/(r[j]*sum);
        os << r[j] << "       " << y[j] << endl; 
     }
    os.close();
}


double norm(double *y){
	double sum=0.;
	for(int i=0; i<N; i++){
		sum += y[i]*y[i]*h; 
	}
	return sqrt(sum);
}


double norm_p(){
	double norm=0., dr=r_max/N;
	for(int i=0; i<N; i++){
		norm += 4.*M_PI*p[i]*r[i]*r[i]*dr;
	}
	return norm;
}

