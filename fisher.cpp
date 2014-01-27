//2D Fisher equation solver

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;


#define SIZEX 500
#define SIZEY 500


long **n, **m; //(discretized) fields
double **phi,**dphi, **p, **dp; //continuous variables to keep track of small changes to m,n resp.
double dt=0.001;
double sqrtdt;
double dx=0.1; 
double v=1.0;
double D=1.0; //diffusion const
double alpha_1=1.5; //growth rates (alpha_1>alpha_2)
double alpha_2=1.0;
double g=1.0; //noise strength

double pmin=1e-3;
double noisemax=1.0;

//double pmin = (log(g*g*dt)*log(g*dt)*g*g*dt)/9000.0;
//double noisemax = fabs(g*g*dt)/3.0;

int totalshift=0;

// random number generator - like drand48() under Linux, but faster
// works only with compilers which have long long int!
long long int xdr48=0x000100010001LL, mndr48=0x0005deece66dLL, doddr48=0xbLL ;

double _drand48(void)  
{
  xdr48=mndr48*xdr48+doddr48 ; xdr48&=0xffffffffffffLL ;
  return (xdr48/281474976710656.0) ;
}

double f(int x) {

    if(x>100 && x<400) return 5*exp( -(x-SIZEX/2)*(x-SIZEX/2)/500.0 );
	else return 0;

}

void init() {

    cout<<pmin<<" "<<noisemax<<endl;

	phi = new double*[SIZEX]; dphi = new double*[SIZEX]; p = new double*[SIZEX]; dp = new double*[SIZEX]; n = new long*[SIZEX]; m = new long *[SIZEX];

	for(int i=0;i<SIZEX; i++) {
		phi[i]=new double[SIZEY]; dphi[i]=new double[SIZEY]; p[i]=new double[SIZEY]; dp[i]=new double[SIZEY]; n[i] = new long[SIZEY]; m[i] = new long[SIZEY];
	}


	for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {

		//phi[i][j]=tanh( (j-SIZE/2)/10.0 - 5*sin(i*2*M_PI/((double)SIZE)))  ;
        m[i][j]=((0.5*tanh( (SIZEY/4-j)/10.0) + 0.5)/pmin);
        n[i][j]=(((0.5*tanh( (SIZEY/4-j)/10.0) + 0.5) * exp (-(i-SIZEX/2)*(i-SIZEX/2)/100.0))/pmin);
        p[i][j]=0;
        phi[i][j]=0;

	}}

    sqrtdt=sqrt(dt);
}

int up(int x) { //functions for PBCs

	if(x==SIZEX-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return SIZEX-1;
	else return x-1;

}


void timestep() {

	/*for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {
            n[i][j]=p[i][j]/pmin;
            m[i][j]=phi[i][j]/pmin;
        }
    }*/

	for(int i=0;i<SIZEX;i++) {
		for(int j=1;j<SIZEY-1;j++) {


            double d2mdx2 = (m[up(i)][j]+m[dwn(i)][j]-2*m[i][j])/(dx*dx); 
            double d2mdy2 = (m[i][j+1]+m[i][j-1]-2*m[i][j])/(dx*dx);
            double d2ndx2 = (n[up(i)][j]+n[dwn(i)][j]-2*n[i][j])/(dx*dx);  
            double d2ndy2 = (n[i][j+1]+n[i][j-1]-2*n[i][j])/(dx*dx);

			dphi[i][j] = dt*alpha_2*m[i][j]*(1-pmin*m[i][j]) + dt*(alpha_1-alpha_2)*n[i][j]*(1-pmin*m[i][j]) + D*dt*(d2mdx2+d2mdy2);

            cout<<dphi[i][j]<<endl;

            dp[i][j] = dt*alpha_1*n[i][j]*(1-pmin*m[i][j]) + D*dt*(d2ndx2+d2ndy2);// + sqrtdt*g*sqrt(n[i][j]*(m[i][j]-n[i][j])) * noisemax*2*(_drand48()-0.5);
            
            if(isnan(dp[i][j])) cout<<"NAN! "<<i<<" "<<j<<" "<<m[i][j]<<" "<<n[i][j]<<" "<<sqrtdt*g*sqrt(n[i][j]*(m[i][j]-n[i][j])) * noisemax*2*(_drand48()-0.5)<<" "<<n[i][j]*(m[i][j]-n[i][j])<<endl;

            //if( isnan(dp[i][j]) ){ cout<<"TEST: "<<p[i][j]<<" "<<n[i][j]<<" "<<pmin*n[i][j]<<" "<<dt*alpha_1*pmin*n[i][j]*(1-phi[i][j])<<" "<<D*dt*(d2pdx2+d2pdy2)<<" "<<2*sqrtdt*g*sqrt(pmin*n[i][j]*(phi[i][j]-pmin*n[i][j]))*(_drand48()-0.5)<<" "<<pmin*n[i][j]*(phi[i][j]-pmin*n[i][j])<<" "<<p[i][j]*(phi[i][j]-p[i][j])<<endl; }
            //dp[i][j] = dt*alpha_1*p[i][j]*(1-phi[i][j]) + D*dt*(d2pdx2+d2pdy2) + 2*sqrtdt*g*p[i][j]*(phi[i][j]-p[i][j])*(_drand48()-0.5);
           
	}}

	for(int i=0;i<SIZEX;i++) {
		for(int j=1;j<SIZEY-1;j++) {
	
			phi[i][j]+=dphi[i][j];	
            p[i][j]+=dp[i][j];	

            m[i][j] = m[i][j] + (int)phi[i][j];
            n[i][j] = n[i][j] + (int)p[i][j];
            phi[i][j] = phi[i][j] - (int)phi[i][j];
            p[i][j] = p[i][j] - (int)p[i][j];

            //n[i][j]+=dp[i][j];
            //m[i][j]+=dphi[i][j];

            //if(n[i][j]<0) n[i][j]=0;
            if(n[i][j]>m[i][j]) n[i][j]=m[i][j];

	}}


	for(int i=0;i<SIZEX;i++) {
	
            //zero flux BC
			phi[i][0]=phi[i][1];
            phi[i][SIZEY-1]=phi[i][SIZEY-2];		
			p[i][0]=p[i][1];
            p[i][SIZEY-1]=p[i][SIZEY-2];	

	}

}

void printgrid(int t) {

	char str[30];
	sprintf(str,"out%i.dat3",t);
	ofstream outp(str);

	for(int i=1;i<SIZEX-1;i++) {
		for(int j=1;j<SIZEY-1;j++) {

			outp<<i<<" "<<j<<" "<<phi[i][j]<<" "<<p[i][j]<<" "<<m[i][j]<<" "<<n[i][j]<<endl;

	}outp<<endl;}

}

void shiftEverythingDown() { //shift everything to keep the interface roughly centred in the simulation

    double tolerance = 0.01;
    int ymin=SIZEY;
    int buffer=100; //number of lattice sites between the bottom of the contour and the bottom of the sim
    
    for(int i=0;i<SIZEX;i++) {
        for(int j=0;j<SIZEY;j++) {
            
            if(fabs(phi[i][j])<tolerance && j<ymin) ymin=j;
            
        }
    }
    
    int shift = ymin-buffer;
    
    /*if(shift<0) {
	printgrid(3);
        cout<<"Error: hitting bottom of sim box!" <<ymin<<endl;
        exit(0);
    }*/

    if(shift>0) {            
	    for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {
		    
		    int x=j+shift;
		    if(x>=SIZEY) {dphi[i][j]=0; dp[i][j]=0;}
		    else {dphi[i][j]=phi[i][x]; dp[i][j]=p[i][x];}
	    
		}
	    }
     }

    else if(shift<0) {

    	for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {
		    
		    int x=j+shift;
		    if(x<0) {dphi[i][j]=1; dp[i][j]=1;}
		    else {dphi[i][j]=phi[i][x]; dp[i][j]=p[i][x];}
	    
		}
	    }

    }

    else return;

    for(int i=0;i<SIZEX;i++) {
	    for(int j=0;j<SIZEY;j++) {
		    
                phi[i][j]=dphi[i][j]; //using dphi as a temp array to update grid simultaneously
                p[i][j]=dp[i][j];
	    
	    }
	}


    
}

int mutantMaxPosition() {

    int maxpos=0;

    for(int i=0;i<SIZEX;i++) {
        for(int j=0;j<SIZEY;j++) {
            if(p[i][j]>0.01 && i>maxpos) maxpos=i;
        }
    }

    return maxpos;

}

double integratedMutantNumber() {

    double total=0;

    for(int i=0;i<SIZEX;i++) {
        for(int j=0;j<SIZEY;j++) {
            total += p[i][j];       
        }
    }

    return total/((double)(SIZEX*SIZEY));

}

double integratedTotal() {

    double total=0;

    for(int i=0;i<SIZEX;i++) {
        for(int j=0;j<SIZEY;j++) {
            total += phi[i][j];       
        }
    }

    return total/((double)(SIZEX*SIZEY));

}


int main() {

	init();
	printgrid(1);

	for(int i=0;i<50001;i++) {

		cout<<i<<" "<<mutantMaxPosition()<<" "<<integratedMutantNumber()<<" "<<integratedTotal()<<endl;

		timestep();
		//shiftEverythingDown();

        if(i%10000==0) printgrid(i);

        //printgrid(i);

	}


	return 0;

}


