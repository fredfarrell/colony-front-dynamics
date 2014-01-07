#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;

#define SIZE 500
#define SGN(x) (double)((x > 0) - (x < 0))
#define SQR(x) x*x


vector<double> rx,ry,deltarx,deltary;
double dt=0.01;
double ds=1.0;
double v=-1.0;
double D=1.0;
double ymin,ymax;



void init() {

	
	for(int i=0;i<SIZE;i++) {

		//rx.push_back(i);
		//ry.push_back(100*exp(-(i-500)*(i-500)/1000.0));
		//ry.push_back(0);

		//CIRCLE
		//rx.push_back( cos(i*i*2*M_PI/((double)SQR(SIZE))) );
		//ry.push_back( sin(i*i*2*M_PI/((double)SQR(SIZE))) );

		rx.push_back( cos(i*2*M_PI/((double)SIZE)) );
		ry.push_back( sin(i*2*M_PI/((double)SIZE)) );

		deltarx.push_back(0);
		deltary.push_back(0);		

		ymin=0;
		ymax=100;

	}

}

void print(int j) {

	char str[30];

	sprintf(str,"out%i.dat",j);
	ofstream outp(str);

	for(int i=0;i<rx.size();i++) {
	
		outp<<i<<" "<<rx[i]<<" "<<ry[i]<<endl;

	}

	outp.close();

}

//functions to calculate i+1 and i-1 with PBCs
int up(int x) {

	if(x==rx.size()-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return rx.size()-1;
	else return x-1;

}


void timestep(int tm) {


	for(int i=0;i<rx.size();i++) {


		double drxds = (2*rx[up(i)] + 3*rx[i] - 6*rx[dwn(i)] + rx[dwn(dwn(i))])/(6.0*ds); //3rd order upwind scheme w/ different sign for the two dimensions
		double dryds = (-ry[up(up(i))] + 6*ry[up(i)] - 3*ry[i] - 2*ry[dwn(i)])/(6.0*ds); 

		double d2rxds2 = (rx[up(i)]+rx[dwn(i)]-2*rx[i])/(ds*ds);
		double d2ryds2 = (ry[up(i)]+ry[dwn(i)]-2*ry[i])/(ds*ds);

		double norm = sqrt(drxds*drxds+dryds*dryds);

		drxds/=norm; dryds/=norm;

		deltarx[i] =  - dt*v*dryds + D*dt*d2rxds2;
		deltary[i] =    dt*v*drxds + D*dt*d2ryds2;

	}


	for(int i=0;i<rx.size();i++) { //move the points

		rx[i]+=deltarx[i]; 
		ry[i]+=deltary[i];
	}


}



int main() {

	init();
	print(1000);

	for (int i=0;i<500;i++) {
		
		timestep(i);
		//if(i%100==0) print(i);
		print(i);
	}

	print(2000);


	return 0;

}
