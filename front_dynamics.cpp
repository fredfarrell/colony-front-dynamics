#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;

#define SIZE 1000
#define SGN(x) (double)((x > 0) - (x < 0))


vector<double> rx,ry,deltarx,deltary;
double dt=0.01;
double dx=0.1;
double v=1.0;
double ymin,ymax;



void init() {

	
	for(int i=0;i<SIZE;i++) {

		rx.push_back(i);
		ry.push_back(1000*exp(-(i-500)*(i-500)/1000.0));
		//ry.push_back(0);

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



void timestep(int tm) {


	for(int i=0;i<rx.size();i++) {


		int iup = (i==rx.size()-1)?0:i+1;
		int idwn = (i==0)?rx.size()-1:i-1;

		double drx = rx[iup]-rx[idwn];
		if (fabs(drx)>SIZE*0.5) drx = drx - SIZE*SGN(drx);
		double dry = ry[iup]-ry[idwn];
		double norm = sqrt(drx*drx+dry*dry);

		drx/=norm; dry/=norm;

		double rxup,rxdwn;
		rxup=(rx[iup]-rx[i]<-SIZE*0.5)?rx[iup]+SIZE:rx[iup]; 
		rxdwn=(rx[idwn]-rx[i]>SIZE*0.5)?rx[idwn]-SIZE:rx[idwn]; 

		deltarx[i] = 0.5*(rxup+rxdwn)  - dt*v*dry/(2*dx); //update with Lax scheme
		deltary[i] =  0.5*(ry[iup]+ry[idwn]) + dt*v*drx/(2*dx);


	}


	for(int i=0;i<rx.size();i++) { //move the points

		int iup = (i==rx.size()-1)?0:i+1;
		int idwn = (i==0)?rx.size()-1:i-1;

		double rxup,rxdwn;
		rxup=(rx[iup]-rx[i]<-SIZE*0.5)?rx[iup]+SIZE:rx[iup]; 
		rxdwn=(rx[idwn]-rx[i]>SIZE*0.5)?rx[idwn]-SIZE:rx[idwn]; 


		rx[i]=deltarx[i]; //update with Lax scheme
		ry[i]=deltary[i];
	}


	for(int i=0;i<rx.size();i++) { //PBCs
		if (rx[i]>SIZE) rx[i]-=SIZE;
		if (rx[i]<0) rx[i]+=SIZE;
	}


}



int main() {

	init();
	print(1000);

	for (int i=0;i<50000;i++) {
		
		timestep(i);
		print(i);
		cout<<i<<" "<<rx.size()<<endl;
	}

	print(2000);


	return 0;

}
