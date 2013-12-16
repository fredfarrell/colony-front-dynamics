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


		//calculate positions of points 1 or 2 up or down from current cell using minimum image
		//sees if the difference between the current cell and (e.g.) i+1 is > SIZE/2, returns a shifted point if so
		double rxup,rxdwn,rx2up,rx2dwn;
		rxup=(rx[up(i)]-rx[i]<-SIZE*0.5)?rx[up(i)]+SIZE:rx[up(i)]; 
		rxdwn=(rx[dwn(i)]-rx[i]>SIZE*0.5)?rx[dwn(i)]-SIZE:rx[dwn(i)]; 
		//rx2up=(rx[up(up(i))]-rx[i]<-SIZE*0.5)?rx[up(up(i))]+SIZE:rx[up(up(i))];  not actually used
		rx2dwn=(rx[dwn(dwn(i))]-rx[i]>SIZE*0.5)?rx[dwn(dwn(i))]-SIZE:rx[dwn(dwn(i))]; 

		double drx = (2*rxup + 3*rx[i] - 6*rxdwn + 2*rx2dwn)/6.0; //3rd order upwind scheme w/ different sign for the two dimensions
		double dry = (-ry[up(up(i))] + 6*ry[up(i)] - 3*ry[i] - 2*ry[dwn(i)])/6.0; 

		double norm = sqrt(drx*drx+dry*dry);

		drx/=norm; dry/=norm;

		deltarx[i] =  - dt*v*dry/dx; //update with upwind scheme
		deltary[i] =    dt*v*drx/dx;

	}


	for(int i=0;i<rx.size();i++) { //move the points

		rx[i]+=deltarx[i]; //update with upwind scheme
		ry[i]+=deltary[i];
	}


	for(int i=0;i<rx.size();i++) { //PBCs
		if (rx[i]>SIZE) rx[i]-=SIZE;
		if (rx[i]<0) rx[i]+=SIZE;
	}


}



int main() {

	init();
	print(1000);

	for (int i=0;i<5000;i++) {
		
		timestep(i);
		print(i);
		cout<<i<<" "<<rx.size()<<endl;
	}

	print(2000);


	return 0;

}
