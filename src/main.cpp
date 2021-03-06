/*
 * main.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: yiminliu
 */

#include "matlap.h"
#include <iostream>
//#include "matlap.cpp"
using namespace::MatLap;
using namespace::std;

int main(){
	Matrix x,y,z;
//	double w = 1.5;
	x = ones(5,2);
	y = ones(2,1);
	x += 1.0;
	z = 5.5 + x;
	z = 10 - z;
//	z = x+y;
//	z += w;
//	z = -z;
//	z = z - 1.0;
//	z -= 2.0;
//	z = z-y;
//	z = x*y;
//	z = rand(3,2);
	z = load("z.txt");
	int tmp;
	tmp = sub2lin({5,3},1,2);
	Matrix U,Sig;
	cout << "z (before svd):" << endl;
	cout << z << endl;
	svd(z,U,Sig);
	cout << "z (after svd):" << endl;
	cout << z << endl << endl;
	cout << "U:" << endl;
	cout << U << endl;
	cout << "Sig:" << endl;
	cout << Sig << endl;
	return 1;
}
