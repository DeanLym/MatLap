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
	double w = 1.5;
	x = Matrix::ones(2,2);
	y = Matrix::ones(2,1);
	x += 1.0;
	z = 5.5 + x;
//	z = x+y;
//	z += w;
//	z = -z;
//	z = z - 1.0;
//	z -= 2.0;
//	z = z-y;
	z = z*y;
	cout << "Number of rows: " << z.nRow_ << endl;
	cout << "Number of cols: " << z.nCol_ << endl;
	cout << "Matrix:" << endl; 
	cout << z.data_[0] << endl <<z.data_[1] << endl;
	return 1;
}
