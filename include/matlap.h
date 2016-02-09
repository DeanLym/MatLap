/*
 * matlap.h
 *
 *  Created on: Jan 26, 2016
 *      Author: yiminliu, zhaoyang jin
 */

#ifndef MATLAP__H_
#define MATLAP__H_

#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <ctime>
#include <fstream>
#include <sstream>
using namespace::std;

namespace MatLap
{
class Matrix
{
public:
	Matrix();
	Matrix(int nRow,int nCol);
	~Matrix();


	Matrix& operator = ( const Matrix& other);
	Matrix operator+ (Matrix& other);
	Matrix operator+ (double other);
	Matrix& operator+= (double other);
	Matrix& operator- ();
	Matrix operator- (Matrix& other);
	Matrix operator- (double other);
	Matrix& operator-= (double other);
	Matrix operator*(Matrix& other);
	friend Matrix operator+(double, Matrix &) ;
	friend Matrix operator-(double, Matrix &) ;
	friend ostream& operator<< (ostream&,Matrix &);
	friend Matrix ones(int nRow,int nCol);
	friend Matrix zeros(int nRow,int nCol);
	friend Matrix rand(int nRow,int nCol);
	friend Matrix size(Matrix& x);
	friend Matrix load(char* fn);
	friend bool svd(Matrix A, Matrix& U, Matrix& Sig);
protected:
	int get_nRow();
	int get_nCol();

private:
	vector<double> data_;
	int nRow_;
	int nCol_;
	bool sysm_;

public:
	static const string ERR_EMPTY_MATRIX;
	static const string ERR_DIM_MISMATCH;
	static const string ERR_SVD_FAIL;
	static const string ERR_NON_POSITIVE_DIM;
	static const string ERR_FILE_NOT_EXIST;
	static const string ERR_INCONSISTENT_DIM_IN_DATAFILE;
};

Matrix ones(int nRow,int nCol);
Matrix zeros(int nRow,int nCol);
Matrix ones(int nRow,int nCol);
Matrix rand(int nRow,int nCol);
Matrix size(Matrix& x);
Matrix load(char* fn);
vector<int> sub2lin(vector<int> dim,vector<int> I, vector<int> J);
int sub2lin(vector<int> dim,int I, int J);
bool svd(Matrix A, Matrix& U, Matrix& Sig);

}
#endif /* MATLAP_H_ */
