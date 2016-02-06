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
using namespace::std;

namespace MatLap
{
class Matrix
{
public:
	Matrix();
	Matrix(int nRow,int nCol);
	~Matrix();
	static Matrix zeros(int nRow,int nCol);
	static Matrix ones(int nRow,int nCol);
	static Matrix size(Matrix x);
	static vector<Matrix> svd(Matrix x);
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
protected:
	int get_nRow();
	int get_nCol();

public:
	vector<double> data_;
	int nRow_;
	int nCol_;
	bool sysm_;

private:
	static const string ERR_EMPTY_MATRIX;
	static const string ERR_DIM_MISMATCH;
};

}
#endif /* MATLAP_H_ */
