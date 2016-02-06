/*
 * matlap.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: yiminliu
 */
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include "matlap.h"


//using namespace MatLap;
using namespace::std;

#ifdef __cplusplus
extern "C"{
#endif

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

#ifdef __cplusplus
}
#endif





namespace MatLap{

const string Matrix::ERR_DIM_MISMATCH = "Dimension mismatch.";
const string Matrix::ERR_EMPTY_MATRIX = "Empty matrix.";


Matrix::Matrix(){
	nRow_ = 0;
	nCol_ = 0;
	sysm_ = false;
}

Matrix::Matrix(int nRow, int nCol){
	nRow_ = nRow;
	nCol_ = nCol;
	data_.reserve(nRow*nCol);
	sysm_ = false;
}


Matrix::~Matrix(){
	data_.clear();
}


int Matrix::get_nRow(){
	return nRow_;
}


int Matrix::get_nCol(){
	return nCol_;
}



Matrix Matrix::ones(int nRow,int nCol){
	Matrix m(nRow,nCol);
	double one = 1.0;
	vector<double> data(nRow*nCol,one);
	m.data_ = data;
	return m;
}


Matrix& Matrix::operator=(const Matrix& other){
	nRow_ = other.nRow_;
	nCol_ = other.nCol_;
	data_ = other.data_;
	return *this;
}


Matrix Matrix::operator+ (Matrix& other){
	try{
		if( this->nRow_ != other.nRow_ || this->nCol_ != other.nCol_)
			throw runtime_error(ERR_DIM_MISMATCH);
		if( this->data_.empty() || other.data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		Matrix result;
		result = *this;
		integer n = this->data_.size();
		integer icr = 1;
	    	// z = x;
		doublereal alpha = 1.0;
		daxpy_(&n, &alpha, &(other.data_[0]), &icr, &(result.data_[0]), &icr);
			// z = y + z;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;

}

Matrix Matrix::operator+ (double other){
	try{
		if( this->data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		Matrix result;
		result = *this;
	    	// z = x;
		for(int i = 0 ; i < this->data_.size(); i++)
			result.data_[i] += other;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;

}

Matrix& Matrix::operator+= (double other){
	try{
		if( this->data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		for(int i = 0 ; i < this->data_.size(); i++)
			this->data_[i] += other;
		return *this;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;
}

Matrix& Matrix::operator- (){
	try{
		if( this->data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		for(int i = 0 ; i < this->data_.size(); i++)
			this->data_[i] = -this->data_[i];
		return *this;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;

}

Matrix Matrix::operator- (Matrix& other){
	try{
		if( this->nRow_ != other.nRow_ || this->nCol_ != other.nCol_)
			throw runtime_error(ERR_DIM_MISMATCH);
		if( this->data_.empty() || other.data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		Matrix result;
		result = (-other);
		integer n = this->data_.size();
		integer icr = 1;
	    	// z = -y;
		doublereal alpha = 1.0;
		daxpy_(&n, &alpha, &(this->data_[0]), &icr, &(result.data_[0]), &icr);
			// z = x + z = x - y;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;

}

Matrix Matrix::operator- (double other){
	try{
		if( this->data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		Matrix result;
		result = *this;
	    	// z = x;
		for(int i = 0 ; i < this->data_.size(); i++)
			result.data_[i] -= other;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;
}



Matrix& Matrix::operator-= (double other){
	try{
		if( this->data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		for(int i = 0 ; i < this->data_.size(); i++)
			this->data_[i] -= other;
		return *this;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;
}

Matrix Matrix::operator* (Matrix& other){
	try{
		if( this->nCol_ != other.nRow_ )
			throw runtime_error(ERR_DIM_MISMATCH);
		if( this->data_.empty() || other.data_.empty())
			throw runtime_error(ERR_EMPTY_MATRIX);

		integer M = this->nRow_;
		integer K = this->nCol_;
		integer N = other.nCol_;
		Matrix result(M,N);
		result.data_.resize(M*N);
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		dgemm_("N","N",&M,&N,&K,&alpha,&(this->data_[0]),&M,&(other.data_[0]),&K,&beta,&(result.data_[0]),&M);
			// z = y + z;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return *this;

}

Matrix operator+ (double value, Matrix& other){
	try{
		if( other.data_.empty())
			throw runtime_error(Matrix::ERR_EMPTY_MATRIX);

		Matrix result;
		result = other;
	    	// z = x;
		for(int i = 0 ; i < other.data_.size(); i++)
			result.data_[i] += value;
		return result;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		terminate();
	}
	return other;
}
}



