/*
 * matlap.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: yiminliu
 */

#include "matlap.h"
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <ctime>

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
const string Matrix::ERR_SVD_FAIL = "SVD failed.";
const string Matrix::ERR_NON_POSITIVE_DIM = "Dimension non-positive.";
Matrix::Matrix(){
	nRow_ = 0;
	nCol_ = 0;
	sysm_ = false;
}

Matrix::Matrix(int nRow, int nCol){
	try{
		if(nRow<=0||nCol<=0)
			throw runtime_error(ERR_NON_POSITIVE_DIM);
		nRow_ = nRow;
		nCol_ = nCol;
		data_.reserve(nRow*nCol);
		sysm_ = false;
	}
	catch(runtime_error& e){
		terminate();
	}
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


Matrix operator- (double value, Matrix& other){
	try{
		if( other.data_.empty())
			throw runtime_error(Matrix::ERR_EMPTY_MATRIX);

		Matrix result;
		result = (-other);
		// z = -x;
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

ostream& operator<< (ostream& os , Matrix &X ){
	for(int i=0 ; i<X.nRow_; i++){
		for(int j=0; j<X.nCol_;j++)
			os << X.data_[X.nRow_*j+i] << " ";
		os<<endl;
	}
	return os;
}



Matrix zeros(int nRow,int nCol){
	Matrix m(nRow,nCol);
	double zero = 0.0;
	vector<double> data(nRow*nCol,zero);
	m.data_ = data;
	return m;
}

Matrix rand(int nRow,int nCol){
	Matrix m(nRow,nCol);
	std::srand(std::time(0)); // use current time as seed for random generator
    double random_number;
	for (int i=0; i<nRow*nCol; ++i) {
		random_number = (double)std::rand()/(double)RAND_MAX;
		m.data_.push_back(random_number);
	}
	return m;
}

Matrix ones(int nRow,int nCol){
	Matrix m(nRow,nCol);
	double one = 1.0;
	vector<double> data(nRow*nCol,one);
	m.data_ = data;
	return m;
}


bool svd(Matrix A, Matrix& U, Matrix& Sig){
	try{
		if( A.data_.empty())
			throw runtime_error(Matrix::ERR_EMPTY_MATRIX);

		integer M = A.nRow_, N = A.nCol_;
		integer NCOL;
		NCOL = min(M,N);
		U.nRow_ = M;
		U.nCol_ = NCOL;
		U.data_.resize(M*NCOL);
		Sig.nRow_ = NCOL;
		Sig.nCol_ = 1;
		Sig.data_.resize(NCOL);

		double *VT;
		VT = NULL;
		integer info = -1;
		doublereal work1;
		integer lwork = -1;
		dgesvd_("S","N",&M,&N,&(A.data_[0]),&M,&(Sig.data_[0]),&(U.data_[0]),&M,VT,&N,&work1,&lwork,&info);
		if (info != 0)
			throw runtime_error(Matrix::ERR_SVD_FAIL);
		lwork = (integer)(work1);
		doublereal *work;
		work = new double[lwork];
		dgesvd_("S","N",&M,&N,&(A.data_[0]),&M,&(Sig.data_[0]),&(U.data_[0]),&M,VT,&N,work,&lwork,&info);

		if (info != 0)
			throw runtime_error(Matrix::ERR_SVD_FAIL);
		return true;
	}
	catch (runtime_error& e){
		//cout << e.what() << endl;
		return false;
	}

}


}



