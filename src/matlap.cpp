/*
 * matlap.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: yiminliu
 */

#include "matlap.h"


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
const string Matrix::ERR_FILE_NOT_EXIST = "File does not exist.";
const string Matrix::ERR_INCONSISTENT_DIM_IN_DATAFILE = "Inconsistent dimension in data file";

Matrix::Matrix(){
	nRow_ = 0;
	nCol_ = 0;
	sysm_ = false;
}

Matrix::Matrix(int nRow, int nCol){
	if(nRow<=0||nCol<=0)
		throw runtime_error(ERR_NON_POSITIVE_DIM);
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


Matrix& Matrix::operator=(const Matrix& other){
	nRow_ = other.nRow_;
	nCol_ = other.nCol_;
	data_ = other.data_;
	return *this;
}


Matrix Matrix::operator+ (Matrix& other){
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

Matrix Matrix::operator+ (double other){
	if( this->data_.empty())
		throw runtime_error(ERR_EMPTY_MATRIX);

	Matrix result;
	result = *this;
	// z = x;
	for(int i = 0 ; i < this->data_.size(); i++)
		result.data_[i] += other;
	return result;

}

Matrix& Matrix::operator+= (double other){
	if( this->data_.empty())
		throw runtime_error(ERR_EMPTY_MATRIX);

	for(int i = 0 ; i < this->data_.size(); i++)
		this->data_[i] += other;
	return *this;
}

Matrix& Matrix::operator- (){
	if( this->data_.empty())
		throw runtime_error(ERR_EMPTY_MATRIX);

	for(int i = 0 ; i < this->data_.size(); i++)
		this->data_[i] = -this->data_[i];
	return *this;

}

Matrix Matrix::operator- (Matrix& other){
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

Matrix Matrix::operator- (double other){
	if( this->data_.empty())
		throw runtime_error(ERR_EMPTY_MATRIX);

	Matrix result;
	result = *this;
	// z = x;
	for(int i = 0 ; i < this->data_.size(); i++)
		result.data_[i] -= other;
	return result;
}



Matrix& Matrix::operator-= (double other){
	if( this->data_.empty())
		throw runtime_error(ERR_EMPTY_MATRIX);

	for(int i = 0 ; i < this->data_.size(); i++)
		this->data_[i] -= other;
	return *this;
}

Matrix Matrix::operator* (Matrix& other){
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

Matrix operator+ (double value, Matrix& other){
	if( other.data_.empty())
		throw runtime_error(Matrix::ERR_EMPTY_MATRIX);

	Matrix result;
	result = other;
	// z = x;
	for(int i = 0 ; i < other.data_.size(); i++)
		result.data_[i] += value;
	return result;
}


Matrix operator- (double value, Matrix& other){
	if( other.data_.empty())
		throw runtime_error(Matrix::ERR_EMPTY_MATRIX);

	Matrix result;
	result = (-other);
	// z = -x;
	for(int i = 0 ; i < other.data_.size(); i++)
		result.data_[i] += value;
	return result;
}

ostream& operator<< (ostream& os , Matrix &X ){
	if( X.data_.empty())
		os <<"[]"<< endl;
	else
		for(int i=0 ; i<X.nRow_; i++){
			for(int j=0; j<X.nCol_;j++)
				os << X.data_[X.nRow_*j+i] << " ";
			os<<endl;
		}
	return os;
}



Matrix zeros(int nRow,int nCol){
	if(nRow<=0||nCol<=0)
		throw runtime_error(Matrix::ERR_NON_POSITIVE_DIM);
	Matrix m(nRow,nCol);
	double zero = 0.0;
	vector<double> data(nRow*nCol,zero);
	m.data_ = data;
	return m;
}

Matrix rand(int nRow,int nCol){
	if(nRow<=0||nCol<=0)
		throw runtime_error(Matrix::ERR_NON_POSITIVE_DIM);
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
	if(nRow<=0||nCol<=0)
		throw runtime_error(Matrix::ERR_NON_POSITIVE_DIM);
	Matrix m(nRow,nCol);
	double one = 1.0;
	vector<double> data(nRow*nCol,one);
	m.data_ = data;
	return m;
}

Matrix load(char* fn){
	Matrix m;
	string line;
	ifstream in(fn);
	int nCol_0 = 0;
	int nCol=0,nRow=0;
	double tmp;
	vector<double> data_tmp;
	if(in.is_open()){
		while(getline(in,line)){
			istringstream str_in(line);
			while(str_in>>tmp){
				data_tmp.push_back(tmp);
				nCol++;
			}
			if((nCol_0 == 0) & (nCol != 0))
				nCol_0 = nCol;
			if((nCol!=0) & (nCol_0!=0) & (nCol==nCol_0)){
				nCol_0 = nCol;
				nCol = 0;
				nRow++;
			}else{
				if((nCol != nCol_0) & (nCol!=0))
					throw runtime_error(Matrix::ERR_INCONSISTENT_DIM_IN_DATAFILE);
			}
		}
		in.close();
		m.nCol_ = nCol_0;
		m.nRow_ = nRow;
		vector<int> dim ={m.nCol_,m.nCol_};
		for(int i=0;i<m.nCol_;i++)
			for(int j=0;j<m.nRow_;j++)
				m.data_.push_back(data_tmp[j*m.nCol_+i]);
	}
	else
		throw runtime_error(Matrix::ERR_FILE_NOT_EXIST);
	return m;
}

bool svd(Matrix A, Matrix& U, Matrix& Sig){
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

int sub2lin(vector<int> dim,int I, int J){
	if(dim.size()!=2)
		throw runtime_error("The first argument must be two dimensional vector.");
	int nRow = dim[0], nCol = dim[1];
	if((I>nRow) | (J>nCol) | (I<0) | (J<0))
		throw runtime_error("Index out of bound.");
	return (J-1)*nRow + I;
}

}



