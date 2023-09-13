/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */


#ifndef MaterialTensorOperations_h
#define MaterialTensorOperations_h


#include <math.h> 
#include <Vector.h>
#include <Matrix.h>


/*
// tensor
class Tensor 
{
public:

	//constructors
	Tensor(void) = default;
	Tensor(int nrows);
	Tensor(int nrows, int ncols);
	Tensor(const Tensor&);
	Tensor(double* data, int size);
	Tensor(double* data, int nrows, int ncols);
#ifdef USE_CXX11   
	Tensor(Tensor&&);
#endif
	
	// destructor
	~Tensor(void) = default;





    // overloaded operators 
    inline double& operator()(int row, int col);
    inline double operator()(int row, int col) const;
    Matrix operator()(const ID& rows, const ID& cols) const;

    Matrix& operator=(const Matrix& M);

#ifdef USE_CXX11
    Matrix& operator=(Matrix&& M);
#endif

    // matrix operations which will preserve the derived type and
    // which can be implemented efficiently without many constructor calls.

    // matrix-scalar operations
    Matrix& operator+=(double fact);
    Matrix& operator-=(double fact);
    Matrix& operator*=(double fact);
    Matrix& operator/=(double fact);

    // matrix operations which generate a new Matrix. They are not the
    // most efficient to use, as constructors must be called twice. They
    // however are useful for matlab like expressions involving Matrices.

    // matrix-scalar operations
    Matrix operator+(double fact) const;
    Matrix operator-(double fact) const;
    Matrix operator*(double fact) const;
    Matrix operator/(double fact) const;

    // matrix-vector operations
    Vector operator*(const Vector& V) const;
    Vector operator^(const Vector& V) const;


    // matrix-matrix operations
    Matrix operator+(const Matrix& M) const;
    Matrix operator-(const Matrix& M) const;
    Matrix operator*(const Matrix& M) const;
    //     Matrix operator/(const Matrix &M) const;    
    Matrix operator^(const Matrix& M) const;
    Matrix& operator+=(const Matrix& M);
    Matrix& operator-=(const Matrix& M);






    // overloaded operators
    inline double operator()(int x) const;
    inline double& operator()(int x);
    double operator[](int x) const;  // these two operator do bounds checks
    double& operator[](int x);
    Vector operator()(const ID& rows) const;
    Vector& operator=(const Vector& V);
#ifdef USE_CXX11   
    Vector& operator=(Vector&& V);
#endif
    Vector& operator+=(double fact);
    Vector& operator-=(double fact);
    Vector& operator*=(double fact);
    Vector& operator/=(double fact);

    Vector operator+(double fact) const;
    Vector operator-(double fact) const;
    Vector operator*(double fact) const;
    Vector operator/(double fact) const;

    Vector& operator+=(const Vector& V);
    Vector& operator-=(const Vector& V);

    Vector operator+(const Vector& V) const;
    Vector operator-(const Vector& V) const;
    double operator^(const Vector& V) const;
    Vector operator/(const Matrix& M) const;

    int operator==(const Vector& V) const;
    int operator==(double) const;
    int operator!=(const Vector& V) const;
    int operator!=(double) const;

    //operator added by Manish @ UB
    Matrix operator%(const Vector& V) const;

private:


};

*/


// tensor operations
class TensorM
{
public:
	TensorM(void) = default;
	~TensorM(void) = default;


	static double macCaulay(const double N);
	static Vector quadratic(const double A, const double B, const double C);


	static Vector I(const int N);
	static Matrix I4(const int N);
	static Matrix IIvol(const int N);
	static Matrix IIdev(const int N);
	

	static double determinant(const Matrix& A);
	static double determinant(const Vector& A);


	static void sqroot(Vector& A);
	static void sqroot(Matrix& A);
	static Matrix transpose(const Matrix& A);


	static double dotdot(const Vector& A, const Vector& B);


	static Matrix inner(const Matrix& A, const Matrix& B);
	static Vector inner(const Matrix& A, const Vector& B);
	static Vector inner(const Vector& A, const Matrix& B);


	static Matrix outer(const Vector& A, const Vector& B);
};
#endif