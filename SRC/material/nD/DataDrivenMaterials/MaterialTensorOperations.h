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


class TensorM
{
public:
	TensorM() = default;
	~TensorM() = default;


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