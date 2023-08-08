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

#include "MaterialTensorOperations.h"


double TensorM::macCaulay(const double N) { double val = (N + fabs(N)) * 0.5; return val; }


Vector TensorM::I(const int N) {
	// return identity matrix in Voigt notation
	Vector kD(N);
	int Dim = 0;
	kD.Zero();
	if (N == 3) {
		Dim = 2;
	}
	else if (N == 6) {
		Dim = 3;
	}
	else {
		opserr << "FATAL: TensorM::I() - invalid material dimension!!\n";
		exit(-1);
	}
	for (int i = 0; i < Dim; ++i) {
		kD(i) = 1.0;
	}
	return kD;
}

Matrix TensorM::I4(const int N) {
	// return 4th order identity tennsor in Voigt notation
	Matrix I4(N, N);
	if (N == 6) {
		I4(0, 0) = 1;
		I4(1, 1) = 1;
		I4(2, 2) = 1;
		I4(3, 3) = 0.5;
		I4(4, 4) = 0.5;
		I4(5, 5) = 0.5;
	}
	else if (N == 3) {
		I4(0, 0) = 0.5;
		I4(1, 1) = 0.5;
		I4(2, 2) = 0.5;
	}
	else {
		opserr << "FATAL: TensorM::I4() - invalid material dimension!!\n";
		exit(-1);
	}
	return I4;
}

Matrix TensorM::IIvol(const int N) {
	// return 4th order Volumetric Tensor in Voigt notation
	// IIvol = I1 tensor I1
	Matrix IIvol(N, N);
	if (N == 6) {
		IIvol(0, 0) = 1;
		IIvol(0, 1) = 1;
		IIvol(0, 2) = 1;
		IIvol(1, 0) = 1;
		IIvol(1, 1) = 1;
		IIvol(1, 2) = 1;
		IIvol(2, 0) = 1;
		IIvol(2, 1) = 1;
		IIvol(2, 2) = 1;
	}
	else if (N == 3) {
		IIvol(0, 0) = 0.5;
		IIvol(0, 1) = 0.5;
		IIvol(1, 0) = 0.5;
		IIvol(1, 1) = 0.5;
	}
	else {
		opserr << "FATAL: TensorM::IIvol() - invalid material dimension!!\n";
		exit(-1);
	}
	return IIvol;
}

Matrix TensorM::IIdev(const int N) {
	//return 4th order Deviatoric Tensor in Voigt notation
	Matrix IIdev(N, N);
	if (N == 6) {
		IIdev(0, 0) = 2.0 / 3.0;
		IIdev(0, 1) = -1.0 / 3.0;
		IIdev(0, 2) = -1.0 / 3.0;
		IIdev(1, 0) = -1.0 / 3.0;
		IIdev(1, 1) = 2.0 / 3.0;
		IIdev(1, 2) = -1.0 / 3.0;
		IIdev(2, 0) = -1.0 / 3.0;
		IIdev(2, 1) = -1.0 / 3.0;
		IIdev(2, 2) = 2.0 / 3.0;
		IIdev(3, 3) = 0.5;
		IIdev(4, 4) = 0.5;
		IIdev(5, 5) = 0.5;
	}
	else if (N == 3) {
		IIdev(0, 0) = 1.0 / 3.0;
		IIdev(0, 1) = -0.5 / 3.0;
		IIdev(1, 0) = -0.5 / 3.0;
		IIdev(1, 1) = 1.0 / 3.0;
		IIdev(2, 2) = 0.25;
	}
	else {
		opserr << "FATAL: TensorM::IIdev() - invalid material dimension!!\n";
		exit(-1);
	}
	return IIdev;
}

double TensorM::determinant(const Matrix& A) {
	double det = 0.0;

	return det;
}

double TensorM::determinant(const Vector& A) {
	double det = 0.0;

	return det;
}

double TensorM::dotdot(const Vector& A, const Vector& B) {
	// c = A^ij (contravariant) * B^ij (contravariant)
	// - double dot operation between two
	//	 square matrices in Voigt notation
	int NA = A.Size();
	int NB = B.Size();
	if (NA != NB) {
		opserr << "FATAL: TensorM::dotdot() - size mismatch!!\n";
		exit(-1);
	}
	Vector Bc(NB);
	if (NB == 3) {
		Bc(0) = B(0); Bc(1) = B(1); Bc(2) = 2 * B(2);
	}
	else if (NB == 6) {
		Bc(0) = B(0); Bc(1) = B(1); Bc(2) = B(2);
		Bc(3) = 2 * B(3); Bc(4) = 2 * B(4); Bc(5) = 2 * B(5);
	}
	else {
		opserr << "FATAL: TensorM::dotdot() - invalid material dimension!!\n";
		exit(-1);
	}

	double scalar = 0;
	for (int i = 0; i < NB; i++) {
		scalar += A(i) * Bc(i);
	}
	return scalar;
}

Matrix TensorM::inner(const Matrix& A, const Matrix& B) {
	// C^ij^kl (contravariant) = A^ij_mn (mixed-variant) * B_mn^kl (covariant)
	// - inner product between two 4th 
	//	 order tensors in Voigt notation
	int NAr = A.noRows();
	int NAc = A.noCols();
	int NBr = B.noRows();
	int NBc = B.noCols();
	if ((NAr != NAc) || (NBr != NBc)) {
		opserr << "FATAL: TensorM::inner() - matrices A and B must be square!!\n";
		exit(-1);
	}
	if ((NBc != NAc) || (NBr != NAr)) {
		opserr << "FATAL: TensorM::inner() - matrices A & B size mismatch!!\n";
		exit(-1);
	}
	Matrix Ac(NAr, NAc);
	if (NAc == 3) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (j == 2) {
					Ac(i, j) = 2 * A(i, j);
				}
				else {
					Ac(i, j) = A(i, j);
				}
			}
		}
	}
	else if (NAc == 6) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (j > 2) {
					Ac(i, j) = 2 * A(i, j);
				}
				else {
					Ac(i, j) = A(i, j);
				}
			}
		}
	}
	else {
		opserr << "FATAL: TensorM::inner() - invalid material dimension!!\n";
		exit(-1);
	}
	Matrix C(NBc, NBc);
	for (int i = 0; i < NBc; i++) {
		for (int j = 0; j < NBc; j++) {
			C(i, j) += Ac(i, j) * B(i, j);
		}
	}
	return C;
}

Vector TensorM::inner(const Matrix& A, const Vector& B) {
	// C^ij (contravariant) = A^ij_kl (mixed-variant) * B_kl (covariant)  
	// - double dot operation between a 4th order tensor
	//	 and a square matrice in Voigt notation
	int NAr = A.noRows();
	int NAc = A.noCols();
	int NB = B.Size();
	if (NAr != NAc) {
		opserr << "FATAL: TensorM::inner() - matrix A must be square!!\n";
		exit(-1);
	}
	if (NAc != NB) {
		opserr << "FATAL: TensorM::inner() - size mismatch!!\n";
		exit(-1);
	}
	//Matrix Ac(NAr, NAc);
	//if (NAc == 3) {
	//	for (int i = 0; i < 3; i++) {
	//		for (int j = 0; j < 3; j++) {
	//			if (j == 2) {
	//				Ac(i, j) = 2 * A(i, j);
	//			}
	//			else {
	//				Ac(i, j) = A(i, j);
	//			}
	//		}
	//	}
	//}
	//else if (NAc == 6) {
	//	for (int i = 0; i < 6; i++) {
	//		for (int j = 0; j < 6; j++) {
	//			if (j > 2) {
	//				Ac(i, j) = 2 * A(i, j);
	//			}
	//			else {
	//				Ac(i, j) = A(i, j);
	//			}
	//		}
	//	}
	//}
	//else {
	//	opserr << "FATAL: TensorM::inner() - invalid material dimension!!\n";
	//	exit(-1);
	//}
	Vector C(NAc);
	for (int i = 0; i < NAc; i++) {
		for (int j = 0; j < NAc; j++) {
			C(i) += A(i, j) * B(j);
		}
	}
	return C;
}

Vector TensorM::inner(const Vector& A, const Matrix& B) {
	// C_kl (covariant) = A^ij (contravariant) * B^ij_kl (mixed-variant)
	// - double dot operation between a 4th order tensor
	//	 and a square matrice in Voigt notation
	int NBr = B.noRows();
	int NBc = B.noCols();
	int NA = A.Size();
	if (NBr != NBc) {
		opserr << "FATAL: TensorM::inner() - matrix A must be square!!\n";
		exit(-1);
	}
	if (NBc != NA) {
		opserr << "FATAL: TensorM::inner() - size mismatch!!\n";
		exit(-1);
	}
	Matrix Bc(NBr, NBc);
	if (NBc == 3) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (j == 2) {
					Bc(i, j) = 2 * B(i, j);
				}
				else {
					Bc(i, j) = B(i, j);
				}
			}
		}
	}
	else if (NBc == 6) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (j > 2) {
					Bc(i, j) = 2 * B(i, j);
				}
				else {
					Bc(i, j) = B(i, j);
				}
			}
		}
	}
	else {
		opserr << "FATAL: TensorM::inner() - invalid material dimension!!\n";
		exit(-1);
	}
	Vector C(NBc);
	for (int i = 0; i < NBc; i++) {
		for (int j = 0; j < NBc; j++) {
			C(j) += A(i) * Bc(i, j);
		}
	}
	return C;
}

Matrix TensorM::outer(const Vector& A, const Vector& B) {
	// Cijkl = Aij * Bkl - dyad operation between two square matrices 
	//					   return a 4th order tensor in Voigt notation
	int NA = A.Size();
	int NB = B.Size();
	Matrix C(NA, NA);
	for (int i = 0; i < NA; i++) {
		for (int j = 0; j < NA; j++) {
			C(i, j) = A(i) * B(j);
		}
	}
	return C;
}

Matrix TensorM::transpose(const Matrix& A) {
	int NAr = A.noRows();
	int NAc = A.noCols();
	Matrix T(NAc, NAr);
	for (int i = 0; i < NAr; i++) {
		for (int j = 0; j < NAc; j++) {
			T(j, i) = A(i, j);
		}
	}
	return T;
}

void TensorM::sqroot(Vector& A) {
	int NA = A.Size();
	for (int i = 0; i < NA; i++) {
		A(i) = sqrt(A(i));
	}
}

void TensorM::sqroot(Matrix& A) {
	int NAr = A.noRows();
	int NAc = A.noCols();
	for (int i = 0; i < NAr; i++) {
		for (int j = 0; j < NAc; j++) {
			A(i, j) = sqrt(A(i, j));
		}
	}
}