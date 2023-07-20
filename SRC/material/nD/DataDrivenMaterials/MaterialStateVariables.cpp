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

#include "MaterialStateVariables.h"

// Public methods
	// constructor
MaterialStateVariables::MaterialStateVariables(const double nOrd)
{
	if (nOrd == 3) {
		eps = Vector(3);
		eps_commit = Vector(3);
		sig = Vector(3);
		sig_commit = Vector(3);
		sig_implex = Vector(3);
		xs = Vector(3);
		xs_commit = Vector(3);
		Ce = Matrix(3, 3);
		Cep = Matrix(3, 3);
	}
	else if (nOrd == 6) {
		eps = Vector(6);
		eps_commit = Vector(6);
		sig = Vector(6);
		sig_commit = Vector(6);
		sig_implex = Vector(6);
		xs = Vector(6);
		xs_commit = Vector(6);
		Ce = Matrix(6, 6);
		Cep = Matrix(6, 6);
	}
	else {
		opserr << "FATAL: MaterialStateVariables::MaterialStateVariables() -> unknown model dimension!\n";
		exit(-1);;
	}
}

	// destructor
MaterialStateVariables::~MaterialStateVariables(void)
{
	/* 
	   do nothing!
	   no space was allocated in the memory-land to be freed... 
	*/
}

// iteration control
int MaterialStateVariables::commitState(void) {

	// commit internal variables
	eps_commit = eps;
	sig_commit = sig;
	xs_commit = xs;
	lambda_commit_old = lambda_commit;
	lambda_commit = lambda;

	return 0;
}

int MaterialStateVariables::revertToLastCommit(void) {

	// restore committed internal variables
	eps = eps_commit;
	sig = sig_commit;
	xs = xs_commit;
	lambda = lambda_commit;
	lambda_commit = lambda_commit_old;

	return 0;
}

	// pack and unpack state variables (vectorize) for message passing
Vector& MaterialStateVariables::pack(void) {
	int vsz = getSize();
	// initialize vectorized state-variables
	int N = 0;
	bool forward = true;
	Vector data(vsz);
	vectorvector(N, eps, data, forward);
	vectorvector(N, eps_commit, data, forward);

	vectorvector(N, sig, data, forward);
	vectorvector(N, sig_commit, data, forward);

	vectorvector(N, xs, data, forward);
	vectorvector(N, xs_commit, data, forward);

	data(N + 3) = dtime_n;
	data(N + 4) = dtime_n_commit;
	data(N + 5) = dtime_is_user_defined;
	data(N + 6) = dtime_first_set;
	return data;
}

void MaterialStateVariables::unpack(Vector& data) {
	int N = 0;
	bool backward = false;

	vectorvector(N, eps, data, backward);
	vectorvector(N, eps_commit, data, backward);

	vectorvector(N, sig, data, backward);
	vectorvector(N, sig_commit, data, backward);

	vectorvector(N, xs, data, backward);
	vectorvector(N, xs_commit, data, backward);

	dtime_n = data(N + 3);
	dtime_n_commit = data(N + 4);
	dtime_is_user_defined = data(N + 5);
	dtime_first_set = data(N + 6);
}

void MaterialStateVariables::printStats(bool detail) {
	if (detail) {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress                  =  " << sig;
		opserr << "Comitted Stress         =  " << sig_commit;
		opserr << "Strain                  =  " << eps;
		opserr << "Comitted Strain         =  " << eps_commit;
		opserr << "-------------------- Pastic Internal Variables --------------------\n";
		opserr << "Plastic Strain          =  " << xs;
		opserr << "Commited Plastic Strain =  " << xs_commit;
		opserr << "Plastic multiplier      =  " << lambda << "\n";
		opserr << "Commited pl. multiplier =  " << lambda_commit << "\n";
		opserr << "Prev. comm. pl. mult.   =  " << lambda_commit_old << "\n";
		opserr << "------------------ Consistent Tangent Operator ------------------\n";
		opserr << "Kt = " << Cep;
	}
	else {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress             =  " << sig;
		opserr << "Strain             =  " << eps;
		opserr << "Plastic multiplier =  " << lambda << "\n";
		opserr << "Plastic Strain     =  " << xs << "\n";
	}
}

// Private methods
	// un/pack-self
int MaterialStateVariables::getSize(void) {

	int size = 9 * 6 + 4;
	return size;
}

void MaterialStateVariables::vectorvector(int& N, Vector& A, Vector& B, bool forward) {
	if (forward) {
		for (int i = 0; i < A.Size(); i++) {
			B(N) = A(i);
			N++;
		}
	}
	else {
		for (int i = 0; i < A.Size(); i++) {
			A(i) = B(N);
			N++;
		}
	}
}

void MaterialStateVariables::vectormatrix(int& N, Matrix& A, Vector& B, bool forward) {
	int nCol = A.noCols();
	int nRow = A.noRows();
	if (forward) {
		for (int i = 0; i < nRow; i++) {
			for (int j = 0; j < nCol; j++) {
				B(N) = A(i, j);
				N++;
			}
		}
	}
	else {
		for (int i = 0; i < nRow; i++) {
			for (int j = 0; j < nCol; j++) {
				A(i, j) = B(N);
				N++;
			}
		}
	}
}