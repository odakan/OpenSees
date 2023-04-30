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
	// destructor
MaterialStateVariables::~MaterialStateVariables(void)
{
	/* 
	   do nothing!
	   no memory was allocated in the memory-land to free... 
	*/
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
	vectorvector(N, eps_commit_old, data, forward);

	vectorvector(N, sig, data, forward);
	vectorvector(N, sig_commit, data, forward);
	vectorvector(N, sig_commit_old, data, forward);

	vectorvector(N, xs, data, forward);
	vectorvector(N, xs_commit, data, forward);
	vectorvector(N, xs_commit_old, data, forward);

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
	vectorvector(N, eps_commit_old, data, backward);

	vectorvector(N, sig, data, backward);
	vectorvector(N, sig_commit, data, backward);
	vectorvector(N, sig_commit_old, data, backward);

	vectorvector(N, xs, data, backward);
	vectorvector(N, xs_commit, data, backward);
	vectorvector(N, xs_commit_old, data, backward);

	dtime_n = data(N + 3);
	dtime_n_commit = data(N + 4);
	dtime_is_user_defined = data(N + 5);
	dtime_first_set = data(N + 6);
}

void MaterialStateVariables::printStats(bool detail) {
	if (detail) {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress                    =  " << sig;
		opserr << "Comitted Stress           =  " << sig_commit;
		opserr << "Strain                    =  " << eps;
		opserr << "Comitted Strain           =  " << eps_commit;
		opserr << "-------------------- Pastic Internal Variables --------------------\n";
		opserr << "Plastic Strain            =  " << xs;
		opserr << "Commited Plastic Strain   =  " << xs_commit;
		opserr << "------------------ Consistent Tangent Operator ------------------\n";
		opserr << "Kt = " << Cep;
	}
	else {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress                    =  " << sig;
		opserr << "Strain                    =  " << eps;
		opserr << "Plastic Strain            =  " << xs << "\n";
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