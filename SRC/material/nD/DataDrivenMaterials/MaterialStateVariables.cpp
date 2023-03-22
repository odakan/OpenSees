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
	// full constructors
MaterialStateVariables::MaterialStateVariables(int alpha_size) 
{
	// resize the alpha matrix
	alpha.resize(6, (alpha_size + 1));
	alpha_commit.resize(6, (alpha_size + 1));
	alpha_commit_old.resize(6, (alpha_size + 1));
	// fill with zeros
	alpha.Zero();
	alpha_commit.Zero();
	alpha_commit_old.Zero();
}

	// destructor
MaterialStateVariables::~MaterialStateVariables(void)
{
	/* 
	   do nothing!
	   no memory was allocated in the memory-land to free... 
	*/
}

	// operator overloading
void MaterialStateVariables::operator+= (const MaterialStateVariables& A)
{
	if (A.alpha.noCols() == alpha.noCols()) {
		MaterialStateVariables B = MaterialStateVariables(alpha.noCols());
		alpha += A.alpha;
		// ...
	}
	else {
		opserr << "MaterialStateVariables::operator+= object alpha matrix sizes must be equal!\n";
		exit(-1);
	}
}

void MaterialStateVariables::operator-= (const MaterialStateVariables& A)
{
	if (A.alpha.noCols() == alpha.noCols()) {
		alpha -= A.alpha;
		// ...
	}
	else {
		opserr << "MaterialStateVariables::operator-= object alpha matrix sizes must be equal!\n";
		exit(-1);
	}
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

	vectormatrix(N, alpha, data, forward);
	vectormatrix(N, alpha_commit, data, forward);
	vectormatrix(N, alpha_commit_old, data, forward);

	data(N) = nYs_active;
	data(N + 1) = nYs_active_commit;
	data(N + 2) = nYs_active_commit_old;

	data(N + 3) = dtime_n;
	data(N + 4) = dtime_n_commit;
	data(N + 5) = dtime_is_user_defined;
	data(N + 6) = dtime_first_set;
	return data;
}

void MaterialStateVariables::unpack(Vector& data) {
	int N = 0;
	int alphasz = getTnys(data);
	bool backward = false;
	alpha = Matrix(6, alphasz);
	alpha_commit = Matrix(6, alphasz);
	alpha_commit_old = Matrix(6, alphasz);

	vectorvector(N, eps, data, backward);
	vectorvector(N, eps_commit, data, backward);
	vectorvector(N, eps_commit_old, data, backward);

	vectorvector(N, sig, data, backward);
	vectorvector(N, sig_commit, data, backward);
	vectorvector(N, sig_commit_old, data, backward);

	vectorvector(N, xs, data, backward);
	vectorvector(N, xs_commit, data, backward);
	vectorvector(N, xs_commit_old, data, backward);

	vectormatrix(N, alpha, data, backward);
	vectormatrix(N, alpha_commit, data, backward);
	vectormatrix(N, alpha_commit_old, data, backward);

	nYs_active = data(N);
	nYs_active_commit = data(N + 1);
	nYs_active_commit_old = data(N + 2);

	dtime_n = data(N + 3);
	dtime_n_commit = data(N + 4);
	dtime_is_user_defined = data(N + 5);
	dtime_first_set = data(N + 6);
}

void MaterialStateVariables::printStats(bool detail) {
	if (detail) {
		opserr << "\n";
		opserr << "MaterialStateVariables::printStats:\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "\n";
		opserr << "Stress                    =  " << sig << "\n";
		opserr << "Comitted Stress           =  " << sig_commit << "\n";
		opserr << "\n";
		opserr << "Strain                    =  " << eps << "\n";
		opserr << "Comitted Strain           =  " << eps_commit << "\n";
		opserr << "\n";
		opserr << "-------------------- Pastic Internal Variables --------------------\n";
		opserr << "\n";
		opserr << "Plastic Strain            =  " << xs << "\n";
		opserr << "Commited Plastic Strain   =  " << xs_commit << "\n";
		opserr << "\n";
		opserr << "Back-stress               =  " << tools::getColumnVector(nYs_active, alpha) << "\n";
		opserr << "Commited Back-stress      =  " << tools::getColumnVector(nYs_active_commit, alpha_commit) << "\n";
		opserr << "\n";
		opserr << "Active Y-Surface          =  " << nYs_active << "\n";
		opserr << "Commited Active Y-Surface =  " << nYs_active_commit << "\n";
		opserr << "\n";
		opserr << "------------------ Consistent Tangent Operator ------------------\n";
		opserr << "\n";
		opserr << "Kt = " << Cep << "\n";
	}
	else {
		opserr << "\n";
		opserr << "MaterialStateVariables::printStats:\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress                    =  " << sig << "";
		opserr << "Strain                    =  " << eps << "";
		opserr << "Plastic Strain            =  " << xs << "";
		opserr << "Back-stress               =  " << tools::getColumnVector(nYs_active, alpha) << "";
		opserr << "Active Y-Surface          =  " << nYs_active << "\n";
		opserr << "\n\n";
	}
}

// Private methods
	// un/pack-self
int MaterialStateVariables::getSize(void) {
	int size = 0;
	int sz_alpha = 0;
	try {
		sz_alpha = alpha.noCols();
	}
	catch (const std::exception&) {
		opserr << "WARNING: MaterialStateVariables are not initialized! --> Assuming sz_alpha = 0!";
	}
	size = 7 + 9 * 6 + 3 * 6 * (sz_alpha);
	return size;
}

int MaterialStateVariables::getTnys(const Vector& data) {
	int size = 0;
	int tnys = 0;
	try {
		size = data.Size();
	}
	catch (const std::exception&) {
		opserr << "WARNING: MaterialStateVariables vector is not initialized! --> Assuming tnys = 0!";
	}
	tnys = tools::macaulay(((size - 8 + 12 * 6) / (4 * 6)));
	return tnys;
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