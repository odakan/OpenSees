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
MaterialStateVariables::MaterialStateVariables(const double nDim)
{
	if (nDim == 2) {
		eps = CTensor(3, 1);
		eps_commit = CTensor(3, 1);
		sig = CTensor(3, 2);
		sig_commit = CTensor(3, 2);
		sig_implex = CTensor(3, 2);
		xs = CTensor(3, 1);
		xs_commit = CTensor(3, 1);
		Ce = CTensor(3, 3, 2);
		Cep = CTensor(3, 3, 2);
	}
	else if (nDim == 3) {
		eps = CTensor(6, 1);			// 1: covariant
		eps_commit = CTensor(6, 1);		// 1: covariant
		sig = CTensor(6, 2);			// 2: contravariant
		sig_commit = CTensor(6, 2);		// 2: contravariant
		sig_implex = CTensor(6, 2);		// 2: contravariant
		xs = CTensor(6, 1);				// 1: covariant
		xs_commit = CTensor(6, 1);		// 1: covariant
		Ce = CTensor(6, 6, 2);			// 2: contravariant
		Cep = CTensor(6, 6, 2);			// 2: contravariant
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
	ksi_commit_old = ksi_commit;
	ksi_commit = ksi;

	return 0;
}

int MaterialStateVariables::revertToLastCommit(void) {

	// restore committed internal variables
	eps = eps_commit;
	sig = sig_commit;
	xs = xs_commit;
	lambda = lambda_commit;
	lambda_commit = lambda_commit_old;
	ksi = ksi_commit;
	ksi_commit = ksi_commit_old;

	return 0;
}

void MaterialStateVariables::printStats(bool detail) {
	if (detail) {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress                  =  " << sig;
		opserr << "Comitted stress         =  " << sig_commit;
		opserr << "Strain                  =  " << eps;
		opserr << "Comitted strain         =  " << eps_commit;
		opserr << "-------------------- Pastic Internal Variables --------------------\n";
		opserr << "Plastic strain          =  " << xs;
		opserr << "Commited plastic strain =  " << xs_commit;
		opserr << "Plastic multiplier      =  " << lambda << "\n";
		opserr << "Commited pl. multiplier =  " << lambda_commit << "\n";
		opserr << "Prev. comm. pl. mult.   =  " << lambda_commit_old << "\n";
		opserr << "Translation multiplier  =  " << ksi << "\n";
		opserr << "Commited tr. multiplier =  " << ksi_commit << "\n";
		opserr << "Prev. comm. tr. mult.   =  " << ksi_commit_old << "\n";
		opserr << "------------------ Consistent Tangent Operator ------------------\n";
		opserr << "Kt = " << Cep;
	}
	else {
		opserr << "MaterialStateVariables::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Stress             =  " << sig;
		opserr << "Strain             =  " << eps;
		opserr << "Plastic multiplier =  " << lambda << "\n";
		opserr << "Translation mult.  =  " << ksi << "\n";
		opserr << "Plastic strain     =  " << xs << "\n";
	}
}