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


#ifndef MaterialStateVariables_h
#define MaterialStateVariables_h


#include <map>
#include <stdio.h>
#include <stdlib.h> 

#include "CTensor.h"


class MaterialStateVariables {
public:
	// state variables
	/*	IMPLEX variables:
			lambda	: lambda-bar, rephrased (linearized) plastic flow, in the explicit step, lambda is
					extrapolated using lambda_commit and lambda_commit_old

			ksi		: ksi_bar, rephrased (linearized) translation direction, in the explicit step, ksi
					is extrapolated using ksi_commit and ksi_commit_old

			dtime_n:

			dtime_is_user_defined:

			dtime_first_set:

			sig_implex: stress tensor computed at the end of the
						explicit step is stored for output
	*/
	CTensor eps = CTensor(6, 1);			// material strain (1: covariant)
	CTensor eps_commit = CTensor(6, 1);		// committed material strain (1: covariant)
	CTensor sig = CTensor(6, 2);			// effective stress (2: contravariant)
	CTensor sig_commit = CTensor(6, 2);		// committed effective stress (2: contravariant)
	CTensor sig_implex = CTensor(6, 2);		// only for output (2: contravariant)
	CTensor xs = CTensor(6, 1);				// equivalent plastic strain (1: covariant)
	CTensor xs_commit = CTensor(6, 1);		// committed equivalent plastic strain (1: covariant)
	double lambda = 0.0;					// plastic multiplier (linearized)
	double lambda_commit = 0.0;				// committed plastic multiplier
	double lambda_commit_old = 0.0;			// previously committed plastic multiplier
	double ksi = 0.0;						// surface translation direction (linearized)
	double ksi_commit = 0.0;				// committed translation direction
	double ksi_commit_old = 0.0;			// previously committed translation direction
	double dtime_n = 0.0;					// time factor
	double dtime_n_commit = 0.0;			// committed time factor
	bool dtime_is_user_defined = false;
	bool dtime_first_set = false;
	// moduli
	CTensor Ce = CTensor(6, 6, 2);			// elastic modulus matrix (2: contravariant)
	CTensor Cep = CTensor(6, 6, 2);			// tangent modulus matrix (2: contravariant)
	double Kmod = 0;						// Current bulk modulus
	double Gmod = 0;						// Current shear modulus
	double Hmod = 0;						// Current plastic modulus

public:
	// null constructor
	MaterialStateVariables(void) = default;

	// full constructors
	MaterialStateVariables(const double nDim);
	MaterialStateVariables(const MaterialStateVariables&) = default;

	// destructor
	~MaterialStateVariables(void);

	// iteration control
	int commitState(void);
	int revertToLastCommit(void);

	// operator overloading
	MaterialStateVariables& operator= (const MaterialStateVariables&) = default;
	void printStats(bool detail);
};
#endif