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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/YieldSurfacePackage.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef YieldSurfacePackage_h
#define YieldSurfacePackage_h
#define _USE_MATH_DEFINES

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2022
//

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include "YieldSurfacePackage.h"


class YieldSurfacePackage {
public:
	// null constructor
	YieldSurfacePackage(void) = default;

	// full constructors
	YieldSurfacePackage(const YieldSurfacePackage&) = default;				// default constructor

	// destructor
	~YieldSurfacePackage(void);

	// operator overloading
	YieldSurfacePackage& operator=(const YieldSurfacePackage&) = default;	// one-to-one assignment

	// operational methods
	YieldSurfacePackage* getCopy(void);										// retun a copy of the object

	// update methods
	void updateHardParams(Vector& var);
	void updateHardParams(double var, int num);
	void updateDilatParams(Vector& var);
	void updateDilatParams(double var, int num);

	// get methods
	double getHref(int num);
	double getHP(int num);
	double getDP(int num);

private:
	// operational parameters
	bool do_online = false;

	// yield surface variables
	double theSize;
	double theModulus;
	Vector theCenter;

};
#endif