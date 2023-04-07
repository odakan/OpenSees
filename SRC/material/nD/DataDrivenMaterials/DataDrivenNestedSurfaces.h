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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/DataDrivenNestedSurfaces.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef DataDrivenNestedSurfaces_h
#define DataDrivenNestedSurfaces_h
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
#include "MultiYieldSurfaceHardeningSoftening.h"


class DataDrivenNestedSurfaces {
public:
	// null constructor
	DataDrivenNestedSurfaces(void) = default;

	// full constructors
	DataDrivenNestedSurfaces(const DataDrivenNestedSurfaces&) = default;											// copy constructor
	DataDrivenNestedSurfaces(double cohesion, double frictionAngle, double peakShearStrain,							// pressure-independent type constructor
		double tnys, double* HModuli, double* HParams);
	DataDrivenNestedSurfaces(double cohesion, double frictionAngle, double dilationAngle, double peakShearStrain,	// pressure-dependent type constructor
		double tnys, double* HModuli, double* HParams);

	// destructor
	~DataDrivenNestedSurfaces(void);

	// operator overloading
	DataDrivenNestedSurfaces& operator=(const DataDrivenNestedSurfaces&) = default;		// one-to-one assignment

	// operational methods
	bool canDelete(void);						// return true if no material is using the object
	void checkin(void);							// increase how_many counter
	void checkout(void);						// decrease how_many counter
	DataDrivenNestedSurfaces* getCopy(void);	// retun a copy of the object

	// setup yield surface
	YieldSurfacePackage* setUpYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial);

private:
	// yield surface paramters
	int TNYS = 0;
	double cohesion = 0.0;
	double frictionAngle = 0.0;
	double dilatancyAngle = 0.0;
	double peakShearStrain = 0.0;
	double residualPressure = 0.0;

	// operational variables
	int how_many = 0;							// number of materials using this surface object

	// representative volume element data
	Vector HModuli;
	Vector HParams;

private:
	// generate automatic yield surfaces
	void generateYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial);

};
#endif