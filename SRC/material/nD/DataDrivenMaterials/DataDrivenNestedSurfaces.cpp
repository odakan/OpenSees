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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/DataDrivenNestedSurfaces.cpp$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2022
//
// Description: This file contains the implementation for the DataDrivenNestedSurfaces class.


#include "DataDrivenNestedSurfaces.h"

#define LARGE_NUMBER 1.0e30


// Public methods
	// full constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(double c0, double phi0,			// pressure-independent type constructor
	double s0, double tnys, double* pm, double* hp)
{
	// initialize variables
	cohesion = c0; frictionAngle = phi0; peakShearStrain = s0;
	TNYS = tnys;

	// allocate space for the lookup table
	HParams = Vector(TNYS + 1);
	HModuli = Vector(TNYS + 1);

	// load data into RAM
	if (pm != nullptr && hp != nullptr)
	{
		for (int i = 0; i < TNYS; i++)
		{
			HParams(i + 1) = hp[i];
			HModuli(i) = pm[i];
		}
		HParams(0) = 0.1;
	}
	else
	{
		opserr << "WARNING: DataDrivenNestedSurfaces::DataDrivenNestedSurfaces: vector parameters returned NULLPTR!\n";
		opserr << "WARNING: DataDrivenNestedSurfaces::DataDrivenNestedSurfaces: Generating automatic surfaces instead!\n";
	}

}


DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(double c0, double phi0,			// pressure-dependent type constructor
	double psi0, double s0, double tnys, double* pm, double* hp)
{

}


	// destructor
DataDrivenNestedSurfaces::~DataDrivenNestedSurfaces(void) 
{
}


	// operational methods
bool DataDrivenNestedSurfaces::canDelete(void) { return (how_many < 2); }
void DataDrivenNestedSurfaces::checkin(void) { how_many++;}
void DataDrivenNestedSurfaces::checkout(void) { how_many--;}


DataDrivenNestedSurfaces* DataDrivenNestedSurfaces::getCopy(void) 
{
	DataDrivenNestedSurfaces* copy = new DataDrivenNestedSurfaces(*this);
	return copy;
}


// setup nested yield surface package
YieldSurfacePackage* DataDrivenNestedSurfaces::setUpYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial) {

	YieldSurfacePackage* theSurface = nullptr;
	theSurface = new YieldSurfacePackage();

	generateYieldSurfaces(theMaterial);

	return theSurface;
}


// Private methods
	// generate methods
void DataDrivenNestedSurfaces::generateYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial) {
// generate plastic modulus and hardening parameter sets
// based on the hyperbolic backbone model
opserr << "generateYieldSurfaces: Works until here!\n";

double  stress1, stress2, strain1, strain2, size, Hep_ref, Href;
double refStrain, peakShear, coneHeight;


	if (theMaterial->getPhi() > 0) {
		double sinPhi = sin(theMaterial->getPhi() * M_PI / 180.);
		double Mnys = 6. * sinPhi / (3. - sinPhi);
		residualPressure = 3. * cohesion / (sqrt(2.) * Mnys);
		coneHeight = -(theMaterial->getPref() - residualPressure);
		peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
		refStrain = (theMaterial->getPeakStrain() * peakShear) / (theMaterial->getGref() * theMaterial->getPeakStrain() - peakShear);
	}

	else if (theMaterial->getPhi() == 0.0) { // cohesion = peakShearStrength
		peakShear = 2 * sqrt(2.) * cohesion / 3;
		refStrain = (theMaterial->getPeakStrain() * peakShear) / (theMaterial->getGref() * theMaterial->getPeakStrain() - peakShear);
		residualPressure = 0.;
	}

	double stressInc = peakShear / TNYS;

	for (int ii = 1; ii <= TNYS; ii++) {
		stress1 = ii * stressInc;
		stress2 = stress1 + stressInc;
		strain1 = stress1 * refStrain / (theMaterial->getGref() * refStrain - stress1);
		strain2 = stress2 * refStrain / (theMaterial->getGref() * refStrain - stress2);
		if (theMaterial->getPhi() > 0.) {
			size = 3. * stress1 / sqrt(2.) / coneHeight;
		}
		else if (theMaterial->getPhi() == 0.) {
			size = 3. * stress1 / sqrt(2.);
		}

		Hep_ref = 2. * (stress2 - stress1) / (strain2 - strain1);

		if ((2. * theMaterial->getGref() - Hep_ref) <= 0) {
			Href = LARGE_NUMBER;
		}
		else {
			Href = (2. * theMaterial->getGref() * Hep_ref) / (2. * theMaterial->getGref() - Hep_ref);
		}

		if (Href < 0) {
			Href = 0;
		}

		if (Href > LARGE_NUMBER) {
			Href = LARGE_NUMBER;
		}

		if (ii == TNYS) {
			Href = 0;
		}

		static Vector temp(6);
		//committedSurfaces[ii] = MultiYieldSurface(temp, size, Href);

	}
}