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
		// traditional hyperbolic surface constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// von Mises constructor
	double n, double c, double peakStrain, double Phi)
{
	// initialize variables
	TNYS = tnys; Kref = k; Gref = g;
	Pref = p; modn = n; cohesion = c;
	Phi = Phi; peakStrain = peakStrain;
	use_custom_surface = false;
	
	// generate yield surfaces
	generateYieldSurfaces();
}


DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// Drucker-Prager constructor
	double n, double c, double peakStrain, double Phi, double dilationAngle)
{

}


DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// Matsuoka-Nakai constructor
	double n, double peakStrain, double Phi)
{

}


		// data-driven surface constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// von Mises constructor
	double n, double c, double* gp, double* hp)
{
	// initialize variables
	TNYS = tnys; Kref = k; Gref = g;
	Pref = p; modn = n; cohesion = c;
	use_custom_surface = true;

	// handle user defined surfaces
	if (TNYS < 0)
	{
		use_custom_surface = true;
		TNYS = abs(TNYS);
		// allocate space for vector parameters
		HardParams = Vector(TNYS + 1);
		Href = Vector(TNYS + 1);
		// store parameters
		if (hp != nullptr && gp != nullptr)
		{
			for (int i = 0; i < TNYS; i++)
			{
				HardParams(i + 1) = hp[i];
				Href(i) = gp[i];
			}
			HardParams(0) = 0.1;
		}
		else
		{
			opserr << "WARNING: DataDrivenNestedSurfaces::DataDrivenNestedSurfaces: vector parameters returned NULLPTR!\n";
			opserr << "WARNING: DataDrivenNestedSurfaces::DataDrivenNestedSurfaces: Generating automatic surfaces instead!\n";
			use_custom_surface = false;
		}
	}

	// generate yield surfaces
	generateYieldSurfaces();
}


DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// Drucker-Prager constructor
	double n, double c)
{

}


DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// Matsuoka-Nakai constructor
	double n)
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


	// update methods
void DataDrivenNestedSurfaces::updateHardParams(Vector& var) { HardParams = var; }
void DataDrivenNestedSurfaces::updateHardParams(double var, int nYs_active) { HardParams(nYs_active) = var; }
void DataDrivenNestedSurfaces::updateDilatParams(Vector& var) { DilatParams = var; }
void DataDrivenNestedSurfaces::updateDilatParams(double var, int nYs_active) { DilatParams(nYs_active) = var; }


// get methods
double DataDrivenNestedSurfaces::getHref(int num) { return Href(num); }
double DataDrivenNestedSurfaces::getHP(int num) { return HardParams(num); }
double DataDrivenNestedSurfaces::getDP(int num) { return DilatParams(num); }


// setup nested yield surface package
YieldSurfacePackage* DataDrivenNestedSurfaces::setUpYieldSurfaces(VonMisesDMM* theMaterial) {

	YieldSurfacePackage* theSurface = nullptr;
	theSurface = new YieldSurfacePackage();
	generateYieldSurfaces(theMaterial);

	return theSurface;
}


YieldSurfacePackage* DataDrivenNestedSurfaces::setUpYieldSurfaces(DruckerPragerDMM* theMaterial) {
	YieldSurfacePackage* theSurface = nullptr;


	return theSurface;
}


YieldSurfacePackage* DataDrivenNestedSurfaces::setUpYieldSurfaces(MatsuokaNakaiDMM* theMaterial) {
	YieldSurfacePackage* theSurface = nullptr;


	return theSurface;
}


// Private methods
	// generate methods
void DataDrivenNestedSurfaces::generateYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial)
{
	// generate plastic modulus and hardening parameter sets
	// based on the hyperbolic backbone model
	opserr << "generateYieldSurfaces: Works until here!\n";

	double  stress1, stress2, strain1, strain2, size, Hep_ref, Href;
	double refStrain, peakShear, coneHeight;

	if (use_custom_surface)
	{

	}
	else
	{
		if (theMaterial->Phi > 0) {
			double sinPhi = sin(Phi * M_PI / 180.);
			double Mnys = 6. * sinPhi / (3. - sinPhi);
			residualPressure = 3. * cohesion / (sqrt(2.) * Mnys);
			coneHeight = -(Pref - residualPressure);
			peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
			refStrain = (peakStrain * peakShear) / (Gref * peakStrain - peakShear);
		}

		else if (Phi == 0.0) { // cohesion = peakShearStrength
			peakShear = 2 * sqrt(2.) * cohesion / 3;
			refStrain = (peakStrain * peakShear) / (Gref * peakStrain - peakShear);
			residualPressure = 0.;
		}

		double stressInc = peakShear / TNYS;

		for (int ii = 1; ii <= TNYS; ii++) {
			stress1 = ii * stressInc;
			stress2 = stress1 + stressInc;
			strain1 = stress1 * refStrain / (Gref * refStrain - stress1);
			strain2 = stress2 * refStrain / (Gref * refStrain - stress2);
			if (Phi > 0.) {
				size = 3. * stress1 / sqrt(2.) / coneHeight;
			}
			else if (Phi == 0.) {
				size = 3. * stress1 / sqrt(2.);
			}

			Hep_ref = 2. * (stress2 - stress1) / (strain2 - strain1);

			if ((2. * Gref - Hep_ref) <= 0) {
				Href = LARGE_NUMBER;
			}
			else {
				Href = (2. * Gref * Hep_ref) / (2. * Gref - Hep_ref);
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

}