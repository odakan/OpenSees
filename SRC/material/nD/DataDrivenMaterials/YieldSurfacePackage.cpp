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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/YieldSurfacePackage.cpp$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2022
//
// Description: This file contains the implementation for the YieldSurfacePackage class.


#include "YieldSurfacePackage.h"

#define LARGE_NUMBER 1.0e30

// Public methods
	// full constructors
		// traditional hyperbolic surface constructors
YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// von Mises constructor
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


YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// Drucker-Prager constructor
	double n, double c, double peakStrain, double Phi, double dilationAngle)
{

}


YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// Matsuoka-Nakai constructor
	double n, double peakStrain, double Phi)
{

}


		// data-driven surface constructors
YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// von Mises constructor
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
			opserr << "WARNING: YieldSurfacePackage::YieldSurfacePackage: vector parameters returned NULLPTR!\n";
			opserr << "WARNING: YieldSurfacePackage::YieldSurfacePackage: Generating automatic surfaces instead!\n";
			use_custom_surface = false;
		}
	}

	// generate yield surfaces
	generateYieldSurfaces();
}


YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// Drucker-Prager constructor
	double n, double c)
{

}


YieldSurfacePackage::YieldSurfacePackage(int tnys, double k, double g, double p,			// Matsuoka-Nakai constructor
	double n)
{

}


	// destructor
YieldSurfacePackage::~YieldSurfacePackage(void) 
{
}

	// operational methods
bool YieldSurfacePackage::canDelete(void) { return (how_many < 2); }
void YieldSurfacePackage::checkin(void) { how_many++;}
void YieldSurfacePackage::checkout(void) { how_many--;}

YieldSurfacePackage* YieldSurfacePackage::getCopy(void) 
{
	YieldSurfacePackage* copy = new YieldSurfacePackage(*this);
	return copy;
}

	// update methods
void YieldSurfacePackage::updateTNYS(int var) { TNYS = var; }
void YieldSurfacePackage::updateKref(double var) { Kref = var; }
void YieldSurfacePackage::updateGref(double var) { Gref = var; }
void YieldSurfacePackage::updatePref(double var) { Pref = var; }
void YieldSurfacePackage::updateModn(double var) { modn = var; }
void YieldSurfacePackage::updatePhi(double var) { Phi = var; }
void YieldSurfacePackage::updatePsi(double var) { Psi = var; }
void YieldSurfacePackage::updateCohesion(double var) { cohesion = var; }
void YieldSurfacePackage::updateHardParams(Vector& var) { HardParams = var; }
void YieldSurfacePackage::updateHardParams(double var, int nYs_active) { HardParams(nYs_active) = var; }
void YieldSurfacePackage::updateDilatParams(Vector& var) { DilatParams = var; }
void YieldSurfacePackage::updateDilatParams(double var, int nYs_active) { DilatParams(nYs_active) = var; }

// get methods
int YieldSurfacePackage::getTNYS(void) { return TNYS; }
double YieldSurfacePackage::getKref(void) { return Kref; }
double YieldSurfacePackage::getGref(void) { return Gref; }
double YieldSurfacePackage::getPref(void) { return Pref; }
double YieldSurfacePackage::getModn(void) { return modn; }
double YieldSurfacePackage::getPhi(void) { return Phi; }
double YieldSurfacePackage::getPsi(void) { return Psi; }
double YieldSurfacePackage::getCohesion(void) { return cohesion; }
double YieldSurfacePackage::getHref(int num) { return Href(num); }
double YieldSurfacePackage::getHP(int num) { return HardParams(num); }
double YieldSurfacePackage::getDP(int num) { return DilatParams(num); }

	// generate methods
void YieldSurfacePackage::generateYieldSurfaces(void)
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
		if (Phi > 0) {
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

// Private methods