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
// Created in:	April 2023
//
// Description: This file contains the implementation for the DataDrivenNestedSurfaces class.


#include "DataDrivenNestedSurfaces.h"


#define LARGE_NUMBER 1.0e30

const bool DEBUG = true;

// Public methods
	// full constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tag, double c0, double phi0,			// full constructor
	double psi0, double s0, double nys, double* pm, double* hp, bool verbosity)
{
	// initialize variables
	beVerbose = verbosity;
	cohesion_init = c0; frictionAngle_init = phi0; peakShearStrain_init = s0;
	tnys_init = nys;

	// evaluate if yield surface generation is possible
		// check data-driven yield surface generation
	if (pm != nullptr && hp != nullptr) {
		isOfflineOK = true;
		isOnlineOK = true;
	}
	else {
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - no material response database is provided. Data-driven yield surfaces are NOT allowed...\n"; }
		isOfflineOK = false;
		isOnlineOK = false;
	}

		// check automatic and offline surface generation
	if (tnys_init <= 0) {
		// if total number of yield surfaces is less than one, prohibit automatic and offline surfaces
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - TNYS <= 0: AUTOMATIC and OFFLINE data-driven yield surfaces are NOT allowed...\n"; }
		isAutomaticOK = false;
		isOfflineOK = false;
	}
	else {
		// check automatic surfaces
		if (peakShearStrain_init <= 0) {
			// if peak shear strain is zero or less than zero, prohibit automatic surfaces
			if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - peakShearStrain <= 0: AUTOMATIC yield surfaces are NOT allowed...\n"; }
			isAutomaticOK = false;;
		}
		else {
			if (frictionAngle_init > 0) {
				isAutomaticOK = true;
			}
			else if (cohesion_init > 0) {
				isAutomaticOK = true;
			}
			else {
				if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - either cohesion or friction angle must be > 0: AUTOMATIC yield surfaces are NOT allowed...\n"; }
				isAutomaticOK = false;
			}
		}
	}

	if (isOfflineOK) {
		// allocate space for the lookup table
		HParams = Vector(tnys_init + 1);
		HModuli = Vector(tnys_init + 1);

		// load data into RAM
		if (pm != nullptr && hp != nullptr)
		{
			for (int i = 0; i < tnys_init; i++)
			{
				HParams(i + 1) = hp[i];
				HModuli(i) = pm[i];
			}
			HParams(0) = 0.1;
		}
	}
}


	// destructor
DataDrivenNestedSurfaces::~DataDrivenNestedSurfaces(void) 
{
}


	// operational methods
bool DataDrivenNestedSurfaces::canDelete(void) { return (how_many < 2); }
void DataDrivenNestedSurfaces::checkin(void) { how_many++;}
void DataDrivenNestedSurfaces::checkout(void) { how_many--;}

bool DataDrivenNestedSurfaces::isAOK(int dataDriver) {
	
	bool aok = false;

	if (dataDriver == 0 && isAutomaticOK) {
		aok = true;
	}
	else if (dataDriver == 1 && isOfflineOK) {
		aok = true;
	}
	else if (dataDriver == 2 && isOnlineOK) {
		aok = true;
	}
	else {
		opserr << "WARNING: DataDrivenNestedSurfaces::isAOK() - none of the yield surface generation approaches is available!\n";
	}

	return aok;
}

DataDrivenNestedSurfaces* DataDrivenNestedSurfaces::getCopy(void)  {

	DataDrivenNestedSurfaces* copy = new DataDrivenNestedSurfaces(*this);
	return copy;
}


// get methods
int DataDrivenNestedSurfaces::getTNYS(void) { return tnys_init; }
double DataDrivenNestedSurfaces::getCohesion(void) { return cohesion_init; }
double DataDrivenNestedSurfaces::getPhi(void) { return frictionAngle_init; }
double DataDrivenNestedSurfaces::getPsi(void) { return dilatancyAngle_init; }
double DataDrivenNestedSurfaces::getPeakStrain(void) { return peakShearStrain_init; }
double DataDrivenNestedSurfaces::getPref(void) { return referencePressure_init; }


// generate nested yield surface package
YieldSurfacePackage DataDrivenNestedSurfaces::generateYieldSurfaces(const int matid, const int dataDriver, const double Pref, const double Gref) {
	
	// initialize the yield surface package
	YieldSurfacePackage yieldSurface;
	
	// set up surfaces
	if (dataDriver == 2) {
		yieldSurface = YieldSurfacePackage(matid);
		setUpOnlineSurfaces(yieldSurface);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	else if (dataDriver == 1) {
		yieldSurface = YieldSurfacePackage(matid, tnys_init);
		setUpOfflineSurfaces(yieldSurface);
		if (beVerbose) { yieldSurface.printStats(true); }
	} 
	else {
		yieldSurface = YieldSurfacePackage(matid, tnys_init, cohesion_init, frictionAngle_init,
			dilatancyAngle_init, peakShearStrain_init, residualPressure_init, referencePressure_init);
		setUpAutomaticSurfaces(yieldSurface, Pref, Gref, tnys_init);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	return yieldSurface;
}


// Private methods
	// set up yield surfaces
void DataDrivenNestedSurfaces::setUpOnlineSurfaces(YieldSurfacePackage& theSurface) {




}

void DataDrivenNestedSurfaces::setUpOfflineSurfaces(YieldSurfacePackage& theSurface) {



}

void DataDrivenNestedSurfaces::setUpAutomaticSurfaces(YieldSurfacePackage& yieldSurface, const double Pref, const double Gref, const int tnys) {
	// set up plastic modulus and hardening parameter sets based on the hyperbolic backbone model

	double  stress1, stress2, strain1, strain2, size, Hep_ref, Href;
	double refStrain, peakShear, coneHeight;
	
	if (yieldSurface.getPhi() > 0) {
		double sinPhi = sin(yieldSurface.getPhi() * M_PI / 180.);
		double MnyieldSurface = 6. * sinPhi / (3. - sinPhi);
		yieldSurface.setPresid(3. * yieldSurface.getCohesion() / (sqrt(2.) * MnyieldSurface));
		coneHeight = -(Pref - yieldSurface.getPresid());
		peakShear = sqrt(2.) * coneHeight * MnyieldSurface / 3.;
		refStrain = (yieldSurface.getPeakStrain() * peakShear) / (Gref * yieldSurface.getPeakStrain() - peakShear);
	}
	else if (yieldSurface.getPhi() == 0.0) {			// cohesion = peakShearStrength
		peakShear = 2 * sqrt(2.) * yieldSurface.getCohesion() / 3;
		refStrain = (yieldSurface.getPeakStrain() * peakShear) / (Gref * yieldSurface.getPeakStrain() - peakShear);
		yieldSurface.setPresid(0.0);
	}

	double stressInc = peakShear / tnys;

	for (int ii = 0; ii <= tnys; ii++) {
		stress1 = ii * stressInc;
		stress2 = stress1 + stressInc;
		strain1 = stress1 * refStrain / (Gref * refStrain - stress1);
		strain2 = stress2 * refStrain / (Gref * refStrain - stress2);
		if (yieldSurface.getPhi() > 0.) {
			size = 3. * stress1 / sqrt(2.) / coneHeight;
		}
		else if (yieldSurface.getPhi() == 0.) {
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

		if (ii == tnys) {
			Href = 0;
		}

		if (ii == 0) {
			size = stressInc * 0.1;
		}

		yieldSurface.setTau(size, ii);
		yieldSurface.setEta(Href, ii);
	}
}