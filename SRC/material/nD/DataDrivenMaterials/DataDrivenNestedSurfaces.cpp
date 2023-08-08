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


// Public methods
	// full constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tag, double c0, double phi0,			// full constructor
	double psi0, double s0, double nys, double deps, Vector hmoduli, Vector hparams, Vector dparams, bool verbosity)
{
	// initialize variables
	beVerbose = verbosity;
	cohesion_init = c0; frictionAngle_init = phi0; dilatancyAngle_init = psi0;  peakShearStrain_init = s0;
	tnys_init = nys;
	if (deps > 0.0) { strain_discretize = deps; }

	// check data-driven yield surface generation 
	if (false) {
		// if database is present
		isOfflineOK = true;
		isOnlineOK = true;
	}
	else {
		// if database is not present
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - no material response database is provided. Data-driven yield surfaces are NOT allowed...\n"; }
		isOfflineOK = true;
		isOnlineOK = true;
	}

		// check automatic, user-custom and offline surface generation
	if (tnys_init <= 0) {
		// if total number of yield surfaces is less than one, prohibit automatic and offline surfaces
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - TNYS <= 0: AUTOMATIC and OFFLINE data-driven yield surfaces are NOT allowed...\n"; }
		isUserCustomOK = false;
		isAutomaticOK = false;
		isOfflineOK = false;
	}
	else {
		// evaluate if yield surface generation is possible
		// check user custom yield surface generation
		int sz_hmod = hmoduli.Size();
		int sz_hps = hparams.Size();
		int sz_dps = dparams.Size();
		if (sz_hmod == nys && sz_hps == nys) {
			// if the user defined custom data
			isUserCustomOK = true;
			HModuli = Vector(nys);
			HParams = Vector(nys);
			for (int i = 0; i < nys; i++) {
				HModuli(i) = hmoduli(i);
				HParams(i) = hparams(i);
			}

			if (sz_dps == nys) {
				isNonassociatedOK = true;
				DParams = Vector(nys);
				for (int i = 0; i < nys; i++) {
					DParams(i) = dparams(i);
				}
			}
		}
		else {
			// if the user hasn't defined custom data
			if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - no custom yield surface data is provided. User custom yield surfaces are NOT allowed...\n"; }
			isUserCustomOK = false;
		}

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

	if (isOfflineOK || isOnlineOK) {
		// load the database into the RAM
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

	if (dataDriver == -1 && isUserCustomOK) {
		aok = true;
	}
	else if (dataDriver == 0 && isAutomaticOK) {
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
YieldSurfacePackage DataDrivenNestedSurfaces::generateYieldSurfaces(const int matid, const int dataDriver, double &Gref, double &Pref, double& Modn) {
	
	// initialize the yield surface package
	YieldSurfacePackage yieldSurface;
	
	// set up surfaces
	if (dataDriver > 1) {
		// online approach
		yieldSurface = YieldSurfacePackage(matid, strain_discretize);
		setUpOnlineSurfaces(yieldSurface, strain_discretize);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	else if (dataDriver == 1) {
		// offline approach
		yieldSurface = YieldSurfacePackage(matid, tnys_init);
		setUpOfflineSurfaces(yieldSurface);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	else if (dataDriver == 0) {
		// automatic surface generation
		yieldSurface = YieldSurfacePackage(matid, tnys_init, cohesion_init, frictionAngle_init,
			dilatancyAngle_init, peakShearStrain_init, residualPressure_init, referencePressure_init);
		setUpAutomaticSurfaces(yieldSurface, Gref, Pref);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	else if (dataDriver == -1) {
		// user custom surface generation
		yieldSurface = YieldSurfacePackage(matid, tnys_init, HParams, HModuli, DParams);
		setUpUserCustomSurfaces(yieldSurface, Gref, Pref, Modn);
		if (beVerbose) { yieldSurface.printStats(true); }
	}
	else {

	}

	return yieldSurface;
}


// Private methods
	// set up yield surfaces
void DataDrivenNestedSurfaces::setUpOnlineSurfaces(YieldSurfacePackage& yieldSurface, double deps) {




}

void DataDrivenNestedSurfaces::setUpOfflineSurfaces(YieldSurfacePackage& yieldSurface) {



}

void DataDrivenNestedSurfaces::setUpUserCustomSurfaces(YieldSurfacePackage& yieldSurface, const double Gref, const double Pref,  double &Modn) {
	// set up plastic modulus and hardening parameter sets based on the user input strain and G/Gmax pairs

	double  stress1(0.0), stress2(0.0), strain1(0.0),   strain2(0.0),   size(0.0),       dilatancy(0.0),
		    Hep_ref(0.0), Href(0.0),    refStrain(0.0), peakShear(0.0), coneHeight(0.0), strain_vol1(0.0), strain_vol2(0.0);


	// do some check and regulations
	if (yieldSurface.getPsi() > 0) {   // ignore user defined friction angle
		double tmax = 0.0; for (int i = 0; i < yieldSurface.getTNYS(); i++) { tmax = fmax(tmax, Gref * HModuli(i) * HParams(i)); }
		double Mnys = -(sqrt(3.0) * tmax - 2. * yieldSurface.getCohesion()) / Pref;
		if (Mnys <= 0) {   // also ignore user defined cohesion
			yieldSurface.setCohesion(sqrt(3.0) / 2.0 * tmax);
			yieldSurface.setPsi(0.0);
			coneHeight = 1.;
			yieldSurface.setPresid(0.0);
		}
		else {
			double sinPhi = 3 * Mnys / (6 + Mnys);
			if (sinPhi < 0. || sinPhi>1.) {
				opserr << "FATAL: DataDrivenNestedSurfaces::setUpUserCustomSurfaces() - invalid friction angle for nDmaterial "
					<< yieldSurface.getTag() << ", please modify reference pressure or G/Gmax curve.\n";
				exit(-1);
			}
			yieldSurface.setPresid(2. * yieldSurface.getCohesion() / Mnys);
			if (yieldSurface.getPresid() < 0.01 * Pref) yieldSurface.setPresid(0.01 * Pref);
			coneHeight = -(Pref - yieldSurface.getPresid());
			yieldSurface.setPsi(asin(sinPhi) * 180 / M_PI);
		}
	}
	else {   // ignore user defined cohesion
		double tmax = 0.0; for (int i = 0; i < yieldSurface.getTNYS(); i++) { tmax = fmax(tmax, Gref * HModuli(i) * HParams(i)); }
		yieldSurface.setCohesion(sqrt(3.) / 2 * tmax);
		coneHeight = 1.0;
		yieldSurface.setPresid(0.0);
		Modn = 0.0; // also ignore user defined pressDependCoeff
	}

	/*
	 *	Required user-input parameters for each yield surface
	 *	-----------------------------------------------------------------------------------------------------------------------------------
	*	HParams: -> Octahedral shear strain = 2.0 * sqrt(2.0 / 3.0 * J2') where J2' is the second invariant of the strain deviator tensor 
	*	HModuli: -> Secant shear modulus reduction curve (Gsec / Gmax)
	*	DParams: -> Secant dilatancy curve (Evol / Eoct) where Evol and Eoct are the volumetric and octahedral shear strains, respectively
	*/

	// first yield surface
	if (yieldSurface.getPhi() > 0.) {
		size = Gref * HModuli(0) * HParams(0) / coneHeight;
	}
	else if (yieldSurface.getPhi() == 0.) {
		size =  Gref * HModuli(0) * HParams(0);
	}
	Href = Gref * HModuli(0) * HParams(0) / HParams(0);
	if (Href > LARGE_NUMBER) Href = LARGE_NUMBER;
	yieldSurface.setTau(size * 0.01, 0);
	yieldSurface.setEta(Href, 0);
	if (yieldSurface.isNonAssociated()) {
		yieldSurface.setBeta(DParams(0), 0);
	}

	// last yield surface
	if (yieldSurface.getPhi() > 0.) {
		size = Gref * HModuli(yieldSurface.getTNYS() - 1) * HParams(yieldSurface.getTNYS() - 1) / coneHeight;
	}
	else if (yieldSurface.getPhi() == 0.) {
		size = Gref * HModuli(yieldSurface.getTNYS() - 1) * HParams(yieldSurface.getTNYS() - 1);
	}
	Href = 0;
	yieldSurface.setTau(size, yieldSurface.getTNYS());
	yieldSurface.setEta(Href, yieldSurface.getTNYS());
	if (yieldSurface.isNonAssociated()) {
		dilatancy = 0;  // assume constant volume deformation
		yieldSurface.setBeta(dilatancy, yieldSurface.getTNYS());
	}

	// other yield surfaces
	for (int i = 1; i < (yieldSurface.getTNYS()); i++) {
		strain1 = HParams(i-1);
		stress1 = Gref * HModuli(i - 1) * HParams(i - 1);
		strain2 = HParams(i);
		stress2 = Gref * HModuli(i) * HParams(i);
		if (yieldSurface.getPhi() > 0.) {
			size = stress1 / coneHeight;
		}
		else if (yieldSurface.getPhi() == 0.) {
			size = stress1;
		}
		Href = (stress2 - stress1) / (strain2 - strain1);
		if (Href > LARGE_NUMBER) Href = LARGE_NUMBER;
		yieldSurface.setTau(size, i);
		yieldSurface.setEta(Href, i);

		if (yieldSurface.isNonAssociated()) {
			strain_vol1 = DParams(i - 1) * strain1;
			strain_vol2 = DParams(i) * strain2;
			dilatancy = (strain_vol2 - strain_vol1) / (strain2 - strain1);
			yieldSurface.setBeta(dilatancy, i);
		}
	}
}

void DataDrivenNestedSurfaces::setUpAutomaticSurfaces(YieldSurfacePackage& yieldSurface, const double Gref, const double Pref) {
	// set up plastic modulus and hardening parameter sets based on the hyperbolic backbone model

	double  stress1(0.0), stress2(0.0), strain1(0.0), strain2(0.0), size(0.0),
		Hep_ref(0.0), Href(0.0), refStrain(0.0), peakShear(0.0), coneHeight(0.0),
		psi_bar(0.0), dilatancy(0.0);

	int tnys = yieldSurface.getTNYS();
	
	if (yieldSurface.isNonAssociated()) {
		psi_bar = tan(yieldSurface.getPsi() * M_PI / 180);	// compute the coefficient of dilation
	}

	if (yieldSurface.getPhi() > 0) {
		double sinPhi = sin(yieldSurface.getPhi() * M_PI / 180.);
		double MnyieldSurface = 6. * sinPhi / (3. - sinPhi);
		yieldSurface.setPresid(3. * yieldSurface.getCohesion() / (sqrt(2.) * MnyieldSurface));
		coneHeight = -(Pref - yieldSurface.getPresid());
		peakShear = sqrt(2.) * coneHeight * MnyieldSurface / 3.;
		refStrain = (yieldSurface.getPeakStrain() * peakShear) / (Gref * yieldSurface.getPeakStrain() - peakShear);
	}
	else if (yieldSurface.getPhi() == 0.0) {			// cohesion = peakShear
		peakShear = 2 * sqrt(2.) * yieldSurface.getCohesion() / 3;
		refStrain = (yieldSurface.getPeakStrain() * peakShear) / (Gref * yieldSurface.getPeakStrain() - peakShear);
		yieldSurface.setPresid(0.0);
	}

	double stressInc = peakShear / (tnys);

	for (int i = 0; i <= tnys; i++) {
		stress1 = i * stressInc;
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

		if (i == tnys) {
			Href = 0;
		}

		if (i == 0) {
			size = 3.0 / sqrt(2.0) * stressInc * 0.01;
		}

		yieldSurface.setTau(size, i);
		yieldSurface.setEta(Href, i);

		if (yieldSurface.isNonAssociated()) {
			// do contracting volumetric flow with hyperbolic decay [from psi_bar to 0]
			dilatancy = fabs(psi_bar) * 2 * (i - tnys) / (i - 2 * tnys);
			yieldSurface.setBeta(dilatancy, i);
		}
	}
}