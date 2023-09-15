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
// $Date: 2023-XX-XX XX:XX:XX $

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
	double psi0, double s0, int nys, double deps, Vector hmoduli, Vector hparams, 
	Vector dparams, const char* dbPathDir, const char* dbMainFile, bool verbosity):
	db(matID, dbPathDir, dbMainFile), matID(tag), beVerbose(verbosity), strain_discretize(deps),
	cohesion(c0), frictionAngle(phi0), dilatancyAngle(psi0), peakShearStrain(s0), tnys(nys)
{
	// check data-driven yield surface generation 
	if (dbPathDir != nullptr && dbMainFile != nullptr) {
		// if database is present
		isPassiveOK = true;
		isActiveOK = true;
	}
	else {
		// if database is not present
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - no material response database is provided. Data-driven yield surfaces are NOT allowed...\n"; }
		isPassiveOK = false;
		isActiveOK = false;
	}

		// check automatic, user-custom and offline surface generation
	if (tnys <= 0) {
		// if total number of yield surfaces is less than one, prohibit automatic and offline surfaces
		if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - TNYS <= 0: AUTOMATIC and PASSIVE data-driven yield surfaces are NOT allowed...\n"; }
		isUserCustomOK = false;
		isAutomaticOK = false;
		isPassiveOK = false;
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
		if (peakShearStrain <= 0) {
			// if peak shear strain is zero or less than zero, prohibit automatic surfaces
			if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - peakShearStrain <= 0: AUTOMATIC yield surfaces are NOT allowed...\n"; }
			isAutomaticOK = false;;
		}
		else {
			if (frictionAngle > 0) {
				isAutomaticOK = true;
			}
			else if (cohesion > 0) {
				isAutomaticOK = true;
			}
			else {
				if (beVerbose) { opserr << "WARNING: DataDrivenNestedSurfaces() - either cohesion or friction angle must be > 0: AUTOMATIC yield surfaces are NOT allowed...\n"; }
				isAutomaticOK = false;
			}
		}
	}
	
	if (isPassiveOK || isActiveOK) {
		if (beVerbose) opserr << db;
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
	else if (dataDriver == 1 && isPassiveOK) {
		aok = true;
	}
	else if (dataDriver == 2 && isActiveOK) {
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
int DataDrivenNestedSurfaces::getTNYS(void) { return tnys; }
double DataDrivenNestedSurfaces::getCohesion(void) { return cohesion; }
double DataDrivenNestedSurfaces::getFrictionAngle(void) { return frictionAngle; }
double DataDrivenNestedSurfaces::getDilatancyAngle(void) { return dilatancyAngle; }
double DataDrivenNestedSurfaces::getPeakStrain(void) { return peakShearStrain; }
double DataDrivenNestedSurfaces::getPref(void) { return referencePressure; }


// set methods
void DataDrivenNestedSurfaces::setTNYS(const int val) { tnys = val; }
void DataDrivenNestedSurfaces::setCohesion(const double val) { cohesion = val; }
void DataDrivenNestedSurfaces::setFrictionAngle(const double val) { frictionAngle = val; }
void DataDrivenNestedSurfaces::setDilatancyAngle(const double val) { dilatancyAngle = val; }
void DataDrivenNestedSurfaces::setPeakStrain(const double val) { peakShearStrain = val; }
void DataDrivenNestedSurfaces::setPref(const double val) { referencePressure = val; }

// set-up yield surfaces
void DataDrivenNestedSurfaces::setUpActiveSurfaces(int& nys, Vector& tau, Vector& eta, Vector& beta,
	Vector& gamma, const Vector& stress, const Vector& strain) {
	// set-up plastic modulus, hardening parameter and dilatancy parameter based on the database, on-the-fly
	
	// initialize vectors
	nys = 2;
	tau = Vector(nys);
	gamma = Vector(nys);	// hardening parameters
	eta = Vector(nys);
	beta = Vector(nys);

	// initialize local variables
	int dim = db.getDim();
	double p_avg = 0.0;
	double e_vol = 0.0;
	Vector e_dev(int(3 * (dim - 1)));
	Vector s_dev(int(3 * (dim - 1)));
	if (dim == 2) {
		p_avg = (stress(0) + stress(1)) * 0.5;
		s_dev = stress; s_dev(0) -= p_avg; s_dev(1) -= p_avg;
		e_vol = strain(0) + strain(1);
		e_dev = strain; e_dev(0) -= e_vol * 0.5; e_dev(1) -= e_vol * 0.5;
	}
	else {
		p_avg = (stress(0) + stress(1) + stress(2)) * 1.0 / 3.0;
		s_dev = stress; s_dev(0) -= p_avg; s_dev(1) -= p_avg; s_dev(2) -= p_avg;
		e_vol = strain(0) + strain(1) + strain(2);
		e_dev = strain; e_dev(0) -= e_vol * 1.0 / 3.0; e_dev(1) -= e_vol * 1.0 / 3.0; e_dev(2) -= e_vol * 1.0 / 3.0;
	}
	
	// initial yield surface (elasticity limit)
	tau(0) = sqrt(1.0 / 3.0 * TensorM::dotdot(s_dev, s_dev));
	gamma(0) = 2 * sqrt(1.0 / 3.0 * TensorM::dotdot(e_dev, e_dev));
	
	// scan the database
	int index = db.seek(p_avg, "strain", strain);
	opserr << "DataDrivenNestedSurfaces::setUpActiveSurfaces() - works until here!\n";
	// next yield surface
	tau(1) = db.getOctahedralStress(index);
	gamma(1) = db.getOctahedralStrain(index);

	// do zero value checks
	if (tau(1) < ZERO_VALUE) { tau(1) = 1; }			// next yield radius
	if (tau(0) < ZERO_VALUE) { tau(0) = 0.1 * tau(1); }	// current yield radius
	if (gamma(1) < ZERO_VALUE) { gamma(1) = 1; }		// next hardening parameter

	// compute moduli [index 0: old (elastic), 1: current]
	eta(1) = (tau(1) - tau(0)) / (gamma(1) - gamma(0));		// kinematic modulus
	// convert octahedral shear strain to the magnitude of the strain deviator tensor
	double sqEE = sqrt(3.0 / 4.0 * gamma(1) * gamma(1));
	beta(1) = (db.getVolumetricStrain(index) - e_vol) / (sqEE - sqrt(TensorM::dotdot(e_dev, e_dev)));	// dilatancy parameter
}

void DataDrivenNestedSurfaces::setUpPassiveSurfaces(int& nys, Vector& tau, Vector& eta, Vector& beta,
	Vector& gamma, const Vector& stress, const Vector& strain) {
	// set-up plastic modulus, hardening parameter and dilatancy parameter sets based on the database, at the begining

	// initialize vectors
	nys = tnys;
	tau = Vector(tnys +1);
	eta = Vector(tnys + 1);
	beta = Vector(tnys + 1);
	gamma = Vector(tnys + 1);

	// 

}

void DataDrivenNestedSurfaces::setUpUserCustomSurfaces(int& nys, bool& nonassociated, Vector& tau, Vector& eta, Vector& beta,
	const double Gref, const double Pref) {
	// set-up plastic modulus, hardening parameter and dilatancy parameter sets based on the user input strain and G/Gmax pairs

	// initialize vectors
	nys = tnys;
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	if (isNonassociatedOK) {
		nonassociated = true;
		beta = Vector(tnys + 1);
	}

	// initialize local variables
	double  stress1(0.0), stress2(0.0), strain1(0.0),   strain2(0.0),   size(0.0),       dilatancy(0.0),
		    Hep_ref(0.0), Href(0.0),    refStrain(0.0), peakShear(0.0), coneHeight(0.0), strain_vol1(0.0), strain_vol2(0.0);

	// do some check and regulations
	if (frictionAngle > 0) {   // ignore user defined friction angle
		double tmax = 0.0; for (int i = 0; i < tnys; i++) { tmax = fmax(tmax, Gref * HModuli(i) * HParams(i)); }
		double Mnys = -(sqrt(3.0) * tmax - 2. * cohesion) / Pref;
		if (Mnys <= 0) {   // also ignore user defined cohesion
			cohesion = sqrt(3.0) / 2.0 * tmax;
			frictionAngle = 0.0;
			coneHeight = 1.;
			residualPressure = 0.0;
		}
		else {
			double sinPhi = 3 * Mnys / (6 + Mnys);
			if (sinPhi < 0. || sinPhi>1.) {
				opserr << "FATAL: DataDrivenNestedSurfaces::setUpUserCustomSurfaces() - invalid friction angle for nDmaterial "
					<< matID << ", please modify reference pressure or G/Gmax curve.\n";
				exit(-1);
			}
			residualPressure = 2. * cohesion / Mnys;
			if (residualPressure < 0.01 * Pref) residualPressure = 0.01 * Pref;
			coneHeight = -(Pref - residualPressure);
			frictionAngle = asin(sinPhi) * 180 / M_PI;
		}
	}
	else {   // ignore user defined cohesion
		double tmax = 0.0; for (int i = 0; i < tnys; i++) { tmax = fmax(tmax, Gref * HModuli(i) * HParams(i)); }
		cohesion = sqrt(3.) / 2 * tmax;
		coneHeight = 1.0;
		residualPressure = 0.0;
	}

	/*
	 *	Required user-input parameters for each yield surface
	 *	-----------------------------------------------------------------------------------------------------------------------------------
	 *	HParams: -> Octahedral shear strain = 2.0 * sqrt(2.0 / 3.0 * J2') where J2' is the second invariant of the strain deviator tensor 
	 *	HModuli: -> Secant shear modulus reduction curve (Gsec / Gmax)
	 *	DParams: -> Secant dilatancy curve (Evol / Eoct) where Evol and Eoct are the volumetric and octahedral shear strains, respectively
	*/

	// first yield surface
	if (frictionAngle > 0) {
		size = Gref * HModuli(0) * HParams(0) / coneHeight;
	}
	else if (frictionAngle == 0) {
		size =  Gref * HModuli(0) * HParams(0);
	}
	Href = Gref * HModuli(0) * HParams(0) / HParams(0);
	if (Href > LARGE_NUMBER) Href = LARGE_NUMBER;
	tau(0) = size * 0.01;
	eta(0) = Href;
	if (nonassociated) {
		beta(0) = DParams(0);
	}

	// last yield surface
	if (frictionAngle > 0) {
		size = Gref * HModuli(tnys - 1) * HParams(tnys - 1) / coneHeight;
	}
	else if (frictionAngle == 0) {
		size = Gref * HModuli(tnys - 1) * HParams(tnys - 1);
	}
	Href = 0;
	tau(tnys) = size;
	eta(tnys) = Href;
	if (nonassociated) {
		beta(tnys) = 0; // assume constant volume deformation
	}

	// other yield surfaces
	for (int i = 1; i < (tnys); i++) {
		strain1 = HParams(i-1);
		stress1 = Gref * HModuli(i - 1) * HParams(i - 1);
		strain2 = HParams(i);
		stress2 = Gref * HModuli(i) * HParams(i);
		if (frictionAngle > 0.) {
			size = stress1 / coneHeight;
		}
		else if (frictionAngle == 0.) {
			size = stress1;
		}
		Href = (stress2 - stress1) / (strain2 - strain1);
		if (Href > LARGE_NUMBER) Href = LARGE_NUMBER;
		tau(i) = size;
		eta(i) = Href;
		if (nonassociated) {
			strain_vol1 = DParams(i - 1) * strain1;
			strain_vol2 = DParams(i) * strain2;
			dilatancy = (strain_vol2 - strain_vol1) / (strain2 - strain1);
			beta(i) = dilatancy;
		}
	}
}

void DataDrivenNestedSurfaces::setUpAutomaticSurfaces(int& nys, bool& nonassociated, Vector& tau, Vector& eta, Vector& beta,
	const double Gref, const double Pref) {
	// set-up plastic modulus, hardening parameter and dilatancy parameter sets based on the hyperbolic backbone model

	// initialize vectors
	nys = tnys;
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	if (dilatancyAngle != 0.0) {
		nonassociated = true;
		beta = Vector(tnys + 1);
	}

	// initialize local variables
	double  stress1(0.0), stress2(0.0), strain1(0.0), strain2(0.0), size(0.0),
		Hep_ref(0.0), Href(0.0), refStrain(0.0), peakShear(0.0), coneHeight(0.0),
		psi_bar(0.0), dilatancy(0.0), Presid(0.0);

	if (frictionAngle > 0) {
		double sinPhi = sin(frictionAngle * M_PI / 180.);
		double MnyieldSurface = 6. * sinPhi / (3. - sinPhi);
		Presid = 3. *cohesion / (sqrt(2.) * MnyieldSurface);
		coneHeight = -(Pref - Presid);
		peakShear = sqrt(2.) * coneHeight * MnyieldSurface / 3.;
		refStrain = (peakShearStrain * peakShear) / (Gref * peakShearStrain - peakShear);
	}
	else if (frictionAngle == 0.0) {			// cohesion = 2 * sqrt(2.) / 3.0 * peakShear
		peakShear = 2.0 * sqrt(2.) *cohesion / 3.0;
		refStrain = (peakShearStrain * peakShear) / (Gref * peakShearStrain - peakShear);
		Presid = 0.0;
	}

	double stressInc = peakShear / tnys;

	// first yield surface
	stress1 = stressInc * 0.01;
	strain1 = stress1 * refStrain / (Gref * refStrain - stress1);
	stress2 = stressInc;
	strain2 = stress2 * refStrain / (Gref * refStrain - stress2);
	if (frictionAngle > 0.) {
		size = 3. * stress1 / sqrt(2.) / coneHeight;
	}
	else if (frictionAngle == 0.) {
		size = 3. * stress1 / sqrt(2.);
	}
	Href = (stress2 - stress1) / (strain2 - strain1);
	tau(0) = size;
	eta(0) = Href;

	if (nonassociated) {
		psi_bar = tan(dilatancyAngle * M_PI / 180);	// compute the coefficient of dilation
		beta(0) = psi_bar;
	}

	for (int i = 1; i <= tnys; i++) {
		stress1 = stressInc * i;
		strain1 = stress1 * refStrain / (Gref * refStrain - stress1);
		stress2 = stress1 + stressInc;
		strain2 = stress2 * refStrain / (Gref * refStrain - stress2);
		if (frictionAngle > 0.) {
			size = 3. * stress1 / sqrt(2.) / coneHeight;
		}
		else if (frictionAngle == 0.) {
			size = 3. * stress1 / sqrt(2.);
		}

		Href = (stress2 - stress1) / (strain2 - strain1);

		if (Href > LARGE_NUMBER) { Href = LARGE_NUMBER; }
		if (i == tnys) { Href = 0; }

		tau(i) = size;
		eta(i) = Href;

		if (nonassociated) {
			// do contracting volumetric flow with hyperbolic decay [from psi_bar to 0]
			dilatancy = fabs(psi_bar) * 2 * (i - tnys) / (i - 2 * tnys);
			beta(i) = dilatancy;
		}
	}
}


// Private methods