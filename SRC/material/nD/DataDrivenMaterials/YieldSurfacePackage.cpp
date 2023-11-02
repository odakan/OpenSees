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
// $Date: 2023-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	April 2023
//
// Description: This file contains the implementation for the YieldSurfacePackage class.

#include "YieldSurfacePackage.h"

// Public methods
	// constructors
YieldSurfacePackage::YieldSurfacePackage(const int mat, const int sub):
	matID(mat), subID(sub)
{
}

YieldSurfacePackage::YieldSurfacePackage(const int mat, const int sub, const int ord, const int driver,
	std::shared_ptr<DataDrivenNestedSurfaces> ptr, const double Gref, const double Pref, const Vector& stress, 
	const Vector& strain, const bool verbosity):
	matID(mat), subID(sub), nOrd(ord), datadriver(driver), library(ptr), beVerbose(verbosity)
{
	// lock the weak_ptr to create a shared_ptr
	if (auto lib = library.lock()) {
		// generate yield surfaces
		if (datadriver == 0) {							// automatic backbone genration
			// auto-backbone constructor
			do_active = false;

			// setup automatic yield surfaces
			lib->setUpAutomaticSurfaces(tnys, nonassociated, tau, eta, beta, Gref, Pref);

			// initialize remaining variables
			alpha = Matrix(nOrd, tnys + 1);
			alpha_commit = Matrix(nOrd, tnys + 1);
		}
		else if (datadriver == 1) {						// passive: do not update once a set once generated
			//passive constructor
			do_active = false;
			nonassociated = true;

			// setup passive yield surfaces
			lib->setUpPassiveSurfaces(tnys, tau, eta, beta, gamma, stress, strain);

			// initialize remaining variables
			alpha = Matrix(nOrd, tnys + 1);
			alpha_commit = Matrix(nOrd, tnys + 1);
		}
		else if (datadriver < -1 || datadriver > 1) {	// active: generate surfaces on-the-fly
			// active constructor
			do_active = true;
			nonassociated = true;

			// setup active yield surfaces
			lib->setUpActiveSurfaces(tnys, tau, eta, beta, gamma, stress, strain);
		
			// initialize remaining variables
			tau_commit = Vector(tau.Size());
			gamma_commit = Vector(gamma.Size());
			eta_commit = Vector(eta.Size());
			if (nonassociated) {
				beta_commit = Vector(beta.Size());
			}
			alpha = Matrix(nOrd, tnys);
			alpha_commit = Matrix(nOrd, tnys);
		}
		else if (datadriver == -1) {					// user custom surface generation
			//custom-backbone constructor
			do_active = false;

			// setup custom yield surfaces
			lib->setUpUserCustomSurfaces(tnys, nonassociated, tau, eta, beta, Gref, Pref);

			// initialize remaining variables
			alpha = Matrix(nOrd, tnys + 1);
			alpha_commit = Matrix(nOrd, tnys + 1);
		} 
		else {
			opserr << "FATAL: YieldSurfacePackage() - unknown yield surface update method!\n";
			exit(-1);
		}
	}
	else {
		opserr << "FATAL: YieldSurfacePackage() - cannot access the database! Database is no longer valid...\n";
		exit(-1);
	}

	// print some detalied statistics
	if (beVerbose) { opserr << "YieldSurfacePackage() -> printStats() ->\n"; printStats(true); }
}

	// destructor
YieldSurfacePackage::~YieldSurfacePackage(void) 
{

}


	// iteration control
int YieldSurfacePackage::commitState(void) {

	alpha_commit = alpha;
	nYs_commit = nYs;
	num_commit = num;

	return 0;
}

int YieldSurfacePackage::revertToLastCommit(void) {

	alpha = alpha_commit;
	nYs = nYs_commit;
	num = num_commit;

	return 0;
}


	// operational methods
YieldSurfacePackage* YieldSurfacePackage::getCopy(void) {

	YieldSurfacePackage* copy = new YieldSurfacePackage(*this);
	return copy;
}

void YieldSurfacePackage::printStats(bool detail) {

	if (detail) {
		opserr << "YieldSurfacePackage::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Attached material ID            =  " << matID << "::" << subID << "\n";
		opserr << "Limit Stresses                  =  " << tau;
		opserr << "Commited limit stresses         =  " << tau_commit;
		opserr << "Plastic Moduli                  =  " << eta;
		opserr << "Commited plastic moduli         =  " << eta_commit;
		opserr << "Dilatancy parameters            =  " << beta;
		opserr << "Commited dilatancy parameters   =  " << beta_commit;
		opserr << "Active no. of surfaces          =  " << nYs << "\n";
		opserr << "Commited active no. of surfaces =  " << nYs_commit << "\n";
		opserr << "Current yield surface           =  " << num << "\n";
		opserr << "Commited current yield surface  =  " << num_commit << "\n";
		opserr << "Back-stress                     =  " << alpha;
		opserr << "Commited Back-stress            =  " << alpha_commit;
	}
	else {
		opserr << "YieldSurfacePackage::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Current yield surface  =  " << num << "\n";
		opserr << "Active no. of surfaces =  " << nYs_commit << "\n";
		opserr << "Back-stress            =  " << alpha_commit;
	}
}

bool YieldSurfacePackage::isNonAssociated(void) { return nonassociated; }

	// get methods
int YieldSurfacePackage::now(void) { return num; }
int YieldSurfacePackage::next(void) { return num + 1; }
int YieldSurfacePackage::prev(void) { return num - 1; }
int YieldSurfacePackage::getTag(void) { return matID; }
int YieldSurfacePackage::getNYS(void) { return nYs; }
int YieldSurfacePackage::getTNYS(void) { return tnys; }
int YieldSurfacePackage::getNYS_commit(void) { return nYs_commit; }

double YieldSurfacePackage::getTau(const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getTau() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		if (index > nYs_commit) {
			return tau(2);	// the next yield surface
		}
		else if (index == nYs_commit) {
			return tau(1);	// the current yield surface
		}
		else if (index == (nYs_commit - 1)) {
			return tau(0);	// the previous yield surface
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::getTau() - a yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			return tau(index);
		}
		else {
			return tau(tnys);
		}
	}
}

double YieldSurfacePackage::getEta(const int index){

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getEta() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		if (index == (nYs_commit + 1)) {
			return eta(2);	// the next yield surface
		}
		else if (index == nYs_commit) {
			return eta(1);	// the current yield surface
		}
		else if (index == (nYs_commit - 1)) {
			return eta(0);	// the previous yield surface
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::getEta() - a yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			return eta(index);
		}
		else {
			return eta(tnys);
		}
	}
}

double YieldSurfacePackage::getBeta(const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getBeta() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		if (index == (nYs_commit + 1)) {
			return beta(2);	// the next yield surface
		}
		else if (index == nYs_commit) {
			return beta(1);	// the current yield surface
		}
		else if (index == (nYs_commit - 1)) {
			return beta(0);	// the previous yield surface
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::getBeta() - a yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			return beta(index);
		}
		else {
			return beta(tnys);
		}
	}
}

double YieldSurfacePackage::getAttraction(const int index) { return cohesion / getTau(index); }

Vector YieldSurfacePackage::getAlpha(const int index) {

	Vector Vect = Vector(nOrd);

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getAlpha() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		if (index == (nYs_commit + 1)) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha(i, 2);
			}
		}
		else if (index == nYs_commit) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha(i, 1);
			}
		}
		else if (index == (nYs_commit - 1)) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha(i, 0);
			}
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::getAlpha() - a yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha(i, index);
			}
		}
		else {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha(i, tnys);
			}
		}
	}

	return Vect;
}

Vector YieldSurfacePackage::getAlpha_commit(const int index) {

	Vector Vect = Vector(nOrd);

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getAlpha_commit() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		if (index == (nYs_commit + 1)) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha_commit(i, 2);
			}
		}
		else if (index == nYs_commit) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha_commit(i, 1);
			}
		}
		else if (index == (nYs_commit - 1)) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha_commit(i, 0);
			}
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::getAlpha_commit() - a yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha_commit(i, index);
			}
		}
		else {
			for (int i = 0; i < nOrd; i++) {
				Vect(i) = alpha_commit(i, tnys);
			}
		}
	}

	return Vect;
}


	// set methods
void YieldSurfacePackage::increment(void) {
	
	if (do_active) {
		num++;
		// call for new surface
		nYs = tau.Size();
	}
	else {
		if (num < tnys) {
			num++;
			nYs = num;
		}
	}
}

void YieldSurfacePackage::setTNYS(int value) {

	tnys = value;
	// ask for re-compuation of yield surfaces with new value
}

/*void YieldSurfacePackage::setPhi(double value) {
	
	frictionAngle = value;
	// ask for re-compuation of yield surfaces with new value
}

void YieldSurfacePackage::setPsi(double value) {
	
	dilatancyAngle = value;
	// ask for re-compuation of yield surfaces with new value
}

void YieldSurfacePackage::setPresid(double value) {

	residualPressure = value;
	// ask for re-compuation of yield surfaces with new value
}

void YieldSurfacePackage::setCohesion(double value) {

	cohesion = value;
	// ask for re-compuation of yield surfaces with new value
}*/

void YieldSurfacePackage::setTau(const double value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setTau() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		tau(index) = value;
	}
	else {
		if (index <= tnys) {
			tau(index) = value;
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::setTau() - a yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setEta(const double value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setEta() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		eta(index) = value;
	}
	else {
		if (index <= tnys) {
			eta(index) = value;
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::setEta() - a yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setBeta(const double value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setBeta() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {
		beta(index) = value;
	}
	else {
		if (index <= tnys) {
			beta(index) = value;
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::setBeta() - a yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setAlpha(const Vector value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setAlpha() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {

	}
	else {
		if (index <= tnys) {
			for (int i = 0; i < nOrd; i++) {
				alpha(i, index) = value(i);
			}
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::setAlpha() - a yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setAlpha_commit(const Vector value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setAlpha_commit() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_active) {

	}
	else {
		if (index <= tnys) {
			for (int i = 0; i < nOrd; i++) {
				alpha_commit(i, index) = value(i);
			}
		}
		else {
			opserr << "FATAL: YieldSurfacePackage::setAlpha_commit() - a yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}


// Private methods
void YieldSurfacePackage::update(void) {

}