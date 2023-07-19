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
// Created in:	April 2023
//
// Description: This file contains the implementation for the YieldSurfacePackage class.


#include "YieldSurfacePackage.h"


// Public methods
	// constructors
YieldSurfacePackage::YieldSurfacePackage(int mat) 
{
	// online constructor

	matID = mat;
	do_online = true;
	tnys = 3;
	nonassociated = true;

	// set material order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}

	// initililze only the  previous, current and next yield surfaces 
	tau = Vector(3);
	tau_commit = Vector(3);
	eta = Vector(3);
	eta_commit = Vector(3);
	beta = Vector(3);
	beta_commit = Vector(3);
	alpha = Matrix(nOrd, 3);
	alpha_commit = Matrix(nOrd, 3);
}

YieldSurfacePackage::YieldSurfacePackage(int mat, int t0)
{
	//offline constructor

	matID = mat;
	do_online = false;
	tnys = t0;
	nonassociated = true;

	// set material order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}

	// initialize all the yield surfaces
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	beta = Vector(tnys + 1);
	alpha = Matrix(nOrd, tnys + 1);
	alpha_commit = Matrix(nOrd, tnys + 1);

}


YieldSurfacePackage::YieldSurfacePackage(int mat, int t0, Vector hStrains, Vector hModuli, Vector hDilation)
{
	//custom-backbone constructor

	matID = mat;
	do_online = false;
	tnys = t0;

	// set material order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}

	// initialize all the yield surfaces
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	if (hDilation.Size() > 1) {
		nonassociated = true;
		beta = Vector(tnys + 1);
	}
	alpha = Matrix(nOrd, tnys + 1);
	alpha_commit = Matrix(nOrd, tnys + 1);
}


YieldSurfacePackage::YieldSurfacePackage(int mat, int t0,
	double c0, double f0, double d0, double p0, double Pres0, double Pref0)
{
	// auto-backbone constructor

	matID = mat;
	tnys = t0;
	cohesion = c0; frictionAngle = f0;
	dilatancyAngle = d0; peakShearStrain = p0; 
	residualPressure = Pres0; referencePressure = Pref0;
	do_online = false;

	// set material order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}


	// initialize all the yield surfaces
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	if (dilatancyAngle != 0.0) {
		nonassociated = true;
		beta = Vector(tnys + 1);
	}
	alpha = Matrix(nOrd, tnys + 1);
	alpha_commit = Matrix(nOrd, tnys + 1);
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
		opserr << "Attached material ID            =  " << matID << "\n";
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
		opserr << "Attached material ID   =  " << matID << "\n";
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
double YieldSurfacePackage::getPhi(void) { return frictionAngle; }
double YieldSurfacePackage::getPsi(void) { return dilatancyAngle; }
double YieldSurfacePackage::getPresid(void) { return residualPressure; }
int YieldSurfacePackage::getNYS_commit(void) { return nYs_commit; }
double YieldSurfacePackage::getCohesion(void) { return cohesion; }
double YieldSurfacePackage::getPeakStrain(void) { return peakShearStrain; }

double YieldSurfacePackage::getTau(const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getTau() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
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

	if (do_online) {
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

	if (do_online) {
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

Vector YieldSurfacePackage::getAlpha(const int index) {

	Vector Vect = Vector(nOrd);

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::getAlpha() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
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

	if (do_online) {
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
	
	if (do_online) {
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

void YieldSurfacePackage::setPhi(double value) { 
	
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
}

void YieldSurfacePackage::setTau(const double value, const int index) {

	if (index < 0) {
		opserr << "FATAL: YieldSurfacePackage::setTau() - a yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
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

	if (do_online) {
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

	if (do_online) {
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

	if (do_online) {

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

	if (do_online) {

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