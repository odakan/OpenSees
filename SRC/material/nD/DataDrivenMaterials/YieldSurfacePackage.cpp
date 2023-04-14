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
	matID = mat;
	do_online = true;
	tnys = 3;
	// initililze only the  committed, current and next yield surfaces 
	tau = Vector(3);
	eta = Vector(3);
	alpha = Matrix(6, 3);
	alpha_commit = Matrix(6, 3);
}

YieldSurfacePackage::YieldSurfacePackage(int mat, int t0)
{
	matID = mat;
	do_online = false;
	tnys = t0;
	// initialize all the yield surfaces
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	alpha = Matrix(6, tnys + 1);
	alpha_commit = Matrix(6, tnys + 1);
}

YieldSurfacePackage::YieldSurfacePackage(int mat, int t0,
	double c0, double f0, double d0, double p0, double Pres0, double Pref0)
{
	matID = mat;
	tnys = t0;
	cohesion = c0; frictionAngle = f0;
	dilatancyAngle = d0; peakShearStrain = p0; 
	residualPressure = Pres0; referencePressure = Pref0;
	do_online = false;
	// initialize all the yield surfaces
	tau = Vector(tnys + 1);
	eta = Vector(tnys + 1);
	alpha = Matrix(6, tnys + 1);
	alpha_commit = Matrix(6, tnys + 1);
}


	// destructor
YieldSurfacePackage::~YieldSurfacePackage(void) 
{

}


	// iteration control
int YieldSurfacePackage::commitState(void) {

	alpha_commit = alpha;
	nYs_commit = nYs;

	return 0;
}

int YieldSurfacePackage::revertToLastCommit(void) {

	alpha = alpha_commit;
	nYs = nYs_commit;

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
		opserr << "Limit Stresses				=  " << tau;
		opserr << "Plastic Moduli				=  " << eta;
		opserr << "Active Y-Surface				=  " << nYs << "\n";
		opserr << "Commited Active Y-Surface	=  " << nYs_commit << "\n";
		opserr << "Back-stress					=  " << alpha;
		opserr << "Commited Back-stress			=  " << alpha_commit;
	}
	else {
		opserr << "YieldSurfacePackage::printStats() ->\n";
		opserr << "-------------------------------------------------------------------\n";
		opserr << "Active Y-Surface	=  " << nYs_commit << "\n";
		opserr << "Back-stress		=  " << alpha_commit;
	}
}


	// get methods
int YieldSurfacePackage::getNYS(void) { return nYs; }
double YieldSurfacePackage::getPhi(void) { return frictionAngle; }
double YieldSurfacePackage::getPsi(void) { return dilatancyAngle; }
double YieldSurfacePackage::getPresid(void) { return residualPressure; }
int YieldSurfacePackage::getNYS_commit(void) { return nYs_commit; }
double YieldSurfacePackage::getCohesion(void) { return cohesion; }
double YieldSurfacePackage::getPeakStrain(void) { return peakShearStrain; }

double YieldSurfacePackage::getTau(const int index, const int num_surface_commit) {

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::getTau: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
		if (index > num_surface_commit) {
			return tau(2);	// the next yield surface
		}
		else if (index == num_surface_commit) {
			return tau(1);	// the current yield surface
		}
		else if (index == (num_surface_commit - 1)) {
			return tau(0);	// the previous yield surface
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::getTau: A yield surface that is older than two steps was requested!";
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

double YieldSurfacePackage::getEta(const int index, const int num_surface_commit){

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
		if (index == (num_surface_commit + 1)) {
			return eta(2);	// the next yield surface
		}
		else if (index == num_surface_commit) {
			return eta(1);	// the current yield surface
		}
		else if (index == (num_surface_commit - 1)) {
			return eta(0);	// the previous yield surface
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface that is older than two steps was requested!";
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

Vector YieldSurfacePackage::getAlpha(const int index, const int num_surface_commit) {

	Vector Vect = Vector(6);

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
		if (index == (num_surface_commit + 1)) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha(i, 2);
			}
		}
		else if (index == num_surface_commit) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha(i, 1);
			}
		}
		else if (index == (num_surface_commit - 1)) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha(i, 0);
			}
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha(i, index);
			}
		}
		else {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha(i, tnys);
			}
		}
	}

	return Vect;
}

Vector YieldSurfacePackage::getAlpha_commit(const int index, const int num_surface_commit) {

	Vector Vect = Vector(6);

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {
		if (index == (num_surface_commit + 1)) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha_commit(i, 2);
			}
		}
		else if (index == num_surface_commit) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha_commit(i, 1);
			}
		}
		else if (index == (num_surface_commit - 1)) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha_commit(i, 0);
			}
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::getZeta: A yield surface that is older than two steps was requested!";
			exit(-1);
		}
	}
	else {
		if (index < tnys) {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha_commit(i, index);
			}
		}
		else {
			for (int i = 0; i < 6; i++) {
				Vect(i) = alpha_commit(i, tnys);
			}
		}
	}

	return Vect;
}


	// set methods
void YieldSurfacePackage::incrementNYS(void) {
	
	if (do_online) {
		nYs++;
	}
	else {
		if (nYs < tnys) {
			nYs++;
		}
	}
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
		opserr << "FATAL:YieldSurfacePackage::setTau: A yield surface with negative number was requested!";
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
			opserr << "FATAL:YieldSurfacePackage::setTau: A yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setEta(const double value, const int index) {

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::setZeta: A yield surface with negative number was requested!";
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
			opserr << "FATAL:YieldSurfacePackage::setZeta: A yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setAlpha(const Vector value, const int index) {

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::setAlpha: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {

	}
	else {
		if (index <= tnys) {
			for (int i = 0; i < 6; i++) {
				alpha(i, index) = value(i);
			}
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::setAlpha: A yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}

void YieldSurfacePackage::setAlpha_commit(const Vector value, const int index) {

	if (index < 0) {
		opserr << "FATAL:YieldSurfacePackage::setAlpha: A yield surface with negative number was requested!";
		exit(-1);
	}

	if (do_online) {

	}
	else {
		if (index <= tnys) {
			for (int i = 0; i < 6; i++) {
				alpha_commit(i, index) = value(i);
			}
		}
		else {
			opserr << "FATAL:YieldSurfacePackage::setAlpha: A yield surface further than TNYS was requested!";
			exit(-1);
		}
	}
}


// Private methods
void YieldSurfacePackage::update(void) {

}