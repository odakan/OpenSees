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


// Public methods
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


}


// Private methods