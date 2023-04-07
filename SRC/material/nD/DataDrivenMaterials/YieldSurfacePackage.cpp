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

YieldSurfacePackage* YieldSurfacePackage::getCopy(void) 
{
	YieldSurfacePackage* copy = new YieldSurfacePackage(*this);
	return copy;
}

	// update methods
void YieldSurfacePackage::updateHardParams(Vector& var) { HardParams = var; }
void YieldSurfacePackage::updateHardParams(double var, int nYs_active) { HardParams(nYs_active) = var; }
void YieldSurfacePackage::updateDilatParams(Vector& var) { DilatParams = var; }
void YieldSurfacePackage::updateDilatParams(double var, int nYs_active) { DilatParams(nYs_active) = var; }

// get methods
double YieldSurfacePackage::getHref(int num) { return Href(num); }
double YieldSurfacePackage::getHP(int num) { return HardParams(num); }
double YieldSurfacePackage::getDP(int num) { return DilatParams(num); }

// Private methods