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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MYSHSmaterials/DataDrivenNestedSurfaces.cpp$
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

// Public methods
	// full constructors
DataDrivenNestedSurfaces::DataDrivenNestedSurfaces(int tnys, double k, double g, double p,			// von Mises constructor
	double n, double c, double* gp, double* hp)
{
	// initialize variables
	TNYS = tnys; Kref = k; Gref = g;
	Pref = p; modn = n; cohesion = c;

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
void DataDrivenNestedSurfaces::updateTNYS(int var) { TNYS = var; }
void DataDrivenNestedSurfaces::updateKref(double var) { Kref = var; }
void DataDrivenNestedSurfaces::updateGref(double var) { Gref = var; }
void DataDrivenNestedSurfaces::updatePref(double var) { Pref = var; }
void DataDrivenNestedSurfaces::updateModn(double var) { modn = var; }
void DataDrivenNestedSurfaces::updatePhi(double var) { Phi = var; }
void DataDrivenNestedSurfaces::updatePsi(double var) { Psi = var; }
void DataDrivenNestedSurfaces::updateCohesion(double var) { cohesion = var; }
void DataDrivenNestedSurfaces::updateHardParams(Vector& var) { HardParams = var; }
void DataDrivenNestedSurfaces::updateHardParams(double var, int nYs_active) { HardParams(nYs_active) = var; }
void DataDrivenNestedSurfaces::updateDilatParams(Vector& var) { DilatParams = var; }
void DataDrivenNestedSurfaces::updateDilatParams(double var, int nYs_active) { DilatParams(nYs_active) = var; }

// get methods
int DataDrivenNestedSurfaces::getTNYS(void) { return TNYS; }
double DataDrivenNestedSurfaces::getKref(void) { return Kref; }
double DataDrivenNestedSurfaces::getGref(void) { return Gref; }
double DataDrivenNestedSurfaces::getPref(void) { return Pref; }
double DataDrivenNestedSurfaces::getModn(void) { return modn; }
double DataDrivenNestedSurfaces::getPhi(void) { return Phi; }
double DataDrivenNestedSurfaces::getPsi(void) { return Psi; }
double DataDrivenNestedSurfaces::getCohesion(void) { return cohesion; }
double DataDrivenNestedSurfaces::getHref(int num) { return Href(num); }
double DataDrivenNestedSurfaces::getHP(int num) { return HardParams(num); }
double DataDrivenNestedSurfaces::getDP(int num) { return DilatParams(num); }

	// generate methods
void DataDrivenNestedSurfaces::generateYieldSurfaces(void)
{

}

// Private methods