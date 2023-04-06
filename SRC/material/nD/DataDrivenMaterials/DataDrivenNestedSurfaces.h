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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/DataDrivenNestedSurfaces.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef DataDrivenNestedSurfaces_h
#define DataDrivenNestedSurfaces_h
#define _USE_MATH_DEFINES

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2022
//

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include "VonMisesDMM.h"
#include "DruckerPragerDMM.h"
#include "MatsuokaNakaiDMM.h"
#include "YieldSurfacePackage.h"


class DataDrivenNestedSurfaces {
public:
	// null constructor
	DataDrivenNestedSurfaces(void) = default;

	// full constructors
	DataDrivenNestedSurfaces(const DataDrivenNestedSurfaces&) = default;											// copy constructor
	DataDrivenNestedSurfaces(double cohesion, double frictionAngle, double peakShearStrain,							// von Mises type constructor
		double* HModuli, double* HParams);	
	DataDrivenNestedSurfaces(double cohesion, double frictionAngle, double dilationAngle, double peakShearStrain,	// Drucker-Prager type constructor
		double* HModuli, double* HParams);	
	DataDrivenNestedSurfaces(double cohesion, double frictionAngle, double dilationAngle, double peakShearStrain,	// Matsuoka-Nakai type constructor
		double* HModuli, double* HParams);	

	// destructor
	~DataDrivenNestedSurfaces(void);

	// operator overloading
	DataDrivenNestedSurfaces& operator=(const DataDrivenNestedSurfaces&) = default;		// one-to-one assignment

	// operational methods
	bool canDelete(void);													// return true if no other material is using the object
	void checkin(void);														// increase how_many counter
	void checkout(void);													// decrease how_many counter
	DataDrivenNestedSurfaces* getCopy(void);								// retun a copy of the object

	// update methods
	void updateTNYS(int var);
	void updateKref(double var);
	void updateGref(double var);
	void updatePref(double var);
	void updateModn(double var);
	void updatePhi(double var);
	void updatePsi(double var);
	void updateCohesion(double var);
	void updateHardParams(Vector& var);
	void updateHardParams(double var, int num);
	void updateDilatParams(Vector& var);
	void updateDilatParams(double var, int num);

	// get methods
	int getTNYS(void);
	double getKref(void);
	double getGref(void);
	double getPref(void);
	double getModn(void);
	double getPhi(void);
	double getPsi(void);
	double getCohesion(void);
	double getHref(int num);
	double getHP(int num);
	double getDP(int num);

	// setup yield surface
	YieldSurfacePackage* setUpYieldSurfaces(VonMisesDMM* theMaterial);
	YieldSurfacePackage* setUpYieldSurfaces(DruckerPragerDMM* theMaterial);
	YieldSurfacePackage* setUpYieldSurfaces(MatsuokaNakaiDMM* theMaterial);


private:
	double residualPressure = 0.0;
	Vector Href;								// plastic shear moduli
	Vector HardParams;							// hardening parameters
	Vector DilatParams;							// dilation parameters

	// surface generation options
	bool use_custom_surface = false;			// use user defined yield surface parameters

	// operational variables
	int how_many = 0;							// number of materials using this surface object

	// generate methods
	void generateYieldSurfaces(MultiYieldSurfaceHardeningSoftening* theMaterial);

};
#endif