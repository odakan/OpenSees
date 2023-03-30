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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MYSHSmaterials/DataDrivenNestedSurfaces.h$
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

class DataDrivenNestedSurfaces {
public:
	// null constructor
	DataDrivenNestedSurfaces(void) = default;

	// full constructors
	DataDrivenNestedSurfaces(const DataDrivenNestedSurfaces&) = default;		// default constructor
		// traditional hyperbolic surface constructors
	DataDrivenNestedSurfaces(int tnys, double Kref, double Gref, double Pref,	// von Mises type constructor
		double modn, double cohesion, double peakShearStrain, double frictionAngle);
	DataDrivenNestedSurfaces(int TNYS, double Kref, double Gref, double Pref,	// Drucker-Prager type constructor
		double modn, double cohesion, double peakShearStrain, double frictionAngle, double dilationAngle);
	DataDrivenNestedSurfaces(int TNYS, double Kref, double Gref, double Pref,	// Matsuoka-Nakai type constructor
		double modn, double peakShearStrain, double frictionAngle);
		// data-driven surface constructors
	DataDrivenNestedSurfaces(int tnys, double Kref, double Gref, double Pref,	// von Mises type constructor
		double modn, double cohesion, double* Href, double* HardParams);
	DataDrivenNestedSurfaces(int TNYS, double Kref, double Gref, double Pref,	// Drucker-Prager type constructor
		double modn, double cohesion);
	DataDrivenNestedSurfaces(int TNYS, double Kref, double Gref, double Pref,	// Matsuoka-Nakai type constructor
		double modn);

	// destructor
	~DataDrivenNestedSurfaces(void);

	// operator overloading
	DataDrivenNestedSurfaces& operator=(const DataDrivenNestedSurfaces&) = default;	// one-to-one assignment

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

private:
	// yield surface input variables
	int TNYS = 0;								// total number of yield surfaces
	double Kref = 0.0;							// Reference bulk modulus
	double Gref = 0.0;							// Reference shear modulus
	double Pref = 0.0;							// Reference pressure
	double modn = 0.0;							// Modulus update power n
	double Phi = 0.0;							// internal friction angle (reference)
	double cohesion = 0.0;						// cohesion (reference)
	double peakStrain = 0.0;					// peak shear strain for the hyperbolic backbone
	double Psi = 0.0;							// dilation angle (reference)
	double residualPressure = 0.0;
	Vector Href;								// plastic shear moduli
	Vector HardParams;							// hardening parameters
	Vector DilatParams;							// dilation parameters

	// surface generation options
	bool use_custom_surface = false;			// use user defined yield surface parameters

	// operational variables
	int how_many = 0;							// number of materials using this surface object

	// generate methods
	void generateYieldSurfaces(void);

};
#endif