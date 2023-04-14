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
// Created in:	April 2023
//

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include "YieldSurfacePackage.h"

class DataDrivenNestedSurfaces {
public:
	// null constructor
	DataDrivenNestedSurfaces(void) = default;

	// full constructors
	DataDrivenNestedSurfaces(const DataDrivenNestedSurfaces&) = default;								// copy constructor
	DataDrivenNestedSurfaces(int tag, double cohesion, double frictionAngle, double dilationAngle,		// full constructor
		double peakShearStrain,	double tnys, double* HModuli, double* HParams, bool verbosity);

	// destructor
	~DataDrivenNestedSurfaces(void);

	// operator overloading
	DataDrivenNestedSurfaces& operator=(const DataDrivenNestedSurfaces&) = default;		// one-to-one assignment

	// operational methods
	bool canDelete(void);						// return true if no material is using the object
	void checkin(void);							// increase how_many counter
	void checkout(void);						// decrease how_many counter
	bool isAOK(int dataDriver);					// whether the desired yield surface package is allowed or not 
	DataDrivenNestedSurfaces* getCopy(void);	// retun a copy of the object

	// get methods
	int getTNYS(void);
	double getCohesion(void);
	double getPhi(void);
	double getPsi(void);
	double getPeakStrain(void);
	double getPref(void);

	// generate yield surface
	YieldSurfacePackage generateYieldSurfaces(const int matid, const int dataDriver, const double Pref, const double Gref, const double TNYS);

private:
	// default yield surface paramters
	int tnys_init = 0;
	double cohesion_init = 0.0;
	double frictionAngle_init = 0.0;
	double dilatancyAngle_init = 0.0;
	double peakShearStrain_init = 0.0;
	double residualPressure_init = 0.0;
	double referencePressure_init = 0.0;

	// operational variables
	bool isAutomaticOK = false;
	bool isOfflineOK = false;
	bool isOnlineOK = false;
	bool beVerbose = false;		// be verbose about internal processes (use for debugging) (no by default)					
	int matID = 0;				// tag of the attached inital material object
	int how_many = 0;			// number of material sub-objects using this nested surface super-object

	// representative volume element data
	Vector HModuli;
	Vector HParams;

private:
	// set up yield surfaces
	void setUpOnlineSurfaces(YieldSurfacePackage& yieldSurface);
	void setUpOfflineSurfaces(YieldSurfacePackage& yieldSurface);
	void setUpAutomaticSurfaces(YieldSurfacePackage& yieldSurface, const double Pref, const double Gref, const int tnys);

};
#endif