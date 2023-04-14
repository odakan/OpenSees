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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/YieldSurfacePackage.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef YieldSurfacePackage_h
#define YieldSurfacePackage_h

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


class YieldSurfacePackage {
public:
	// null constructor
	YieldSurfacePackage(void) = default;

	// full constructors
	YieldSurfacePackage(int matID);													// online constructor
	YieldSurfacePackage(int matID, int tnys);										// offline constructor
	YieldSurfacePackage(const YieldSurfacePackage&) = default;							// default constructor
	YieldSurfacePackage(int matID, int tnys,											// auto-backbone constructor
		double cohesion, double frictionAngle, double dilatancyAngle, 
		double peakShearStrain, double residualPressure, double referencePressure);

	// destructor
	~YieldSurfacePackage(void);

	// operator overloading
	YieldSurfacePackage& operator = (const YieldSurfacePackage&) = default;				// one-to-one assignment

	// iteration control
	int commitState(void);
	int revertToLastCommit(void);

	// operational methods
	YieldSurfacePackage* getCopy(void);													// retun a copy of the object
	void printStats(bool detail);

	// get methods
	int getNYS(void);
	double getPhi(void);
	double getPsi(void);
	double getPresid(void);
	int getNYS_commit(void);
	double getCohesion(void);
	double getPeakStrain(void);
	double getTau(const int index, const int num_surface_commit);
	double getEta(const int index, const int num_surface_commit);
	Vector getAlpha(const int index, const int num_surface_commit);
	Vector getAlpha_commit(const int index, const int num_surface_commit);

	// set methods
	void incrementNYS(void);									// increment number of active yield surface (nYs)
	void setPhi(double value);
	void setPsi(double value);
	void setPresid(double value);
	void setCohesion(double value);
	void setTau(double value, const int index);
	void setEta(double value, const int index);
	void setAlpha(Vector value, const int index);
	void setAlpha_commit(Vector value, const int index);	


private:
	// yield surface paramters (local copy)
	double cohesion = 0.0;
	double frictionAngle = 0.0;
	double dilatancyAngle = 0.0;
	double peakShearStrain = 0.0;
	double residualPressure = 0.0;
	double referencePressure = 0.0;

	// operational parameters
	int matID = 0;
	int tnys = 0;
	bool do_online = false;

	// yield surface state variables
	Vector tau = Vector(1);					// limit isotropic stress
	Vector eta = Vector(1);					// plastic modulus
	Matrix alpha = Matrix(6, 1);			// backstress (TNYS + 1)
	Matrix alpha_commit = Matrix(6, 1);		// commited backstress (TNYS + 1)
	int nYs = 0;							// number of active yield surface
	int nYs_commit = 0;						// committed number of active yield surface

	// online update
	void update(void);
};
#endif