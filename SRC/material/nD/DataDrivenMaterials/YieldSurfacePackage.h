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
// $Date: 2023-XX-XX XX:XX:XX $

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
#include <elementAPI.h>
#include "DataDrivenNestedSurfaces.h"

class YieldSurfacePackage {
public:
	// null constructor
	YieldSurfacePackage(void) = default;

	// full constructors	
	YieldSurfacePackage(const int matID, const int driver, DataDrivenNestedSurfaces* ptr, 
		const double Gref, const double Pref, const Vector& stress, const Vector& strain, const bool verbosity);

	// destructor
	~YieldSurfacePackage(void);

	// operator overloading
	YieldSurfacePackage& operator = (const YieldSurfacePackage&) = default;				// one-to-one assignment

	// iteration control
	int commitState(void);
	int revertToLastCommit(void);

	// operational methods
	YieldSurfacePackage* getCopy(void);		// retun a copy of the object
	void printStats(bool detail);
	bool isNonAssociated(void);				// return true if the flow rule non-associated

	// get methods
	int now(void);
	int next(void);
	int prev(void);
	int getTag(void);
	int getNYS(void);
	int getTNYS(void);
	int getNYS_commit(void);
	double getTau(const int index);
	double getEta(const int index);
	double getBeta(const int index);
	Vector getAlpha(const int index);
	Vector getAlpha_commit(const int index);

	// set methods
	void increment(void);									// increment current yield surface (ys)
	void setTNYS(int value);
	void setTau(double value, const int index);
	void setEta(double value, const int index);
	void setBeta(double value, const int index);
	void setAlpha(Vector value, const int index);
	void setAlpha_commit(Vector value, const int index);	


private:
	// operational parameters
	int nOrd = 6;								// material order
	int matID = 0;								// tag of the attached material
	int tnys = 0;								// total number of yield surfaces
	int datadriver = 0;							// yield surface update method
	bool do_active = false;						// activate on-the-fly update of the data-driven surfaces
	bool beVerbose = false;						// be verbose about internal processes (use for debugging) (no by default)
	bool nonassociated = false;					// associated or nonassociated flow
	DataDrivenNestedSurfaces* lib = nullptr;	// pointer to the yield surface library

	// yield surface state variables
	/*	Note: nYs != ys
			nYs indicates how many surfaces needs to be translated between iterations
			whereas ys stores the index of the current yield surface.
			
		IMPLEX variables:
			ys and ys_commit are used during the explicit stage of the IMPLEX iterations.
			The extrapolation is done by keeping ys constant -> ys = ys_commit.
	*/
	Vector tau = Vector(1);					// limit isotropic stress
	Vector tau_commit = Vector(1);			// committed limit isotropic stress (used in active update)
	Vector eta = Vector(1);					// plastic modulus
	Vector eta_commit = Vector(1);			// committed plastic modulus (used in active update)
	Vector beta = Vector(1);				// dilatancy parameter
	Vector beta_commit = Vector(1);		    // committed dilatancy parameter (used in active update)
	Vector gamma = Vector(1);				// hardening parameter (for output) (used in active and passive update)
	Vector gamma_commit = Vector(1);		// committed hardening parameter (for output) (used in active update)
	Matrix alpha = Matrix(6, 1);			// backstress (TNYS + 1)
	Matrix alpha_commit = Matrix(6, 1);		// commited backstress (TNYS + 1)
	int nYs = 0;							// number of active yield surfaces
	int nYs_commit = 0;						// committed number of active yield surfaces
	int num = 0;							// current yield surface (!= nYs if online)
	int num_commit = 0;						// committed yield surface

	// online update
	void update(void);
};
#endif