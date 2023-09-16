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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialDataPoint.h$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

#ifndef MaterialDataPoint_h
#define MaterialDataPoint_h

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//

#include <math.h>
#include <Vector.h>
#include "MaterialTensorOperations.h"

constexpr double SMALL_VALUE = 1e-18;
constexpr double ZERO_VALUE = 1e-6;
constexpr double LARGE_VALUE = 1e18;

class DataPoint {
public:
	// variables
	double Toct = 0.0;	// octahedral shear stress
	double Pavg = 0.0;	// mean average pressure
	double Goct = 0.0;	// octahedral shear strain
	double Evol = 0.0;	// volumetric strain
	double Lode = 0.0;	// lode angle
	double I1 = 0.0;	// first invariant of the strain tensor
	double J2 = 0.0;	// second invariant of the strain deviator tensor
	double J3 = 0.0;	// third invariant of the strain deviator tensor
	double Gsec = 0.0;	// secant shear modulus
	double Bsec = 0.0;	// secant dilatancy

	// constructors
	DataPoint(void) = default;							// null constructor
	DataPoint(const DataPoint&) = default;				// copy constructor
	DataPoint(const Vector& rawdata, Vector& info);		// full constructor
		
	// destructor
	~DataPoint();

	// operator overloading
	DataPoint& operator= (const DataPoint&) = default;
	friend OPS_Stream& operator<<(OPS_Stream& s, const DataPoint& obj);

	// operational functions
	int getDim(void);
	double getVR(void);
	Vector getStress(void);
	Vector getStrain(void);
	Vector getFabric(void);
	Vector getTestInfo(void);

private:
	// material state data
	int nDim = 0;
	Vector stress = Vector(6);
	Vector strain = Vector(6);
	Vector fabric = Vector(6);
	double voidRatio = 0.0;
	double testAngle = 0.0;
	double testPref = 0.0;


	// variable determination
	void compute(void);
};
#endif