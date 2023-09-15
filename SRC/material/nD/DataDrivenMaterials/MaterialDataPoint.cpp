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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialDataPoint.cpp$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//
// Description: This file contains the implementation for the DataPoint class.

#include "MaterialDataPoint.h"

DataPoint::DataPoint(const Vector& rawdata, Vector& info)
{
	// find out data dimension
	int size = rawdata.Size();

	// recieve test info 
	testPref = info(0);
	testAngle = info(1);

	if (size == 14) {
		// incoming 2D vector legend:
		// <sig11 sig22 sig12 sig21> <eps11 eps22 eps12 eps21> <phi11 phi22 phi12 phi21> vR null
		nDim = 2;
		// initialize data vectors
		stress.resize(3);
		strain.resize(3);
		fabric.resize(3);
		for (int i = 0; i < 3; i++) { 
			stress(i) = rawdata(i);
			strain(i) = rawdata(i + 4);
			fabric(i) = rawdata(i + 8);
		}
		voidRatio = rawdata(12);
	}
	else if (size == 29 ) {
		// incoming 2D vector legend:
		// <sig11 sig22 sig33 sig23 sig13 sig12 sig32 sig31 sig21> <eps11 eps22 eps33 eps23 eps13 eps12 eps32 eps31 eps21>  
		// <phi11 phi22 phi33 phi23 phi13 phi12 phi32 phi31 phi21> vR null
		nDim = 3;
		// initialize data vectors
		stress.resize(6);
		strain.resize(6);
		fabric.resize(6);
		for (int i = 0; i < 6; i++) {
			stress(i) = rawdata(i);
			strain(i) = rawdata(i + 9);
			fabric(i) = rawdata(i + 18);
		}
	} 
	else {
		opserr << "FATAL: DataPoint() - unsupported material dimension! must be either 2 or 3...";
		exit(-1);
	}

	// compute material variables
	compute();
}

DataPoint::~DataPoint()
{
}

OPS_Stream& operator<<(OPS_Stream& s, const DataPoint& obj) {
	
	s << "Data-point summary:             \n";
	s << "--------------------------------\n";
	s << "Oct. shear stress : " << obj.Toct << "\n";
	s << "Isotropic pressure: " << obj.Pavg << "\n";
	s << "Oct. shear strain : " << obj.Goct << "\n";
	s << "Volumetric strain : " << obj.Evol << "\n";

	return s << endln;
}

// operational functions
int DataPoint::getDim(void) { return nDim; }
double DataPoint::getVR(void) { return voidRatio; }
Vector DataPoint::getStress(void) { return stress; }

Vector DataPoint::getStrain(void) {

	opserr << "DataPoint::getStrain() - works until here!\n"; 
	opserr << this;
	return strain;

}

Vector DataPoint::getFabric(void) { return fabric; }
Vector DataPoint::getTestInfo(void) { Vector info(2); info(0) = testPref; info(1) = testAngle;  return info; }

// variable determination
void DataPoint::compute(void) {
	if (nDim == 2) {
		// do computation for 2D material
			// using stress tensor
		Pavg = 0.5 * (stress(0) + stress(1));
		Vector deviator = stress; deviator(0) -= Pavg; deviator(1) -= Pavg;
		Toct = sqrt(1.0 / 3.0 * (deviator(0) * deviator(0) + deviator(1) * deviator(1) + 2 * deviator(2) * deviator(2)));
		Lode = 0.0;

			// using strain tensor
		Evol = strain(0) + strain(1);
		deviator = strain; deviator(0) -= Evol * 0.5; deviator(1) -= Evol * 0.5;
		Goct = 2 * sqrt(1.0 / 3.0 * (deviator(0) * deviator(0) + deviator(1) * deviator(1) + 2 * deviator(2) * deviator(2)));
		I1 = 0.0;
		J2 = 0.0;
		J3 = 0.0;

			// moduli
		Bsec = 0.0;
		Gsec = 0.0;
	}
	else {
		// do computation for 3D material

	}
}