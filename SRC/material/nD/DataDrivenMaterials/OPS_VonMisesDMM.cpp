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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/OPS_VonMisesDMM.cpp $
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	June 2022
//
// Description: This file contains the implementation for the OPS_VonMisesDMM function.

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                     OPS_VonMisesDMM initializer for the                        |
 |                von Mises data-driven multi-scale nD material                   |
 |                                                                                |
 +--------------------------------------------------------------------------------+
 |                                                                                |
 |         Authors: Onur Deniz Akan (IUSS),                                       |
 +                  Guido Camata, Enrico Spacone (UNICH)                          +
 |                  and Carlo G. Lai (UNIPV)                                      |
 |                                                                                |
 +      Istituto Universitario di Studi Superiori di Pavia          (IUSS)        +
 |		Universita degli Studi 'G. d'Annunzio' Chieti - Pescara	    (UNICH)       |
 |      Universita degli Studi di Pavia                             (UNIPV)       |
 +			                                                                      +
 |                                                                                |
 |           Email: onur.akan@iusspavia.it                                        |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- June 2022                                                    |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

 /*
 | Argument               | Type             | Description
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | $tag                   | integer          |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | $Kref, $Gref, $Pref    | 3 double         |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -R $Rho                | string + double  |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -M $Modn               | string + double  |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -t $tnys               | string + integer |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -c $cohesion           | string + double  |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -f $frictionAngle      | string + double  |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -s $peakShearAngle     | string + double  |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -implex                | string           |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -ddType $flag          | string + integer |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -hRatio $HModuli1      | string + list    |
 |         $HModuli2 ...  |       of doubles |
 |                        |                  |
 |------------------------|------------------|--------------------------------------------------------------------------------------
 | -hStrain $HParams1     | string + list    |
 |          $HParams2 ... |       of doubles |
 |                        |                  |
 */

#include "VonMisesDMM.h"

#define OPS_Export 
OPS_Export void*

OPS_VonMisesDMM(void)
{
	// display kudos
	static int kudos = 0;
	if (++kudos == 1)
		opserr << "von Mises Data-Driven Multi-scale nDmaterial - Written: OD.Akan, G.Camata, E.Spacone, CG.Lai, IUSS Pavia \n";

	// initialize pointers
	NDMaterial* theMaterial = nullptr;				// Pointer to an nD material to be returned
	DataDrivenNestedSurfaces* theData = nullptr;	// Pointer to a nested yield surface object

	// initialize material parameters
	int maxTNYS = 1000;
	int tag = 0; int TNYS = 0; double rho = 0.;
	double Kref = 0.; double Gref = 0.; double Pref = 0.; double modn = 0.;
	double cohesion = 0.; double frictionAngle = 0.; double peakShearStrain= 0.;

	// user defined yield surfaces
	double* HParams = nullptr; 
	double* HModuli = nullptr;

	// other inputs
	int dataDriverType = 0;		// use hyperbolic backbone by default
	int integrationType = 0;	// use implicit integration by default


	// begin recieving
	int numArgs = OPS_GetNumRemainingInputArgs();
	int numData = 1;

	// recieve mandatory inputs
	if (numArgs < 4) {	
		opserr << "OPS_VonMisesDMM: \n";
		opserr << "Please define at least: nDMaterial VonMisesDMM tag? Kref? Gref? Pref? for linear elastic analysis...\n\n";
		opserr << "Want: nDMaterial VonMisesDMM tag? Kref? Gref? Pref? <-R Rho? -M Modn? -t tnys? -c cohesion? -f frictionAngle? -s peakShearStrain? -implex?>\n";
		opserr << "         and other optional parameters for data-driven material modeling:\n";
		opserr << "         <-ddType flag?> <-hRatio $HModuli1 $HModuli2 $HModuli3 ... -hStrain $HParams1 $HParams2 $HParams3 ...>\n";
		return theMaterial;
	}

	// input #1 - recieve unique material tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "OPS_VonMisesDMM: WARNING: invalid tag? (must be an integer)" << endln;
		return theMaterial;
	}

	// input #2 - recieve material bulk modulus at reference pressure
	if (OPS_GetDoubleInput(&numData, &Kref) < 0) {
		opserr << "OPS_VonMisesDMM: WARNING: invalid Kref value after flag: -K (must be double)\n";
		return theMaterial;
	}

	// input #3 - recieve material shear modulus at reference pressure
	if (OPS_GetDoubleInput(&numData, &Gref) < 0) {
		opserr << "OPS_VonMisesDMM: WARNING: invalid Gref value after flag: -G (must be double)\n";
		return theMaterial;
	}

	// input #4 - recieve material reference pressure
	if (OPS_GetDoubleInput(&numData, &Pref) < 0) {
		opserr << "OPS_VonMisesDMM: WARNING: invalid Pref value after flag: -P (must be double)\n";
		return theMaterial;
	}

	// continue with the optional inputs
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();

		// input #5 - recieve material mass density
		if (strcmp(inputstring, "-R") == 0) {
			if (OPS_GetDoubleInput(&numData, &rho) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid rho value after flag: -r (must be double)\n";
				return theMaterial;
			}
		}

		// input #6 - recieve material elastic modulus update power
		if (strcmp(inputstring, "-M") == 0) {
			if (OPS_GetDoubleInput(&numData, &modn) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid modn value after flag: -m (must be double)\n";
				return theMaterial;
			}
		}

		// input #7 - recieve total number of yield surfaces
		if (strcmp(inputstring, "-t") == 0) {
			if (OPS_GetIntInput(&numData, &TNYS) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid TNYS value after flag: -T (must be integer)\n";
				return theMaterial;
			}
			if (abs(TNYS) > maxTNYS || TNYS == 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid TNYS value! (must be an integer -" << maxTNYS << " <= TNYS <= " << maxTNYS << " && TNYS != 0)\n";
				return theMaterial;
			}
		}

		// input #8 - recieve yield surface cohesion
		if (strcmp(inputstring, "-c") == 0) {
			if (OPS_GetDoubleInput(&numData, &cohesion) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid cohesion value after flag: -c (must be double)\n";
				return theMaterial;
			}
		}

		// input #9 - recieve yield surface friction angle
		if (strcmp(inputstring, "-f") == 0) {
			if (OPS_GetDoubleInput(&numData, &frictionAngle) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid frictionAngle value after flag: -f (must be double)\n";
				return theMaterial;
			}
		}

		// input #10 - recieve yield surface peak shear strain
		if (strcmp(inputstring, "-s") == 0) {
			if (OPS_GetDoubleInput(&numData, &peakShearStrain) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid peakShearStrain value after flag: -s (must be double)\n";
				return theMaterial;
			}
		}

		// input #11 - recieve material integration type
		if (strcmp(inputstring, "-implex") == 0) {
			integrationType = 1;
		}

		// input #12 - recieve yield surface data driver type
		if (strcmp(inputstring, "-ddType") == 0) {
			if (OPS_GetIntInput(&numData, &dataDriverType) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid data driver type after flag: -ddType (must be integer)\n"; 
				opserr << "OPS_VonMisesDMM: WARNING:                                  [0: backbone, 1: offline, 2:online]\n";
				return theMaterial;
			}
		}

		// input #13 - recieve yield surface G/Gmax ratios
		if (strcmp(inputstring, "-hRatio") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "Want: nDMaterial VonMisesDMM: please enter $HModuli vector for each yield surface after flag: -hRatio (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HModuli = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HModuli) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid HModuli! (must be " << TNYS << " many doubles  after flag: -hRatio)\n";
				return theMaterial;
			}

		}

		// input #14 - recieve yield surface hardening parameters (strain discretization)
		if (strcmp(inputstring, "-hStrain") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "Want: nDMaterial VonMisesDMM: please enter HParams vector for each yield surface after flag: -hStrain (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HParams = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HParams) < 0) {
				opserr << "OPS_VonMisesDMM: WARNING: invalid HParams! (must be " << TNYS << " many doubles  after flag: -hStrain)\n";
				return theMaterial;
			}
		}

	} // all inputs recieved!
	

	// create a nested yield surface object
	theData = new DataDrivenNestedSurfaces(cohesion, frictionAngle, peakShearStrain, HModuli, HParams);

	if (theData == nullptr) {
		opserr << "FATAL: OPS_VonMisesDMM: cannot create VonMisesDMM material with tag: " << tag << "\n";
		opserr << "FATAL: OPS_VonMisesDMM: yield surface data yielded bad result...";
		exit(-1);
	}

	// create a VonMisesDMM nDmaterial object
	theMaterial = new VonMisesDMM(tag, rho, Kref, Gref, Pref, modn, TNYS, theData, dataDriverType, integrationType);

	if (theMaterial == nullptr) {
		opserr << "FATAL: OPS_VonMisesDMM: cannot create VonMisesDMM material with tag: " << tag << "\n";
		exit(-1);
	}

	// free the memory
	if (HModuli != nullptr)
	{
		delete[] HModuli;
		HModuli = nullptr;
	}
	if (HParams != nullptr)
	{
		delete[] HParams;
		HParams = nullptr;
	}

	theData = nullptr;

	return theMaterial;
}