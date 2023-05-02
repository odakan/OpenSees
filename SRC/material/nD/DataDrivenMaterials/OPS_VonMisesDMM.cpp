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

#include "VonMisesDMM.h"

void Help(void) {
	opserr << "\n";
	opserr << "vonMisesDMM Help: \n";
	opserr << "----------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	opserr << "nDMaterial VonMisesDMM tag? -K Kref? -G Gref? -P Pref? <-R Rho? -M Modn? -t tnys? -c cohesion? -f frictionAngle? -d dilatancyAngle -s peakShearStrain? -implex?>\n";
	opserr << "         and other optional parameters for data-driven material modeling:\n";
	opserr << "         <-ddType flag?> <-HModuli $HModuli1 $HModuli2 $HModuli3 ... -hStrain $HParams1 $HParams2 $HParams3 ...>\n";
	opserr << "\n\n";
	opserr << "	| Argument               | Type             | Description                                                                            |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| $tag                   | integer          |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| $Kref, $Gref, $Pref    | 3 double         |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -R $Rho                | string + double  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -M $Modn               | string + double  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -t $tnys               | string + integer |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -c $cohesion           | string + double  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -f $frictionAngle      | string + double  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -s $peakShearAngle     | string + double  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -implex                | string           |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -ddType $flag          | string + double  | $flag is a double in the domain : [0, 1] . 0 : hyperbolic backbone, 1 : offline wheras |\n";
	opserr << "	|                        |                  | a value is between 0 and 1 is the strain dicretization step for the online method      |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -hRatio $HModuli1      | string + list    |                                                                                        |\n";
	opserr << "	|         $HModuli2 ...  | of doubles       |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -hStrain $HParams1     | string + list    |                                                                                        |\n";
	opserr << "	|          $HParams2 ... | of doubles       |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "\n\n";
}

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
	bool beVerbose = false;
	int maxTNYS = 250;
	int tag = 0; int TNYS = 0; double rho = 0.0;
	double Kref = 0.0; double Gref = 0.0; double Pref = 0.0; double modn = 0.0;
	double cohesion = 0.0; double frictionAngle = 0.0; double dilatancyAngle = 0.0; double peakShearStrain = 0.0;

	// user defined yield surfaces
	double* HParams = nullptr; Vector hps(1);
	double* HModuli = nullptr; Vector hmod(1);
	double* HDilatancies = nullptr; Vector dps(1);

	// other inputs
	int integrationType = 0;		// use implicit integration by default
	double dataDriverType = 0;		// use hyperbolic backbone by default


	// begin recieving
	int numArgs = OPS_GetNumRemainingInputArgs();
	int numData = 1;

	// check if help is requested
	if (numArgs < 3) {
		while (OPS_GetNumRemainingInputArgs() > 0) {
			const char* inputstring = OPS_GetString();
			// recive print help flag, print material command use and quit
			if (strcmp(inputstring, "-help") == 0) {
				Help();
				exit(-1);
			}
		}
	}

	// check mandatory inputs
	if (numArgs < 7) {
		opserr << "FATAL: OPS_VonMisesDMM() - please define at least -> nDMaterial VonMisesDMM tag? -K Kref? -G Gref? -P Pref? for linear elastic analysis...\n";
		Help();
		return theMaterial;
	}

	// input #1 - recieve unique material tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "FATAL: OPS_VonMisesDMM() - invalid tag? (must be an integer)\n";
		return theMaterial;
	}

	// continue with the optional inputs
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();

		// input #2 - recieve material bulk modulus at reference pressure
		if (strcmp(inputstring, "-K") == 0) {
			if (OPS_GetDoubleInput(&numData, &Kref) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid Kref value after flag: -K (must be double)\n";
				return theMaterial;
			}
		}

		// input #3 - recieve material shear modulus at reference pressure
		if (strcmp(inputstring, "-G") == 0) {
			if (OPS_GetDoubleInput(&numData, &Gref) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid Gref value after flag: -G (must be double)\n";
				return theMaterial;
			}
		}

		// input #4 - recieve material reference pressure
		if (strcmp(inputstring, "-P") == 0) {
			if (OPS_GetDoubleInput(&numData, &Pref) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid Pref value after flag: -P (must be double)\n";
				return theMaterial;
			}
		}

		// input #5 - recieve material mass density
		if (strcmp(inputstring, "-R") == 0) {
			if (OPS_GetDoubleInput(&numData, &rho) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid rho value after flag: -r (must be double)\n";
				return theMaterial;
			}
		}

		// input #6 - recieve material elastic modulus update power
		if (strcmp(inputstring, "-M") == 0) {
			if (OPS_GetDoubleInput(&numData, &modn) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid modn value after flag: -m (must be double)\n";
				return theMaterial;
			}
		}

		// input #7 - recieve total number of yield surfaces
		if (strcmp(inputstring, "-t") == 0) {
			if (OPS_GetIntInput(&numData, &TNYS) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid TNYS value after flag: -T (must be integer)\n";
				return theMaterial;
			}
			if (abs(TNYS) > maxTNYS || TNYS == 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid TNYS value! (must be an integer TNYS <= " << maxTNYS << " and TNYS cannot be 0)\n";
				return theMaterial;
			}
			TNYS = abs(TNYS);	// make sure TNYS is a positive integer
		}

		// input #8 - recieve yield surface cohesion
		if (strcmp(inputstring, "-c") == 0) {
			if (OPS_GetDoubleInput(&numData, &cohesion) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid cohesion value after flag: -c (must be double)\n";
				return theMaterial;
			}
		}

		// input #9 - recieve yield surface friction angle
		if (strcmp(inputstring, "-f") == 0) {
			if (OPS_GetDoubleInput(&numData, &frictionAngle) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid frictionAngle value after flag: -f (must be double)\n";
				return theMaterial;
			}
		}

		// input #10 - recieve yield surface dilatancy angle
		if (strcmp(inputstring, "-d") == 0) {
			if (OPS_GetDoubleInput(&numData, &dilatancyAngle) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid dilatancyAngle value after flag: -d (must be double)\n";
				return theMaterial;
			}
		}

		// input #11 - recieve yield surface peak shear strain
		if (strcmp(inputstring, "-s") == 0) {
			if (OPS_GetDoubleInput(&numData, &peakShearStrain) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid peakShearStrain value after flag: -s (must be double)\n";
				return theMaterial;
			}
		}

		// input #12 - recieve material integration type
		if (strcmp(inputstring, "-implex") == 0) {
			integrationType = 1;
		}

		// input #13 - recieve yield surface data driver type
		if (strcmp(inputstring, "-ddType") == 0) {
			if (OPS_GetDoubleInput(&numData, &dataDriverType) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid data driver type after flag: -ddType (must be integer)\n";
				opserr << "FATAL: OPS_VonMisesDMM() - [$val is a double. -1: user custom surface, 0 : hyperbolic backbone, 1 : offline wheras\n";
				opserr << "FATAL: OPS_VonMisesDMM() -         a value between 0 and 1 is the strain dicretization step for the online method]\n";
				return theMaterial;
			}
		}

		// input #14 - recieve yield surface Gsec/Gmax ratios
		if (strcmp(inputstring, "-hModuli") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "FATAL: OPS_VonMisesDMM() - please enter G Modulus vector for each yield surface after flag: -hModuli (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HModuli = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HModuli) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid G Moduli! (must be " << TNYS << " many doubles  after flag: -hModuli)\n";
				return theMaterial;
			}
			hmod = Vector(TNYS);
			for (int i = 0; i < TNYS; i++) { hmod(i) = HModuli[i]; }
		}

		// input #15 - recieve yield surface dilatancy
		if (strcmp(inputstring, "-hDilatancies") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "FATAL: OPS_VonMisesDMM() - please enter Dilatancy vector for each yield surface after flag: -hDilatancies (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HDilatancies = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HDilatancies) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid Dilatancies! (must be " << TNYS << " many doubles  after flag: -hDilatancies)\n";
				return theMaterial;
			}
			dps = Vector(TNYS);
			for (int i = 0; i < TNYS; i++) { dps(i) = HDilatancies[i]; }
		}

		// input #16 - recieve yield surface hardening parameters (strain discretization)
		if (strcmp(inputstring, "-hStrains") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "FATAL: OPS_VonMisesDMM() - please enter vector of Hardening Parameters for each yield surface after flag: -hStrains (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HParams = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HParams) < 0) {
				opserr << "FATAL: OPS_VonMisesDMM() - invalid Hardening Parameters! (must be " << TNYS << " many doubles  after flag: -hStrains)\n";
				return theMaterial;
			}
			hps = Vector(TNYS);
			for (int i = 0; i < TNYS; i++) { hps(i) = HParams[i]; }
		}

		// operational flags
			// recive print help flag, print material command use, make verbose and continue
		if (strcmp(inputstring, "-help") == 0) {
			Help();
		}

			// recive verbosity flag
		if (strcmp(inputstring, "-debug") == 0) {
			beVerbose = true;
		}

	} // all inputs recieved!

	// check mandatory parameters
	if (Kref == 0) {
		opserr << "FATAL: OPS_VonMisesDMM() - please enter a valid Kref value! Current Kref = " << Kref << ". Kref cannot be 0!\n";
		Help();
		return theMaterial;
	}

	if (Gref == 0) {
		opserr << "FATAL: OPS_VonMisesDMM() - please enter a valid Gref value! Current Gref = " << Gref << ". Gref cannot be 0!\n";
		Help();
		return theMaterial;
	}

	if (Pref < 101) {
		// Pref is kPa in most cases, though it is not a must. Just invite the user to double check the input value  
		opserr << "WARNING: OPS_VonMisesDMM() - please double check the Pref value in use! Current Pref = " << Pref << "!\n";
	}

	if (beVerbose) {
		opserr << "\nOPS_VonMisesDMM:\n";
		opserr << "-----------------\n";
		opserr << "tag             : " << tag << "\n";
		opserr << "rho             : " << rho << "\n";
		opserr << "Kref            : " << Kref << "\n";
		opserr << "Gref            : " << Gref << "\n";
		opserr << "Pref            : " << Pref << "\n";
		opserr << "Modn            : " << modn << "\n";
		opserr << "TNYS            : " << TNYS << "\n";
		opserr << "cohesion        : " << cohesion << "\n";
		opserr << "frictionAngle   : " << frictionAngle << "\n";
		opserr << "dilatancyAngle  : " << dilatancyAngle << "\n";
		opserr << "peakShearStrain : " << peakShearStrain << "\n";
		opserr << "driverType      : " << dataDriverType << "\n";
		opserr << "integrationType : " << integrationType << "\n";
		if (HModuli != nullptr) {
			opserr << "HModuli         : " << hmod;
		}
		if (HParams != nullptr) {
			opserr << "HParams         : " << hps;
		}
		if (HDilatancies != nullptr) {
			opserr << "HDilatancies    : " << dps;
		}
		opserr << "\n";
	}
		
	// create a nested yield surface object
	theData = new DataDrivenNestedSurfaces(tag, cohesion, frictionAngle, dilatancyAngle, peakShearStrain, TNYS, dataDriverType, hmod, hps, dps, beVerbose);
	
	if (theData == nullptr) {
		opserr << "FATAL: OPS_VonMisesDMM() - cannot create VonMisesDMM material with tag: " << tag << "\n";
		opserr << "FATAL: OPS_VonMisesDMM() - yield surface data yielded bad result...";
		return theMaterial;
	}
	
	// create a VonMisesDMM nDmaterial object
	theMaterial = new VonMisesDMM(tag, rho, Kref, Gref, Pref, modn, theData, dataDriverType, integrationType, beVerbose);

	if (theMaterial == nullptr) {
		opserr << "FATAL: OPS_VonMisesDMM() - cannot create VonMisesDMM material with tag: " << tag << "\n";
		return theMaterial;
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

	// break link after passing theData to theMaterial
	theData = nullptr;
	
	return theMaterial;
}