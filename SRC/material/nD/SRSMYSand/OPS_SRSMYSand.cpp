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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/OPS_SRSMYSand.cpp $
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $
//
// Description: This file contains the implementation for the OPS_SRSMYSand function.

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                       SRSMYSand nD material initializer                        |
 +                                                                                +
 |--------------------------------------------------------------------------------|
 |                                                                                |
 +         Authors: Onur Deniz Akan (IUSS),                                       +
 |                  Guido Camata, Enrico Spacone (UNICH),                         |
 |                  Carlo G. Lai (UNIPV) and                                      |
 +                  Claudio Tamagnini (UNIPG)                                     +
 |                                                                                |
 |         Istituto Universitario di Studi Superiori di Pavia (IUSS)              |
 +		   Universita degli Studi Chieti - Pescara	          (UNICH)             +
 |         Universita degli Studi di Pavia                    (UNIPV)             |
 |         Università degli Studi di Perugia                  (UNIPG)             |
 +			                                                                      +
 |                                                                                |
 |         Email: onur.akan@iusspavia.it                                          |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- November 2023                                                |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

#include "SRSMYSand.h"




void* OPS_SRSMYSand()
{
	const int numParam = 13;
	const int totParam = 26;
	int tag;
	double param[totParam];
	param[numParam] = 20;
	param[numParam + 1] = 5.0;
	param[numParam + 2] = 3.;
	param[numParam + 3] = 1.;
	param[numParam + 4] = 0.;
	param[numParam + 5] = 0.6;
	param[numParam + 6] = 0.9;
	param[numParam + 7] = 0.02;
	param[numParam + 8] = 0.7;
	param[numParam + 9] = 101.;
	param[numParam + 10] = 0.1;
	param[numParam + 11] = 0.;
	param[numParam + 12] = 1.;

	int argc = OPS_GetNumRemainingInputArgs() + 2;
	const char* arg[] = { "nd", "rho", "refShearModul",
			"refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe",
			"phaseTransformAngle", "contractionParam1",
			"contractionParam3","dilationParam1","dilationParam3",
			"numberOfYieldSurf (=20)",
			"contractionParam2=5.0", "dilationParam2=3.0",
			"liquefactionParam1=1.0", "liquefactionParam2=0.0",
			"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
			"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.1)",
			"Hv (=0)", "Pv (=1.)" };
	if (argc < (3 + numParam)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: nDMaterial PressureDependMultiYield02 tag? " << arg[0];
		opserr << "? " << "\n";
		opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? " << "\n";
		opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? " << "\n";
		opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
		opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? " << "\n";
		opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? " << "\n";
		opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? " << "\n";
		opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << "\n";
		opserr << arg[22] << "? " << arg[23] << "? " << "\n";
		return 0;
	}

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid PressureDependMultiYield02 tag" << "\n";
		return 0;
	}

	int in = 17;
	for (int i = 3; (i < argc && i < in); i++)
		if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
			opserr << "WARNING invalid " << arg[i - 3] << "\n";
			opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
			return 0;
		}

	static double* gredu = 0;

	// user defined yield surfaces
	if (param[numParam] < 0 && param[numParam] > -100) {
		param[numParam] = -int(param[numParam]);
		gredu = new double[int(2 * param[numParam])];

		for (int i = 0; i < 2 * param[numParam]; i++)
			if (OPS_GetDoubleInput(&numdata, &gredu[i]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
				return 0;
			}
	}

	if (gredu != 0) {
		for (int i = in + int(2 * param[numParam]); i < argc; i++)
			if (OPS_GetDoubleInput(&numdata, &param[i - 3 - int(2 * param[numParam])]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
				return 0;
			}
	}
	else {
		for (int i = in; i < argc; i++)
			if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial PressureDependMultiYield02: " << tag << "\n";
				return 0;
			}
	}


	SRSMYSand* temp =
		new SRSMYSand(tag, param[0], param[1], param[2],
			param[3], param[4], param[5],
			param[6], param[7], param[8],
			param[9], param[10], param[11],
			param[12], param[13], gredu, param[14],
			param[15], param[16], param[17],
			param[18], param[19], param[20], param[21],
			param[22], param[23], param[24], param[25]);

	if (gredu != 0) {
		delete[] gredu;
		gredu = 0;
	}

	return temp;
}











/*
void HelpVM(void) {

}

#define OPS_Export 
OPS_Export void*

OPS_SRSMYSand(void)
{
	// display kudos
	static int kudos = 0;
	if (++kudos == 1)
		opserr << "";

	// initialize pointers
	NDMaterial* theMaterial = nullptr;								// Pointer to an nD material to be returned

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


	// begin recieving
	int numArgs = OPS_GetNumRemainingInputArgs();
	int numData = 1;

	// check if help is requested
	if (numArgs < 3) {
		while (OPS_GetNumRemainingInputArgs() > 0) {
			const char* inputstring = OPS_GetString();
			// recive print help flag, print material command use and quit
			if (strcmp(inputstring, "-help") == 0) {
				HelpVM();
				opserr << "FATAL: OPS_SRSMYSand() - Program terminated since '-help' is the first option in the SRSMYSand nDmaterial declaration.\n";
				opserr << "                           Remove or move the '-help' option further along the declaration to continue...\n";
				exit(-1);
			}
		}
	}

	// check mandatory inputs
	if (numArgs < 7) {
		opserr << "FATAL: OPS_SRSMYSand() - please define at least -> nDMaterial SRSMYSand tag? -K Kref? -G Gref? -P Pref? for linear elastic analysis...\n";
		HelpVM();
		return theMaterial;
	}

	// input #1 - recieve unique material tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "FATAL: OPS_SRSMYSand() - invalid tag? (must be an integer)\n\n";
		HelpVM();
		return theMaterial;
	}

	// continue with the optional inputs
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();

		// input #2 - recieve material bulk modulus at reference pressure
		if (strcmp(inputstring, "-K") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Kref) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid Kref value after flag: -K (must be double)\n";
				return theMaterial;
			}
		}

		// input #3 - recieve material shear modulus at reference pressure
		if (strcmp(inputstring, "-G") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Gref) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid Gref value after flag: -G (must be double)\n";
				return theMaterial;
			}
		}

		// input #4 - recieve material reference pressure
		if (strcmp(inputstring, "-P") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Pref) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid Pref value after flag: -P (must be double)\n";
				return theMaterial;
			}
		}

		// input #5 - recieve material mass density
		if (strcmp(inputstring, "-R") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &rho) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid rho value after flag: -r (must be double)\n";
				return theMaterial;
			}
		}

		// input #6 - recieve material elastic modulus update power
		if (strcmp(inputstring, "-M") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &modn) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid modn value after flag: -m (must be double)\n";
				return theMaterial;
			}
		}

		// input #7 - recieve total number of yield surfaces
		if (strcmp(inputstring, "-t") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, &TNYS) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid TNYS value after flag: -T (must be integer)\n";
				return theMaterial;
			}
			if (abs(TNYS) > maxTNYS || TNYS == 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid TNYS value! (must be an integer TNYS <= " << maxTNYS << " and TNYS cannot be 0)\n";
				return theMaterial;
			}
			TNYS = abs(TNYS);	// make sure TNYS is a positive integer
		}

		// input #8 - recieve yield surface cohesion
		if (strcmp(inputstring, "-c") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &cohesion) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid cohesion value after flag: -c (must be double)\n";
				return theMaterial;
			}
		}

		// input #9 - recieve yield surface friction angle
		if (strcmp(inputstring, "-f") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &frictionAngle) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid frictionAngle value after flag: -f (must be double)\n";
				return theMaterial;
			}
		}

		// input #10 - recieve yield surface dilatancy angle
		if (strcmp(inputstring, "-d") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dilatancyAngle) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid dilatancyAngle value after flag: -d (must be double)\n";
				return theMaterial;
			}
		}

		// input #11 - recieve yield surface peak shear strain
		if (strcmp(inputstring, "-s") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &peakShearStrain) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid peakShearStrain value after flag: -s (must be double)\n";
				return theMaterial;
			}
		}

		// input #12 - recieve material integration type
		if (strcmp(inputstring, "-implex") == 0) {
			integrationType = 1;
		}

		// input #13 - recieve yield surface Gsec/Gmax ratios
		if (strcmp(inputstring, "-hModuli") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "FATAL: OPS_SRSMYSand() - please enter G Modulus vector for each yield surface after flag: -hModuli (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HModuli = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HModuli) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid G Moduli! (must be " << TNYS << " many doubles  after flag: -hModuli)\n";
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
				opserr << "FATAL: OPS_SRSMYSand() - please enter Dilatancy vector for each yield surface after flag: -hDilatancies (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HDilatancies = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HDilatancies) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid Dilatancies! (must be " << TNYS << " many doubles  after flag: -hDilatancies)\n";
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
				opserr << "FATAL: OPS_SRSMYSand() - please enter vector of Hardening Parameters for each yield surface after flag: -hStrains (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HParams = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HParams) < 0) {
				opserr << "FATAL: OPS_SRSMYSand() - invalid Hardening Parameters! (must be " << TNYS << " many doubles  after flag: -hStrains)\n";
				return theMaterial;
			}
			hps = Vector(TNYS);
			for (int i = 0; i < TNYS; i++) { hps(i) = HParams[i]; }
		}

		// operational flags
			// recive print help flag, print material command use, make verbose and continue
		if (strcmp(inputstring, "-help") == 0) {
			HelpVM();
		}

		// recive verbosity flag
		if (strcmp(inputstring, "-debug") == 0) {
			beVerbose = true;
		}

	} // all inputs recieved!

	// check mandatory parameters
	if (Kref == 0) {
		opserr << "FATAL: OPS_SRSMYSand() - please enter a valid Kref value! Current Kref = " << Kref << ". Kref cannot be 0!\n";
		HelpVM();
		return theMaterial;
	}

	if (Gref == 0) {
		opserr << "FATAL: OPS_SRSMYSand() - please enter a valid Gref value! Current Gref = " << Gref << ". Gref cannot be 0!\n";
		HelpVM();
		return theMaterial;
	}

	if (Pref < 101) {
		// Pref is kPa in most cases, though it is not a must. Just invite the user to double check the input value  
		opserr << "WARNING: OPS_SRSMYSand() - please double check the Pref value in use! Current Pref = " << Pref << "!\n";
	}

	// create a SRSMYSand nDmaterial object
	theMaterial = new SRSMYSand();

	return theMaterial;
}
*/