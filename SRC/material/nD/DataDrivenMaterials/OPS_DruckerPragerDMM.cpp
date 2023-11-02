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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/OPS_DruckerPragerDMM.cpp $
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	June 2022
//
// Description: This file contains the implementation for the OPS_DruckerPragerDMM function.

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                  OPS_DruckerPragerDMM initializer for the                      |
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

#include <fstream>
#include <string>


#include "DruckerPragerDMM.h"


void HelpDP(void) {
	opserr << "\n";
	opserr << "DruckerPragerDMM Help: \n";
	opserr << "----------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	opserr << "nDMaterial DruckerPragerDMM tag? -K Kref? -G Gref? -P Pref? <-R Rho? -M Modn? -t tnys? -c cohesion? -f frictionAngle? -d dilatancyAngle -s peakShearStrain? -implex?>\n";
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
	opserr << "	| -ddType $type $step    | string + string  | $flag is a double in the domain : [0, 1] . 0 : hyperbolic backbone, 1 : offline wheras |\n";
	opserr << "	|                        | + double         | a value is between 0 and 1 is the strain dicretization step for the online method      |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -hRatio $HModuli1      | string + list    |                                                                                        |\n";
	opserr << "	|         $HModuli2 ...  | of doubles       |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -hStrain $HParams1     | string + list    |                                                                                        |\n";
	opserr << "	|          $HParams2 ... | of doubles       |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "	| -dbPath $Path          | string           | absolute path of the database folder with experiments                                  |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|                        |                  |                                                                                        |\n";
	opserr << "	|------------------------------------------------------------------------------------------------------------------------------------|\n";
	opserr << "\n\n";
}

#define OPS_Export 
OPS_Export void*

OPS_DruckerPragerDMM(void)
{
	// display kudos
	static int kudos = 0;
	if (++kudos == 1)
		opserr << "Drucker-Prager Data-Driven Multi-scale nDmaterial - Written: OD.Akan, G.Camata, E.Spacone, CG.Lai, IUSS Pavia \n";

	// initialize pointers
	NDMaterial* theMaterial = nullptr;								// Pointer to an nD material to be returned
	std::shared_ptr<DataDrivenNestedSurfaces> theData = nullptr;	// Pointer to a nested yield surface object

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
	int nFiles = 0;
	int integrationType = 0;		// use implicit integration by default
	double dataDriverType = 0;		// use hyperbolic backbone by default
	const char* dbPathDir = nullptr;		// pointer to the database directory
	const char* dbMainFile = nullptr;		// pointer to the database main file
	char* fullPath = nullptr;
	bool isDatabase = false;


	// begin recieving
	int numArgs = OPS_GetNumRemainingInputArgs();
	int numData = 1;

	// check if help is requested
	if (numArgs < 3) {
		while (OPS_GetNumRemainingInputArgs() > 0) {
			const char* inputstring = OPS_GetString();
			// recive print help flag, print material command use and quit
			if (strcmp(inputstring, "-help") == 0) {
				HelpDP();
				opserr << "FATAL: OPS_DruckerPragerDMM() - Program terminated since '-help' is the first option in the DruckerPragerDMM nDmaterial declaration.\n";
				opserr << "                           Remove or move the '-help' option further along the declaration to continue...\n";
				exit(-1);
			}
		}
	}

	// check mandatory inputs
	if (numArgs < 7) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - please define at least -> nDMaterial DruckerPragerDMM tag? -K Kref? -G Gref? -P Pref? for linear elastic analysis...\n";
		HelpDP();
		return theMaterial;
	}

	// input #1 - recieve unique material tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - invalid tag? (must be an integer)\n\n";
		HelpDP();
		return theMaterial;
	}

	// continue with the optional inputs
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();

		// input #2 - recieve material bulk modulus at reference pressure
		if (strcmp(inputstring, "-K") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Kref) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid Kref value after flag: -K (must be double)\n";
				return theMaterial;
			}
		}

		// input #3 - recieve material shear modulus at reference pressure
		if (strcmp(inputstring, "-G") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Gref) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid Gref value after flag: -G (must be double)\n";
				return theMaterial;
			}
		}

		// input #4 - recieve material reference pressure
		if (strcmp(inputstring, "-P") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Pref) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid Pref value after flag: -P (must be double)\n";
				return theMaterial;
			}
		}

		// input #5 - recieve material mass density
		if (strcmp(inputstring, "-R") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &rho) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid rho value after flag: -r (must be double)\n";
				return theMaterial;
			}
		}

		// input #6 - recieve material elastic modulus update power
		if (strcmp(inputstring, "-M") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &modn) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid modn value after flag: -m (must be double)\n";
				return theMaterial;
			}
		}

		// input #7 - recieve total number of yield surfaces
		if (strcmp(inputstring, "-t") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, &TNYS) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid TNYS value after flag: -T (must be integer)\n";
				return theMaterial;
			}
			if (abs(TNYS) > maxTNYS || TNYS == 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid TNYS value! (must be an integer TNYS <= " << maxTNYS << " and TNYS cannot be 0)\n";
				return theMaterial;
			}
			TNYS = abs(TNYS);	// make sure TNYS is a positive integer
		}

		// input #8 - recieve yield surface cohesion
		if (strcmp(inputstring, "-c") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &cohesion) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid cohesion value after flag: -c (must be double)\n";
				return theMaterial;
			}
		}

		// input #9 - recieve yield surface friction angle
		if (strcmp(inputstring, "-f") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &frictionAngle) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid frictionAngle value after flag: -f (must be double)\n";
				return theMaterial;
			}
		}

		// input #10 - recieve yield surface dilatancy angle
		if (strcmp(inputstring, "-d") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dilatancyAngle) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid dilatancyAngle value after flag: -d (must be double)\n";
				return theMaterial;
			}
		}

		// input #11 - recieve yield surface peak shear strain
		if (strcmp(inputstring, "-s") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &peakShearStrain) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid peakShearStrain value after flag: -s (must be double)\n";
				return theMaterial;
			}
		}

		// input #12 - recieve material integration type
		if (strcmp(inputstring, "-implex") == 0) {
			integrationType = 1;
		}

		// input #13 - recieve yield surface data driver type
		if (strcmp(inputstring, "-ddType") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dataDriverType) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid data driver type after flag: -ddType (must be integer)\n";
				opserr << "FATAL: OPS_DruckerPragerDMM() - [$val is a double. -1: user custom surface, 0 : hyperbolic backbone, 1 : offline wheras\n";
				opserr << "FATAL: OPS_DruckerPragerDMM() -         a value between 0 and 1 is the strain dicretization step for the online method]\n";
				return theMaterial;
			}
		}

		// input #14 - recieve yield surface Gsec/Gmax ratios
		if (strcmp(inputstring, "-hModuli") == 0) {
			// do some checks
			int numArgs = OPS_GetNumRemainingInputArgs();
			if (numArgs < abs(TNYS)) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - please enter G Modulus vector for each yield surface after flag: -hModuli (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HModuli = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HModuli) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid G Moduli! (must be " << TNYS << " many doubles  after flag: -hModuli)\n";
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
				opserr << "FATAL: OPS_DruckerPragerDMM() - please enter Dilatancy vector for each yield surface after flag: -hDilatancies (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HDilatancies = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HDilatancies) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid Dilatancies! (must be " << TNYS << " many doubles  after flag: -hDilatancies)\n";
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
				opserr << "FATAL: OPS_DruckerPragerDMM() - please enter vector of Hardening Parameters for each yield surface after flag: -hStrains (must be doubles)\n";
				return theMaterial;
			}
			// recive
			numData = int(abs(TNYS));
			HParams = new double[abs(TNYS)];
			if (OPS_GetDoubleInput(&numData, HParams) < 0) {
				opserr << "FATAL: OPS_DruckerPragerDMM() - invalid Hardening Parameters! (must be " << TNYS << " many doubles  after flag: -hStrains)\n";
				return theMaterial;
			}
			hps = Vector(TNYS);
			for (int i = 0; i < TNYS; i++) { hps(i) = HParams[i]; }
		}

		// input #17 - recieve the database path
		if (strcmp(inputstring, "-dbPath") == 0) {
			if (OPS_GetNumRemainingInputArgs() >= 4) {
				for (int i = 0; i < 2; i++) {
					inputstring = OPS_GetString();
					if (strcmp(inputstring, "directory") == 0) {
						dbPathDir = OPS_GetString();
					}
					else if (strcmp(inputstring, "filename") == 0)
					{
						dbMainFile = OPS_GetString();
					}
				}
			}

			// Calculate the length of the concatenated string
			size_t totalLength = strlen(dbPathDir) + strlen(dbMainFile) + 1; // +1 for the null terminator

			if (totalLength > 1) {
				// Allocate memory for the full path
				fullPath = new char[totalLength];
				// Copy the directory to the full path
				strcpy(fullPath, dbPathDir);
				// Concatenate the filename to the full path
				strcat(fullPath, dbMainFile);
				try {
					// Create an ifstream object to read the file
					std::ifstream file_stream(fullPath);
					// Database is present
					isDatabase = true;
					// Check if the file opened successfully
					if (!file_stream.is_open()) {
						opserr << "FATAL: OPS_DruckerPragerDMM() - failed to open the file at: " << fullPath << "\n";
						exit(-1); // Return an error code
					}
					// Read and print the contents of the file line by line
					std::string line;
					while (std::getline(file_stream, line)) {
						nFiles++;
					}
					// Close the file when done
					file_stream.close();
				}
				catch (const std::exception&) {
					continue;
				}

				// Free fullPath
				delete[] fullPath;
				fullPath = nullptr;
			}

			// Load the response database in to the memory as a linked list

		}

		// operational flags
			// recive print help flag, print material command use, make verbose and continue
		if (strcmp(inputstring, "-help") == 0) {
			HelpDP();
		}

		// recive verbosity flag
		if (strcmp(inputstring, "-debug") == 0) {
			beVerbose = true;
		}

	} // all inputs recieved!

	// check mandatory parameters
	if (Kref == 0) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - please enter a valid Kref value! Current Kref = " << Kref << ". Kref cannot be 0!\n";
		HelpDP();
		return theMaterial;
	}

	if (Gref == 0) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - please enter a valid Gref value! Current Gref = " << Gref << ". Gref cannot be 0!\n";
		HelpDP();
		return theMaterial;
	}

	if (Pref < 101) {
		// Pref is kPa in most cases, though it is not a must. Just invite the user to double check the input value  
		opserr << "WARNING: OPS_DruckerPragerDMM() - please double check the Pref value in use! Current Pref = " << Pref << "!\n";
	}

	if (beVerbose) {
		opserr << "\nOPS_DruckerPragerDMM\n";
		opserr << "---------------------------------------------\n";
		opserr << "tag              : " << tag << "\n";
		opserr << "rho              : " << rho << "\n";
		opserr << "Kref             : " << Kref << "\n";
		opserr << "Gref             : " << Gref << "\n";
		opserr << "Pref             : " << Pref << "\n";
		opserr << "Modn             : " << modn << "\n";
		opserr << "TNYS             : " << TNYS << "\n";
		opserr << "cohesion         : " << cohesion << "\n";
		opserr << "frictionAngle    : " << frictionAngle << "\n";
		opserr << "dilatancyAngle   : " << dilatancyAngle << "\n";
		opserr << "peakShearStrain  : " << peakShearStrain << "\n";
		opserr << "driverType       : " << dataDriverType << "\n";
		opserr << "impl-ex          : " << integrationType << "\n";
		if (HModuli != nullptr) {
			opserr << "HModuli          : " << hmod;
		}
		if (HParams != nullptr) {
			opserr << "HParams          : " << hps;
		}
		if (HDilatancies != nullptr) {
			opserr << "HDilatancies     : " << dps;
		}
		if (isDatabase) {
			size_t totalLength = strlen(dbPathDir) + strlen(dbMainFile) + 1;
			fullPath = new char[totalLength];
			strcpy(fullPath, dbPathDir);
			strcat(fullPath, dbMainFile);
			opserr << "Path to database : " << fullPath << "\n";
			opserr << "Number of cases  : " << nFiles << "\n";
			delete[] fullPath;
			fullPath = nullptr;
		}
		opserr << "\n";
	}

	// create a nested yield surface object
	theData = std::make_shared<DataDrivenNestedSurfaces>(tag, cohesion, frictionAngle, dilatancyAngle, peakShearStrain, TNYS, dataDriverType, hmod, hps, dps, dbPathDir, dbMainFile, beVerbose);

	if (theData == nullptr) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - cannot create DruckerPragerDMM material with tag: " << tag << "\n";
		opserr << "FATAL: OPS_DruckerPragerDMM() - yield surface data yielded bad result...";
		return theMaterial;
	}

	// create a DruckerPragerDMM nDmaterial object
	theMaterial = new DruckerPragerDMM(tag, rho, Kref, Gref, Pref, modn, theData, dataDriverType, integrationType, beVerbose);

	if (theMaterial == nullptr) {
		opserr << "FATAL: OPS_DruckerPragerDMM() - cannot create DruckerPragerDMM material with tag: " << tag << "\n";
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