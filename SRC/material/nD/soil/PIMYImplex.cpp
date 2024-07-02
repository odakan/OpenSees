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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PIMYImplex.cpp,v $
// $Revision: 1.0 $
// $Date: 2024-01-05 12:00:00 $

// Written by:	    Onur Deniz Akan		(onur.akan@iusspavia.it)
// Based on:        ZHY's PressureIndependMultiYield.cpp
// Created:         December 2023
// Last Modified:
//
// Description: This file contains the implementation for the PIMYImplex function.

#include <math.h>
#include <stdlib.h>
#include "PIMYImplex.h"
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <MultiYieldSurface.h>

static const Vector& I1()
{
	// return 2nd order identity tensor 
	// basis of c is independent of the matrix representation
	static Vector c(6);
	c(0) = c(1) = c(2) = 1.0;
	return c;
}

static const Matrix& IIvol()
{
	// return 4th order volumetric operator 
	// basis of c is independent of the matrix representation
	static Matrix c(6, 6);
	c(0, 0) = 1.0;
	c(0, 1) = 1.0;
	c(0, 2) = 1.0;
	c(1, 0) = 1.0;
	c(1, 1) = 1.0;
	c(1, 2) = 1.0;
	c(2, 0) = 1.0;
	c(2, 1) = 1.0;
	c(2, 2) = 1.0;
	return c;
}

static const Matrix& IIdev()
{
	// return 4th order deviatoric operator 
	// basis of c is contravariant (i.e., operates on strain and leads to stress)
	static Matrix c(6, 6);
	c(0, 0) = 2.0 / 3.0;
	c(0, 1) = -1.0 / 3.0;
	c(0, 2) = -1.0 / 3.0;
	c(1, 0) = -1.0 / 3.0;
	c(1, 1) = 2.0 / 3.0;
	c(1, 2) = -1.0 / 3.0;
	c(2, 0) = -1.0 / 3.0;
	c(2, 1) = -1.0 / 3.0;
	c(2, 2) = 2.0 / 3.0;
	c(3, 3) = 0.5;
	c(4, 4) = 0.5;
	c(5, 5) = 0.5;
	return c;
}

Matrix PIMYImplex::theTangent(6, 6);
T2Vector PIMYImplex::subStrainRate;
int PIMYImplex::matCount = 0;
int* PIMYImplex::loadStagex = 0;  //=0 if elastic; =1 if plastic
int* PIMYImplex::ndmx = 0;        //num of dimensions (2 or 3)
int* PIMYImplex::numOfSurfacesx = 0;
bool* PIMYImplex::doImplex = false;
bool* PIMYImplex::beVerbose = false;
double* PIMYImplex::rhox = 0;
double* PIMYImplex::frictionAnglex = 0;
double* PIMYImplex::peakShearStrainx = 0;
double* PIMYImplex::refPressurex = 0;
double* PIMYImplex::cohesionx = 0;
double* PIMYImplex::pressDependCoeffx = 0;
double* PIMYImplex::residualPressx = 0;

// workspace
static Matrix workM66(6, 6);
static Vector workV6(6);

void* OPS_PIMYImplex()
{
	// display kudos
	static int kudos = 0;
	if (++kudos == 1)
		opserr << "PIMYImplex nDmaterial - Written: OD.Akan, G.Camata, E.Spacone, CG.Lai, C.Tamagnini, IUSS Pavia \n";

	// initialize pointers
	NDMaterial* theMaterial = nullptr;      // pointer to an nD material to be returned
	static double* custom_curve = nullptr;  // pointer to user input G reduction curve

	// initialize material parameters
	int tag = -1; double Gref = -1.0; double Kref = -1.0;
	double cohesion = -1.0; double gPeak = -1.0;

	// optional parameters with default values
	int maxTNYS = 80; int TNYS = 20;
	double rho = 0.0; double fAngle = 0.0;
	double Pref = 101.3; double mPow = 0.0;
	int implexFlag = 0; int verboseFlag = 0;

	// begin recieving
	int numArgs = OPS_GetNumRemainingInputArgs();
	int numData = 1;

	// check if help is requested
	if (numArgs < 3) {
		while (OPS_GetNumRemainingInputArgs() > 0) {
			const char* inputstring = OPS_GetString();
			// recive print help flag, print material command use and quit
			if (strcmp(inputstring, "-help") == 0) {
				opserr << "FATAL: OPS_PIMYImplex() - Code terminated since '-help' is the first option in the nDmaterial PIMYImplex declaration.\n";
				opserr << "                            Remove or move the '-help' option further along the declaration to continue...\n";
				exit(-1);
			}
		}
	}

	// check mandatory inputs
	if (numArgs < 6) {
		opserr << "FATAL: OPS_PIMYImplex() - please define the following minimum parameters:\n";
		opserr << "nDMaterial PIMYImplex tag? -G Gref? -K Kref? -gPeak peakShearStrain? -cohesion cohe?\n";
		return theMaterial;
	}

	// input #1 - recieve unique material tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "FATAL: OPS_PIMYImplex() - invalid tag? (must be followed by an integer)\n\n";
		return theMaterial;
	}

	// continue with the remaining inputs
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();

		// input #2 - recieve material bulk modulus at reference pressure
		if (strcmp(inputstring, "-Kref") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Kref) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid Kref value after flag: -K (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #3 - recieve material shear modulus at reference pressure
		if (strcmp(inputstring, "-Gref") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Gref) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid Gref value after flag: -G (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #4 - recieve material reference pressure
		if (strcmp(inputstring, "-Pref") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &Pref) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid Pref value after flag: -P (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #5 - recieve material mass density
		if (strcmp(inputstring, "-rho") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &rho) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid rho value after flag: -r (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #6 - recieve material elastic modulus update power
		if (strcmp(inputstring, "-mPow") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &mPow) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid pressDependCoe value after flag: -mPow (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #7 - recieve total number of yield surfaces
		if (strcmp(inputstring, "-tnys") == 0 || strcmp(inputstring, "-TNYS") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, &TNYS) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid tnys value after flag: -tnys (must be followed by an integer)\n";
				return theMaterial;
			}
			if (abs(TNYS) > maxTNYS || TNYS == 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid tnys value! (must be an integer tnys <= " << maxTNYS << " and tnys cannot be 0)\n";
				return theMaterial;
			}
			if (TNYS < 0) {
				// recieve user input yield surfaces
				TNYS = abs(TNYS);	// make sure TNYS is a positive integer
				custom_curve = new double[int(2 * TNYS)];
				for (int i = 0; i < int(2 * TNYS); i++) {
					if (OPS_GetDoubleInput(&numData, &custom_curve[i]) < 0) {
						opserr << "WARNING invalid " << " double" << "\n";
						opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
						return 0;
					}
				}
			}

		}

		// input #8 - recieve yield surface cohesion
		if (strcmp(inputstring, "-cohesion") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &cohesion) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid cohesion value after flag: -cohesion (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #9 - recieve yield surface friction angle
		if (strcmp(inputstring, "-fAngle") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &fAngle) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid frictionAngle value after flag: -fAngle (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #11 - recieve yield surface peak shear strain
		if (strcmp(inputstring, "-gPeak") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &gPeak) < 0) {
				opserr << "FATAL: OPS_PIMYImplex() - invalid peakShearStrain value after flag: -gPeak (must be followed by a double)\n";
				return theMaterial;
			}
		}

		// input #12 - recieve material integration type
		if (strcmp(inputstring, "-implex") == 0) {
			implexFlag = 1; // rephrased-flow formulation
		}

		// input #26 - recieve verbosity paramter
		if (strcmp(inputstring, "-info") == 0) {
			verboseFlag = 1;
		}
	}

	// do a value check
	if (tag < 0) { opserr << "FATAL: OPS_PIMYImplex() - tag < 0!\n"; return theMaterial; }
	if (rho < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - rho < 0!\n"; return theMaterial; }
	if (Gref < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - Gref < 0!\n"; return theMaterial; }
	if (Kref < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - Kref < 0!\n"; return theMaterial; }
	if (Pref < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - Pref < 0!\n"; return theMaterial; }
	if (mPow < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - mPow < 0!\n"; return theMaterial; }
	if (fAngle < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - fAngle < 0!\n"; return theMaterial; }
	if (gPeak < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - gPeak < 0!\n"; return theMaterial; }
	if (TNYS < 0) { opserr << "FATAL: OPS_PIMYImplex() - TNYS < 0!\n"; return theMaterial; }
	if (cohesion < 0.0) { opserr << "FATAL: OPS_PIMYImplex() - cohesion < 0!\n"; return theMaterial; }
	if (implexFlag < 0) { opserr << "FATAL: OPS_PIMYImplex() - integrationType < 0!\n"; return theMaterial; }

	if (verboseFlag == 1) {
		opserr << "\nOPS_PIMYImplex\n";
		opserr << "------------------------\n";
		opserr << "tag      : " << tag << "\n";
		opserr << "rho      : " << rho << "\n";
		opserr << "Gref     : " << Gref << "\n";
		opserr << "Kref     : " << Kref << "\n";
		opserr << "cohesion : " << cohesion << "\n";
		opserr << "gPeak    : " << gPeak << "\n";
		opserr << "fAngle   : " << fAngle << "\n";
		opserr << "Pref     : " << Pref << "\n";
		opserr << "mPow     : " << mPow << "\n";
		opserr << "TNYS     : " << TNYS << "\n";
		if (implexFlag == 1) {
			opserr << "solution : IMPL-EX\n";
		}
		else {
			opserr << "solution : implicit\n";
		}
		if (custom_curve != nullptr) {
			opserr << "G/Gmax   : custom curve\n";
		}
		else {
			opserr << "G/Gmax   : hyperbolic\n";
		}
		opserr << "\n";
	}

	// create a PIMYImplex nDmaterial object
	theMaterial = new PIMYImplex(tag, rho, Gref, Kref, cohesion, gPeak, fAngle, Pref,
		mPow, TNYS, custom_curve, implexFlag, verboseFlag);

	if (theMaterial == nullptr) {
		opserr << "FATAL: OPS_PIMYImplex() - cannot create PIMYImplex material with tag: " << tag << "\n";
		return theMaterial;
	}

	// free the memory
	if (custom_curve != nullptr) {
		delete[] custom_curve;
		custom_curve = nullptr;
	}

	return theMaterial;
}

PIMYImplex::PIMYImplex(int tag,
	double r, double refShearModul,
	double refBulkModul,
	double cohesi, double peakShearStra,
	double frictionAng, double refPress, double pressDependCoe,
	int numberOfYieldSurf, double* gredu,
	int implexFlag, int verboseFlag)
	: NDMaterial(tag, ND_TAG_PIMYImplex), currentStress(), currentStressImplex(),
	trialStress(), currentStrain(), strainRate()
{
	// handle solution options
		// handle dimension
	int nd = 0;
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nd = 2;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nd = 3;
	}
	else {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - unknown model dimension...\n";
		exit(-1);
	}

	if (nd != 2 && nd != 3) {
		opserr << "FATAL:PIMYImplex:: dimension error" << endln;
		opserr << "Dimension has to be 2 or 3, you give nd= " << nd << endln;
		exit(-1);
	}
	if (refShearModul <= 0) {
		opserr << "FATAL:PIMYImplex::PIMYImplex: refShearModulus <= 0" << endln;
		exit(-1);
	}
	if (refBulkModul <= 0) {
		opserr << "FATAL:PIMYImplex::PIMYImplex: refBulkModulus <= 0" << endln;
		exit(-1);
	}
	if (frictionAng < 0.) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: frictionAngle < 0" << endln;
		opserr << "Will reset frictionAngle to zero." << endln;
		frictionAng = 0.;
	}
	if (frictionAng == 0. && cohesi <= 0.) {
		opserr << "FATAL:PIMYImplex::PIMYImplex: frictionAngle && cohesion <= 0." << endln;
		exit(-1);
	}
	if (cohesi <= 0) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: cohesion <= 0" << endln;
		opserr << "Will reset cohesion to zero." << endln;
		cohesi = 0.;
	}
	if (peakShearStra <= 0) {
		opserr << "FATAL:PIMYImplex::PIMYImplex: peakShearStra <= 0" << endln;
		exit(-1);
	}
	if (refPress <= 0) {
		opserr << "FATAL:PIMYImplex::PIMYImplex: refPress <= 0" << endln;
		exit(-1);
	}
	if (pressDependCoe < 0) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: pressDependCoe < 0" << endln;
		opserr << "Will reset pressDependCoe to zero." << endln;
		pressDependCoe = 0.;
	}
	if (pressDependCoe > 0 && frictionAng == 0) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: pressDependCoe > 0 while frictionAngle = 0" << endln;
		opserr << "Will reset pressDependCoe to zero." << endln;
		pressDependCoe = 0.;
	}
	if (numberOfYieldSurf <= 0) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: numberOfSurfaces <= 0" << endln;
		opserr << "Will use 10 yield surfaces." << endln;
		numberOfYieldSurf = 10;
	}
	if (numberOfYieldSurf > 100) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: numberOfSurfaces > 100" << endln;
		opserr << "Will use 100 yield surfaces." << endln;
		numberOfYieldSurf = 100;
	}
	if (r < 0) {
		opserr << "WARNING:PIMYImplex::PIMYImplex: mass density < 0" << endln;
		opserr << "Will use rho = 0." << endln;
		r = 0.;
	}

	int* temp1 = loadStagex;
	int* temp2 = ndmx;
	double* temp3 = rhox;
	double* temp6 = frictionAnglex;
	double* temp7 = peakShearStrainx;
	double* temp8 = refPressurex;
	double* temp9 = cohesionx;
	double* temp10 = pressDependCoeffx;
	int* temp11 = numOfSurfacesx;
	double* temp12 = residualPressx;
	bool* temp13 = doImplex;
	bool* temp15 = beVerbose;

	int newCount = matCount + 1;
	loadStagex = new int[newCount];
	ndmx = new int[newCount];
	rhox = new double[newCount];
	frictionAnglex = new double[newCount];
	peakShearStrainx = new double[newCount];
	refPressurex = new double[newCount];
	cohesionx = new double[newCount];
	pressDependCoeffx = new double[newCount];
	numOfSurfacesx = new int[newCount];
	residualPressx = new double[newCount];
	doImplex = new bool[newCount];
	beVerbose = new bool[newCount];

	for (int i = 0; i < matCount; i++) {
		loadStagex[i] = temp1[i];
		ndmx[i] = temp2[i];
		rhox[i] = temp3[i];
		frictionAnglex[i] = temp6[i];
		peakShearStrainx[i] = temp7[i];
		refPressurex[i] = temp8[i];
		cohesionx[i] = temp9[i];
		pressDependCoeffx[i] = temp10[i];
		numOfSurfacesx[i] = temp11[i];
		residualPressx[i] = temp12[i];
		doImplex[i] = temp13[i];
		beVerbose[i] = temp15[i];
	}

	if (matCount > 0) {
		delete[] temp1; delete[] temp2; delete[] temp3;
		delete[] temp6; delete[] temp7; delete[] temp8;
		delete[] temp9; delete[] temp10; delete[] temp11;
		delete[] temp12; delete[] temp13;
		delete[] temp15;
	}

	ndmx[matCount] = nd;
	loadStagex[matCount] = 0;   //default
	refShearModulus = refShearModul;
	refBulkModulus = refBulkModul;
	frictionAnglex[matCount] = frictionAng;
	peakShearStrainx[matCount] = peakShearStra;
	refPressurex[matCount] = -refPress;  //compression is negative
	cohesionx[matCount] = cohesi;
	pressDependCoeffx[matCount] = pressDependCoe;
	numOfSurfacesx[matCount] = numberOfYieldSurf;
	rhox[matCount] = r;
	doImplex[matCount] = (implexFlag > 0);
	beVerbose[matCount] = (verboseFlag == 1);

	e2p = 0;
	matN = matCount;
	matCount = newCount;

	theSurfaces = new MultiYieldSurface[numberOfYieldSurf + 1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numberOfYieldSurf + 1];
	activeSurfaceNum = committedActiveSurf = 0;

	mGredu = gredu;
	setUpSurfaces(gredu);  // residualPress is calculated inside.
}


PIMYImplex::PIMYImplex()
	: NDMaterial(0, ND_TAG_PIMYImplex),
	currentStress(), trialStress(), currentStrain(), currentStressImplex(),
	strainRate(), theSurfaces(0), committedSurfaces(0)
{
	//does nothing
}


PIMYImplex::PIMYImplex(const PIMYImplex& a)
	: NDMaterial(a.getTag(), ND_TAG_PIMYImplex),
	currentStress(a.currentStress), trialStress(a.trialStress),
	currentStressImplex(a.currentStressImplex),
	currentStrain(a.currentStrain), strainRate(a.strainRate)
{
	matN = a.matN;
	e2p = a.e2p;
	refShearModulus = a.refShearModulus;
	refBulkModulus = a.refBulkModulus;

	int numOfSurfaces = numOfSurfacesx[matN];

	committedActiveSurf = a.committedActiveSurf;
	activeSurfaceNum = a.activeSurfaceNum;

	theSurfaces = new MultiYieldSurface[numOfSurfaces + 1];  //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces + 1];
	for (int i = 1; i <= numOfSurfaces; i++) {
		committedSurfaces[i] = a.committedSurfaces[i];
		theSurfaces[i] = a.theSurfaces[i];
	}
}


PIMYImplex::~PIMYImplex()
{
	if (theSurfaces != 0) delete[] theSurfaces;
	if (committedSurfaces != 0) delete[] committedSurfaces;
}


void PIMYImplex::elast2Plast(void)
{
	int loadStage = loadStagex[matN];
	double frictionAngle = frictionAnglex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	if (loadStage != 1 || e2p == 1) return;
	e2p = 1;

	if (currentStress.volume() > 0. && frictionAngle > 0.) {
		//opserr << "WARNING:PIMYImplex::elast2Plast(): material in tension." << endln;
		currentStress.setData(currentStress.deviator(), 0);
	}

	paramScaling();  // scale surface parameters corresponding to initial confinement

	// Active surface is 0, return
	if (currentStress.deviatorLength() == 0.) return;

	// Find active surface
	while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
		if (committedActiveSurf == numOfSurfaces) {
			//opserr <<"WARNING:PIMYImplex::elast2Plast(): stress out of failure surface"<<endln;
			deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
			initSurfaceUpdate();
			return;
		}
	}
	committedActiveSurf--;
	initSurfaceUpdate();
}


int PIMYImplex::setTrialStrain(const Vector& strain)
{
	//opserr << "setTrialStrain() - in\n";
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	// handle implex time step
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!dtime_first_set) {
			dtime_n_commit = dtime_n;
			dtime_first_set = true;
		}
	}

	static Vector temp(6);
	if (ndm == 3 && strain.Size() == 6)
		temp = strain;
	else if (ndm == 2 && strain.Size() == 3) {
		temp[0] = strain[0];
		temp[1] = strain[1];
		temp[2] = 0.0;
		temp[3] = strain[2];
		temp[4] = 0.0;
		temp[5] = 0.0;
	}
	else {
		opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);
	}

	//strainRate.setData(temp-currentStrain.t2Vector(1),1);
	temp -= currentStrain.t2Vector(1);
	strainRate.setData(temp, 1);
	//opserr << "Strain increment = " << strainRate.t2Vector();
	return 0;
}


int PIMYImplex::setTrialStrain(const Vector& strain, const Vector& rate)
{
	return setTrialStrain(strain);
}


int PIMYImplex::setTrialStrainIncr(const Vector& strain)
{
	//opserr << "setTrialStrainIncr() - in\n";
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	// handle implex time step
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!dtime_first_set) {
			dtime_n_commit = dtime_n;
			dtime_first_set = true;
		}
	}

	static Vector temp(6);
	if (ndm == 3 && strain.Size() == 6)
		temp = strain;
	else if (ndm == 2 && strain.Size() == 3) {
		temp[0] = strain[0];
		temp[1] = strain[1];
		temp[3] = strain[2];
	}
	else {
		opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);
	}

	strainRate.setData(temp, 1);
	return 0;
}


int PIMYImplex::setTrialStrainIncr(const Vector& strain, const Vector& rate)
{
	return setTrialStrainIncr(strain);
}


const Matrix& PIMYImplex::getTangent(void)
{
	int loadStage = loadStagex[matN];
	bool implex = doImplex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 3;

	if (loadStage == 1 && e2p == 0) elast2Plast();

	if (loadStage != 1) {  //linear elastic
		theTangent = 2 * refShearModulus * IIdev() + refBulkModulus * IIvol();
	}
	else {
		if (!implex) {
			double coeff;
			static Vector devia(6);

			if (activeSurfaceNum > 0) {
				devia = trialStress.deviator();
				devia -= theSurfaces[activeSurfaceNum].center();

				double size = theSurfaces[activeSurfaceNum].size();
				double plastModul = theSurfaces[activeSurfaceNum].modulus();
				coeff = 6. * refShearModulus * refShearModulus / (2. * refShearModulus + plastModul) / size / size;
			}

			else coeff = 0.;

			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++) {
					theTangent(i, j) = -coeff * devia[i] * devia[j];
					if (i == j) theTangent(i, j) += refShearModulus;
					if (i < 3 && j < 3 && i == j) theTangent(i, j) += refShearModulus;
					if (i < 3 && j < 3) theTangent(i, j) += (refBulkModulus - 2. * refShearModulus / 3.);
				}
		}
		else {
			// time factor for explicit extrapolation
			double time_factor = 1.0;
			if (dtime_n_commit > 0.0)
				time_factor = dtime_n / dtime_n_commit;
			// note: the implex method just wants the ratio of the new to the old time step
			// not the real time step, so it is just fine to assume it to 1.
			// otherwise we have to deal with the problem of the opensees pseudo-time step
			// being the load multiplier in continuation methods...
			time_factor = 0.0;

			// extrapolate state variables
			lambda = lambda_commit + time_factor * (lambda_commit - lambda_commit_old);
			chi = chi_commit; //+ time_factor * (chi_commit - chi_commit_old);
			ksi = ksi_commit; //+ time_factor * (ksi_commit - ksi_commit_old);
			double lambda_bar = (chi > 0.0) ? (lambda / chi) : (0.0);

			// compute consistent tangent
			theTangent.addMatrix(0.0, IIdev(), 2.0 * refShearModulus * (1.0 - lambda_bar * ksi * 2.0 * refShearModulus));
			theTangent.addMatrix(1.0, IIvol(), refBulkModulus);
		}
	}

	if (ndm == 3)
		return theTangent;
	else {
		static Matrix workM(3, 3);
		workM(0, 0) = theTangent(0, 0);
		workM(0, 1) = theTangent(0, 1);
		workM(0, 2) = theTangent(0, 3);
		workM(1, 0) = theTangent(1, 0);
		workM(1, 1) = theTangent(1, 1);
		workM(1, 2) = theTangent(1, 3);
		workM(2, 0) = theTangent(3, 0);
		workM(2, 1) = theTangent(3, 1);
		workM(2, 2) = theTangent(3, 3);
		return workM;
	}
}


const Matrix& PIMYImplex::getInitialTangent(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 3;

	theTangent = 2.0 * refShearModulus * IIdev() + refBulkModulus * IIvol();

	if (ndm == 3)
		return theTangent;
	else {
		static Matrix workM(3, 3);
		workM(0, 0) = theTangent(0, 0);
		workM(0, 1) = theTangent(0, 1);
		workM(0, 2) = theTangent(0, 3);
		workM(1, 0) = theTangent(1, 0);
		workM(1, 1) = theTangent(1, 1);
		workM(1, 2) = theTangent(1, 3);
		workM(2, 0) = theTangent(3, 0);
		workM(2, 1) = theTangent(3, 1);
		workM(2, 2) = theTangent(3, 3);
		return workM;
	}
}


const Vector& PIMYImplex::getStress(void)
{
	int loadStage = loadStagex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	bool implex = doImplex[matN];
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 3;

	if (loadStage == 1 && e2p == 0) elast2Plast();

	if (loadStage != 1) {  //linear elastic
		theTangent = 2.0 * refShearModulus * IIdev() + refBulkModulus * IIvol();
		workV6.addVector(0.0, currentStress.t2Vector(), 1.0);
		workV6.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
		trialStress.setData(workV6);
	}
	else {
		if (!implex) {
			implicitSress();
		}
		else {
			// time factor for explicit extrapolation
			double time_factor = 1.0;
			if (dtime_n_commit > 0.0)
				time_factor = dtime_n / dtime_n_commit;
			// note: the implex method just wants the ratio of the new to the old time step
			// not the real time step, so it is just fine to assume it to 1.
			// otherwise we have to deal with the problem of the opensees pseudo-time step
			// being the load multiplier in continuation methods...
			time_factor = 0.0;

			// extrapolate state variables
			lambda = lambda_commit + time_factor * (lambda_commit - lambda_commit_old);
			chi = chi_commit;// +time_factor * (chi_commit - chi_commit_old);
			ksi = ksi_commit;// +time_factor * (ksi_commit - ksi_commit_old);
			double lambda_bar = (chi > 0.0) ? (lambda / chi) : (0.0);

			// compute implex stress
				// trial stress
			theTangent = 2.0 * refShearModulus * IIdev() + refBulkModulus * IIvol();
			workV6.addVector(0.0, currentStress.t2Vector(), 1.0);
			workV6.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
			trialStress.setData(workV6);

			// compute plastic strain	
			workV6.addVector(0.0, theSurfaces[activeSurfaceNum].center(), 1.0);		// get the center
			workV6.addVector(-1.0, trialStress.deviator(), 1.0);					// get the shifted deviator
			workV6 *= lambda_bar * ksi;												// compute the contact stress and then the plastic strain

			// correct stress
			workV6.addVector(-2.0 * refShearModulus, trialStress.deviator(), 1.0);	// correct the stress deviator
			trialStress.setData(workV6, trialStress.volume());						// update trial stress
		}
	}

	if (ndm == 3)
		return trialStress.t2Vector();
	else {
		static Vector workV(3);
		workV[0] = trialStress.t2Vector()[0];
		workV[1] = trialStress.t2Vector()[1];
		workV[2] = trialStress.t2Vector()[3];
		return workV;
	}
}


const Vector& PIMYImplex::getStrain(void)
{
	return getCommittedStrain();
}


int PIMYImplex::commitState(void)
{
	int loadStage = loadStagex[matN];
	bool implex = doImplex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	if (loadStage == 1) {

		if (implex) {
			// do an implicit iteration before commit
			currentStressImplex = trialStress;		// store for output
			implicitSress();						// do implicit correction
		}

		// commit implex variables
		chi_commit_old = chi_commit;
		chi_commit = chi;
		lambda_commit_old = lambda_commit;
		lambda_commit = lambda;
		ksi_commit_old = ksi_commit;
		ksi_commit = ksi;

		committedActiveSurf = activeSurfaceNum;
		for (int i = 1; i <= numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
	}

	currentStress = trialStress;

	//currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
	static Vector temp(6);
	temp = currentStrain.t2Vector();
	temp += strainRate.t2Vector();
	currentStrain.setData(temp);
	temp.Zero();
	strainRate.setData(temp);

	return 0;
}


int PIMYImplex::revertToLastCommit(void)
{
	int loadStage = loadStagex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	if (loadStage == 1) {
		// revert implex variables
		chi_commit = chi_commit_old;
		chi = chi_commit;
		lambda_commit = lambda_commit_old;
		lambda = lambda_commit;
		ksi_commit = ksi_commit_old;
		ksi = ksi_commit;

		activeSurfaceNum = committedActiveSurf;
		for (int i = 1; i <= numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
	}

	trialStress = currentStress;
	//currentStrain is always the lastest comitted (step n)

	return 0;
}


NDMaterial* PIMYImplex::getCopy(void)
{
	PIMYImplex* copy = new PIMYImplex(*this);
	return copy;
}


NDMaterial* PIMYImplex::getCopy(const char* code)
{
	if (strcmp(code, "PlaneStrain") == 0 || strcmp(code, "ThreeDimensional") == 0) {
		PIMYImplex* copy = new PIMYImplex(*this);
		return copy;
	}
	else {
		opserr << "ERROR PIMYImplex::getCopy -- cannot make copy for type " << code << endln;
		return 0;
	}
}


const char* PIMYImplex::getType(void) const
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int PIMYImplex::getOrder(void) const
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	return (ndm == 2) ? 3 : 6;
}


int PIMYImplex::setParameter(const char** argv, int argc, Parameter& param)
{
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	// check for material tag
	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0], "shearModulus") == 0) {
			return param.addObject(10, this);
		}
		else if (strcmp(argv[0], "bulkModulus") == 0) {
			return param.addObject(11, this);
		}
		else if (strcmp(argv[0], "frictionAngle") == 0) {
			return param.addObject(12, this);
		}
		else if (strcmp(argv[0], "cohesion") == 0) {
			return param.addObject(13, this);
		}
	}
	return -1;
}

int PIMYImplex::updateParameter(int responseID, Information& info)
{
	if (responseID == 1) {
		static int switch2Implex = 0;
		static int switch2Implicit = 0;
		static int msgtag = 0;

		// reset message if the material has changed
		if (msgtag != this->getTag()) {
			switch2Implex = 0;
			switch2Implicit = 0;
			msgtag = this->getTag();
		}

		if (info.theInt == 11) {
			loadStagex[matN] = 1;
			doImplex[matN] = true;
			// display message
			switch2Implicit = 0;
			if (++switch2Implex == 1)
				if (beVerbose[matN]) { opserr << "PIMYImplex::updateParameter() - material " << msgtag << ": switching to IMPL-EX solution...\n"; }
		}
		else if (info.theInt == 10) {
			loadStagex[matN] = 1;
			doImplex[matN] = false;
			// display message
			switch2Implex = 0;
			if (++switch2Implicit == 1)
				if (beVerbose[matN]) { opserr << "PIMYImplex::updateParameter() - material " << msgtag << ": switching to IMPLICIT solution...\n"; }
		}
		else {
			loadStagex[matN] = info.theInt;
		}
	}
	else if (responseID == 10) {
		refShearModulus = info.theDouble;
	}
	else if (responseID == 11) {
		refBulkModulus = info.theDouble;
	}
	else if (responseID == 12) {
		frictionAnglex[matN] = info.theDouble;
		double* g = 0;
		setUpSurfaces(g);
		paramScaling();
		initSurfaceUpdate();
	}
	else if (responseID == 13) {
		cohesionx[matN] = info.theDouble;
		double* g = 0;
		setUpSurfaces(g);
		paramScaling();
		initSurfaceUpdate();
	}

	// used by BBarFourNodeQuadUP element
	else if (responseID == 20 && ndmx[matN] == 2)
		ndmx[matN] = 0;

	return 0;
}


int PIMYImplex::sendSelf(int commitTag, Channel& theChannel)
{
	int loadStage = loadStagex[matN];
	int ndm = ndmx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double rho = rhox[matN];
	double frictionAngle = frictionAnglex[matN];
	double peakShearStrain = peakShearStrainx[matN];
	double refPressure = refPressurex[matN];
	double cohesion = cohesionx[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	double residualPress = residualPressx[matN];
	bool implex = doImplex[matN];
	bool info = beVerbose[matN];

	int i, res = 0;

	static ID idData(8);
	idData(0) = this->getTag();
	idData(1) = numOfSurfaces;
	idData(2) = loadStage;
	idData(3) = ndm;
	idData(4) = matN;
	idData(5) = matCount;
	idData(6) = (implex) ? (1.0) : (0.0);
	idData(7) = (info) ? (1.0) : (0.0);

	res += theChannel.sendID(this->getDbTag(), commitTag, idData);
	if (res < 0) {
		opserr << "PIMYImplex::sendSelf -- could not send ID\n";
		return res;
	}

	Vector data(34 + numOfSurfaces * 8);
	static Vector temp(6);
	data(0) = rho;
	data(1) = refShearModulus;
	data(2) = refBulkModulus;
	data(3) = frictionAngle;
	data(4) = peakShearStrain;
	data(5) = refPressure;
	data(6) = cohesion;
	data(7) = pressDependCoeff;
	data(8) = residualPress;
	data(9) = e2p;
	data(10) = committedActiveSurf;
	data(11) = activeSurfaceNum;

	temp = currentStress.t2Vector();
	for (i = 0; i < 6; i++) data(i + 12) = temp[i];

	temp = currentStrain.t2Vector();
	for (i = 0; i < 6; i++) data(i + 18) = temp[i];

	data(24) = lambda;
	data(25) = lambda_commit;
	data(26) = lambda_commit_old;
	data(27) = chi;
	data(28) = chi_commit;
	data(29) = chi_commit_old;
	data(30) = ksi;
	data(31) = ksi_commit;
	data(32) = ksi_commit_old;
	data(33) = dtime_n;
	data(34) = dtime_n_commit;
	data(35) = (dtime_is_user_defined) ? (1.0) : (0.0);
	data(36) = (dtime_first_set) ? (1.0) : (0.0);

	for (i = 0; i < numOfSurfaces; i++) {
		int k = 37 + i * 8;
		data(k) = committedSurfaces[i + 1].size();
		data(k + 1) = committedSurfaces[i + 1].modulus();
		temp = committedSurfaces[i + 1].center();
		data(k + 2) = temp(0);
		data(k + 3) = temp(1);
		data(k + 4) = temp(2);
		data(k + 5) = temp(3);
		data(k + 6) = temp(4);
		data(k + 7) = temp(5);
	}

	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "PIMYImplex::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}


int PIMYImplex::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int i, res = 0;

	static ID idData(8);

	res += theChannel.recvID(this->getDbTag(), commitTag, idData);
	if (res < 0) {
		opserr << "PIMYImplex::recvSelf -- could not recv ID\n";
		return res;
	}

	this->setTag((int)idData(0));
	int numOfSurfaces = idData(1);
	int loadStage = idData(2);
	int ndm = idData(3);
	matN = idData(4);

	int matCountSendSide = idData(5);

	bool implex = (idData(6) == 1.0);
	bool info = (idData(7) == 1.0);

	Vector data(34 + idData(1) * 8);
	static Vector temp(6);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "PIMYImplex::recvSelf -- could not recv Vector\n";
		return res;
	}

	double rho = data(0);
	refShearModulus = data(1);
	refBulkModulus = data(2);
	double frictionAngle = data(3);
	double peakShearStrain = data(4);
	double refPressure = data(5);
	double cohesion = data(6);
	double pressDependCoeff = data(7);
	double residualPress = data(8);
	e2p = data(9);
	committedActiveSurf = data(10);
	activeSurfaceNum = data(11);

	for (i = 0; i < 6; i++)
		temp[i] = data(i + 12);
	currentStress.setData(temp);

	for (i = 0; i < 6; i++)
		temp[i] = data(i + 18);
	currentStrain.setData(temp);

	lambda = data(24);
	lambda_commit = data(25);
	lambda_commit_old = data(26);
	chi = data(27);
	chi_commit = data(28);
	chi_commit_old = data(29);
	ksi = data(30);
	ksi_commit = data(31);
	ksi_commit_old = data(32);
	dtime_n = data(33);
	dtime_n_commit = data(34);
	dtime_is_user_defined = (data(35) == 1.0);
	dtime_first_set = (data(36) == 1.0);

	if (committedSurfaces != 0) {
		delete[] committedSurfaces;
		delete[] theSurfaces;
	}

	theSurfaces = new MultiYieldSurface[numOfSurfaces + 1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces + 1];

	for (i = 0; i < numOfSurfaces; i++) {
		int k = 37 + i * 8;
		temp(0) = data(k + 2);
		temp(1) = data(k + 3);
		temp(2) = data(k + 4);
		temp(3) = data(k + 5);
		temp(4) = data(k + 6);
		temp(5) = data(k + 7);
		committedSurfaces[i + 1].setData(temp, data(k), data(k + 1));
	}

	int* temp1, * temp2, * temp11;
	double* temp3, * temp6, * temp7, * temp8, * temp9, * temp10, * temp12;
	bool* temp13, * temp14, * temp15;

	if (matCountSendSide > matCount) {

		temp1 = loadStagex;
		temp2 = ndmx;
		temp3 = rhox;
		temp6 = frictionAnglex;
		temp7 = peakShearStrainx;
		temp8 = refPressurex;
		temp9 = cohesionx;
		temp10 = pressDependCoeffx;
		temp11 = numOfSurfacesx;
		temp12 = residualPressx;
		temp13 = doImplex;
		temp15 = beVerbose;

		loadStagex = new int[matCountSendSide];
		ndmx = new int[matCountSendSide];
		rhox = new double[matCountSendSide];
		frictionAnglex = new double[matCountSendSide];
		peakShearStrainx = new double[matCountSendSide];
		refPressurex = new double[matCountSendSide];
		cohesionx = new double[matCountSendSide];
		pressDependCoeffx = new double[matCountSendSide];
		numOfSurfacesx = new int[matCountSendSide];
		residualPressx = new double[matCountSendSide];
		doImplex = new bool[matCountSendSide];
		beVerbose = new bool[matCountSendSide];

		for (int i = 0; i < matCount; i++) {
			loadStagex[i] = temp1[i];
			ndmx[i] = temp2[i];
			rhox[i] = temp3[i];
			frictionAnglex[i] = temp6[i];
			peakShearStrainx[i] = temp7[i];
			refPressurex[i] = temp8[i];
			cohesionx[i] = temp9[i];
			pressDependCoeffx[i] = temp10[i];
			numOfSurfacesx[i] = temp11[i];
			residualPressx[i] = temp12[i];
			doImplex[i] = temp13[i];
			beVerbose[i] = temp15[i];
		}
		if (matCount > 0) {
			delete[] temp1; delete[] temp2; delete[] temp3;
			delete[] temp6; delete[] temp7; delete[] temp8;
			delete[] temp9; delete[] temp10; delete[] temp11;
			delete[] temp12; delete[] temp13; delete[] temp15;
		}
		matCount = matCountSendSide;
	}

	loadStagex[matN] = loadStage;
	ndmx[matN] = ndm;
	numOfSurfacesx[matN] = numOfSurfaces;
	rhox[matN] = rho;
	frictionAnglex[matN] = frictionAngle;
	peakShearStrainx[matN] = peakShearStrain;
	refPressurex[matN] = refPressure;
	cohesionx[matN] = cohesion;
	pressDependCoeffx[matN] = pressDependCoeff;
	residualPressx[matN] = residualPress;
	doImplex[matN] = implex;
	beVerbose[matN] = info;

	return res;
}


Response*
PIMYImplex::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// begin change by Alborz Ghofrani - UW --- get only 6 components of stress
	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		if ((argc > 1) && (atoi(argv[1]) > 2) && (atoi(argv[1]) < 8)) {
			return new MaterialResponse(this, 2 + atoi(argv[1]), this->getStressToRecord(atoi(argv[1])));
		}
		else {
			return new MaterialResponse(this, 1, this->getCommittedStress());
		}
	// end change by Alborz Ghofrani - UW

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getCommittedStrain());

	else if (strcmp(argv[0], "tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());

	else if (strcmp(argv[0], "backbone") == 0) {
		int numOfSurfaces = numOfSurfacesx[matN];
		static Matrix curv(numOfSurfaces + 1, (argc - 1) * 2);
		for (int i = 1; i < argc; i++)
			curv(0, (i - 1) * 2) = atoi(argv[i]);
		return new MaterialResponse(this, 4, curv);
	}
	else if (strcmp(argv[0], "stressImplex") == 0 || strcmp(argv[0], "stressesImplex") == 0) {
		if ((argc > 1) && (atoi(argv[1]) > 2) && (atoi(argv[1]) < 8)) {
			return new MaterialResponse(this, 19 + atoi(argv[1]), this->getStressToRecord(atoi(argv[1]), true));
		}
		else {
			return new MaterialResponse(this, 21, this->getCommittedStress(true));
		}
	}
	else if (strcmp(argv[0], "stateVarsImplex") == 0) {
		return new MaterialResponse(this, 27, this->getImplexSateVariables());
	}
	else
		return 0;
}


int PIMYImplex::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getCommittedStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getCommittedStrain();
		return 0;
	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = getTangent();
		return 0;
	case 4:
		if (matInfo.theMatrix != 0)
			getBackbone(*(matInfo.theMatrix));
		return 0;
		// begin change by Alborz Ghofrani UW --- get 6 components of stress
	case 5:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(3);
		return 0;
	case 6:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(4);
		return 0;
	case 7:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(5);
		return 0;
	case 8:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(6);
		return 0;
	case 9:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(7);
		return 0;
		// end change by Alborz Ghofrani UW
	case 21:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getCommittedStress(true);
		return 0;
	case 22:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(3, true);
		return 0;
	case 23:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(4, true);
		return 0;
	case 24:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(5, true);
		return 0;
	case 25:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(6, true);
		return 0;
	case 26:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStressToRecord(7);
		return 0;
	case 27:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getImplexSateVariables();
		return 0;
	default:
		return -1;
	}
}


void PIMYImplex::getBackbone(Matrix& bb)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	double vol, conHeig, scale, factor, shearModulus, stress1,
		stress2, strain1, strain2, plastModulus, elast_plast, gre;

	for (int k = 0; k < bb.noCols() / 2; k++) {
		vol = bb(0, k * 2);
		if (vol <= 0.) {
			opserr << k << "\nNDMaterial " << this->getTag()
				<< ": invalid confinement for backbone recorder, " << vol << endln;
			continue;
		}
		conHeig = vol + residualPress;
		scale = -conHeig / (refPressure - residualPress);
		factor = pow(scale, pressDependCoeff);
		shearModulus = factor * refShearModulus;
		for (int i = 1; i <= numOfSurfaces; i++) {
			if (i == 1) {
				stress2 = committedSurfaces[i].size() * factor / sqrt(3.0);
				strain2 = stress2 / shearModulus;
				bb(1, k * 2) = strain2; bb(1, k * 2 + 1) = shearModulus;
			}
			else {
				stress1 = stress2; strain1 = strain2;
				plastModulus = factor * committedSurfaces[i - 1].modulus();
				elast_plast = 2 * shearModulus * plastModulus / (2 * shearModulus + plastModulus);
				stress2 = factor * committedSurfaces[i].size() / sqrt(3.0);
				strain2 = 2 * (stress2 - stress1) / elast_plast + strain1;
				gre = stress2 / strain2;
				bb(i, k * 2) = strain2; bb(i, k * 2 + 1) = gre;
			}
		}
	}

}

void PIMYImplex::Print(OPS_Stream& s, int flag)
{
	s << "PIMYImplex - loadStage: " << loadStagex[matN] << endln;
}


const Vector& PIMYImplex::getCommittedStress(bool isImplex)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;
	int numOfSurfaces = numOfSurfacesx[matN];
	static Vector temp6(6);
	T2Vector& stress = currentStress;
	T2Vector& stressImplex = currentStressImplex;

	if (isImplex)
		temp6 = stressImplex.t2Vector();
	else
		temp6 = stress.t2Vector();

	double scale = sqrt(3. / 2.) * currentStress.deviatorLength() / committedSurfaces[numOfSurfaces].size();
	if (loadStagex[matN] != 1) scale = 0.;
	if (ndm == 3) {
		static Vector temp7(7);
		temp7[0] = temp6[0];
		temp7[1] = temp6[1];
		temp7[2] = temp6[2];
		temp7[3] = temp6[3];
		temp7[4] = temp6[4];
		temp7[5] = temp6[5];
		temp7[6] = scale;
		return temp7;
	}
	else {
		static Vector temp5(5);
		temp5[0] = temp6[0];
		temp5[1] = temp6[1];
		temp5[2] = temp6[2];
		temp5[3] = temp6[3];
		temp5[4] = scale;
		return temp5;
	}
}

// begin change by Alborz Ghofrani - UW --- get 6 components of stress
const Vector& PIMYImplex::getStressToRecord(int numOutput, bool isImplex)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3) {
		static Vector temp7(7);
		temp7 = this->getCommittedStress(isImplex);
		if (numOutput == 6)
		{
			static Vector temp6(6);
			temp6[0] = temp7[0];
			temp6[1] = temp7[1];
			temp6[2] = temp7[2];
			temp6[3] = temp7[3];
			temp6[4] = temp7[4];
			temp6[5] = temp7[5];
			return temp6;
		}
		else if (numOutput == 7)
		{
			return temp7;
		}
		else {
			opserr << "Wrong number of stress components to record!" << endln;
			return temp7;
		}
	}

	else {
		static Vector temp5(5);
		temp5 = this->getCommittedStress(isImplex);
		if (numOutput == 3)
		{
			static Vector temp3(3);
			temp3[0] = temp5[0];
			temp3[1] = temp5[1];
			temp3[2] = temp5[3];
			return temp3;
		}
		else if (numOutput == 4)
		{
			static Vector temp4(4);
			temp4[0] = temp5[0];
			temp4[1] = temp5[1];
			temp4[2] = temp5[2];
			temp4[3] = temp5[3];
			return temp4;
		}
		else if (numOutput == 5)
		{
			return temp5;
		}
		else {
			opserr << "Wrong number of stress components to record!" << endln;
			return temp5;
		}
	}
}
// end change by Alborz Ghofrani - UW 

const Vector& PIMYImplex::getImplexSateVariables(void)
{
	static Vector sv(3);
	sv(0) = lambda;
	sv(1) = chi;
	sv(2) = ksi;
	return sv;
}

const Vector& PIMYImplex::getCommittedStrain(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3)
		return currentStrain.t2Vector(1);
	else {
		static Vector workV(3), temp6(6);
		temp6 = currentStrain.t2Vector(1);
		workV[0] = temp6[0];
		workV[1] = temp6[1];
		workV[2] = temp6[3];
		return workV;
	}
}

// NOTE: surfaces[0] is not used
void PIMYImplex::setUpSurfaces(double* gredu)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
	double peakShearStrain = peakShearStrainx[matN];

	double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double pi = 3.14159265358979;
	double refStrain, peakShear, coneHeight;

	if (gredu == 0) {  //automatic generation of surfaces
		if (frictionAngle > 0) {
			double sinPhi = sin(frictionAngle * pi / 180.);
			double Mnys = 6. * sinPhi / (3. - sinPhi);
			residualPress = 3. * cohesion / (sqrt(2.) * Mnys);
			coneHeight = -(refPressure - residualPress);
			peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
			refStrain = (peakShearStrain * peakShear)
				/ (refShearModulus * peakShearStrain - peakShear);
		}

		else if (frictionAngle == 0.) { // cohesion = peakShearStrength
			peakShear = 2 * sqrt(2.) * cohesion / 3;
			refStrain = (peakShearStrain * peakShear)
				/ (refShearModulus * peakShearStrain - peakShear);
			residualPress = 0.;
		}

		double stressInc = peakShear / numOfSurfaces;

		for (int ii = 1; ii <= numOfSurfaces; ii++) {
			stress1 = ii * stressInc;
			stress2 = stress1 + stressInc;
			strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
			strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
			if (frictionAngle > 0.) size = 3. * stress1 / sqrt(2.) / coneHeight;
			else if (frictionAngle == 0.) size = 3. * stress1 / sqrt(2.);

			elasto_plast_modul = 2. * (stress2 - stress1) / (strain2 - strain1);

			if ((2. * refShearModulus - elasto_plast_modul) <= 0)
				plast_modul = UP_LIMIT;
			else
				plast_modul = (2. * refShearModulus * elasto_plast_modul) /
				(2. * refShearModulus - elasto_plast_modul);
			if (plast_modul < 0) plast_modul = 0;
			if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
			if (ii == numOfSurfaces) plast_modul = 0;

			static Vector temp(6);
			committedSurfaces[ii] = MultiYieldSurface(temp, size, plast_modul);
		}  // ii
	}
	else {  //user defined surfaces
		if (frictionAngle > 0) {   // ignore user defined frictionAngle
			int ii = 2 * (numOfSurfaces - 1);
			double tmax = refShearModulus * gredu[ii] * gredu[ii + 1];
			double Mnys = -(sqrt(3.) * tmax - 2. * cohesion) / refPressure;
			if (Mnys <= 0) {   // also ignore user defined cohesion
				cohesion = sqrt(3.) / 2 * tmax;
				frictionAngle = 0.;
				coneHeight = 1.;
				residualPress = 0.;
			}
			else {
				double sinPhi = 3 * Mnys / (6 + Mnys);
				if (sinPhi < 0. || sinPhi>1.) {
					opserr << "\nNDMaterial " << this->getTag() << ": Invalid friction angle, please modify ref. pressure or G/Gmax curve." << endln;
					exit(-1);
				}
				residualPress = 2. * cohesion / Mnys;
				if (residualPress < 0.01 * refPressure) residualPress = 0.01 * refPressure;
				coneHeight = -(refPressure - residualPress);
				frictionAngle = asin(sinPhi) * 180 / pi;
			}
		}
		else if (frictionAngle == 0.) {   // ignore user defined cohesion
			int ii = 2 * (numOfSurfaces - 1);
			double tmax = refShearModulus * gredu[ii] * gredu[ii + 1];
			cohesion = sqrt(3.) / 2 * tmax;
			coneHeight = 1.;
			residualPress = 0.;
		}

		/*
		opserr << "\nNDMaterial " <<this->getTag()<<": Friction angle = "<<frictionAngle
			<<", Cohesion = "<<cohesion<<"\n"<<endln;
		*/


		if (frictionAngle == 0.) pressDependCoeff = 0.; // ignore user defined pressDependCoeff

		for (int i = 1; i < numOfSurfaces; i++) {
			int ii = 2 * (i - 1);
			strain1 = gredu[ii];
			stress1 = refShearModulus * gredu[ii + 1] * strain1;
			strain2 = gredu[ii + 2];
			stress2 = refShearModulus * gredu[ii + 3] * strain2;

			size = sqrt(3.) * stress1 / coneHeight;
			elasto_plast_modul = 2. * (stress2 - stress1) / (strain2 - strain1);
			if ((2. * refShearModulus - elasto_plast_modul) <= 0)
				plast_modul = UP_LIMIT;
			else
				plast_modul = (2. * refShearModulus * elasto_plast_modul) /
				(2. * refShearModulus - elasto_plast_modul);
			if (plast_modul <= 0) {
				opserr << "\nNDMaterial " << this->getTag() << ": Surface " << i
					<< " has plastic modulus < 0.\n Please modify G/Gmax curve.\n" << endln;
				exit(-1);
			}
			if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;

			static Vector temp(6);
			committedSurfaces[i] = MultiYieldSurface(temp, size, plast_modul);

			if (i == (numOfSurfaces - 1)) {
				plast_modul = 0;
				size = sqrt(3.) * stress2 / coneHeight;
				committedSurfaces[i + 1] = MultiYieldSurface(temp, size, plast_modul);
			}
		}
	}

	residualPressx[matN] = residualPress;
	frictionAnglex[matN] = frictionAngle;
	cohesionx[matN] = cohesion;
}


double PIMYImplex::yieldFunc(const T2Vector& stress,
	const MultiYieldSurface* surfaces, int surfaceNum)
{
	static Vector temp(6);
	//temp = stress.deviator() - surfaces[surfaceNum].center();
	temp = stress.deviator();
	temp -= surfaces[surfaceNum].center();

	double sz = surfaces[surfaceNum].size();
	return 3. / 2. * (temp && temp) - sz * sz;
}


void PIMYImplex::deviatorScaling(T2Vector& stress, const MultiYieldSurface* surfaces,
	int surfaceNum, int count)
{
	count++;
	int numOfSurfaces = numOfSurfacesx[matN];

	double diff = yieldFunc(stress, surfaces, surfaceNum);

	if (surfaceNum < numOfSurfaces && diff < 0.) {
		double sz = surfaces[surfaceNum].size();
		double deviaSz = sqrt(sz * sz + diff);
		static Vector devia(6);
		devia = stress.deviator();
		static Vector temp(6);
		temp = devia - surfaces[surfaceNum].center();
		double coeff = (sz - deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
		devia.addVector(1.0, temp, coeff);
		stress.setData(devia, stress.volume());
		deviatorScaling(stress, surfaces, surfaceNum, count);  // recursive call
	}

	if (surfaceNum == numOfSurfaces && fabs(diff) > LOW_LIMIT) {
		double sz = surfaces[surfaceNum].size();
		static Vector newDevia(6);
		newDevia.addVector(0.0, stress.deviator(), sz / sqrt(diff + sz * sz));
		stress.setData(newDevia, stress.volume());
	}
}


void PIMYImplex::initSurfaceUpdate()
{
	if (committedActiveSurf == 0) return;

	int numOfSurfaces = numOfSurfacesx[matN];

	static Vector devia(6);
	devia = currentStress.deviator();
	double Ms = sqrt(3. / 2. * (devia && devia));
	static Vector newCenter(6);

	if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
		//newCenter = devia * (1. - committedSurfaces[activeSurfaceNum].size() / Ms);
		newCenter.addVector(0.0, devia, 1.0 - committedSurfaces[committedActiveSurf].size() / Ms);
		committedSurfaces[committedActiveSurf].setCenter(newCenter);
	}

	for (int i = 1; i < committedActiveSurf; i++) {
		newCenter = devia * (1. - committedSurfaces[i].size() / Ms);
		committedSurfaces[i].setCenter(newCenter);
	}
}


void PIMYImplex::paramScaling(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];

	if (frictionAngle == 0.) return;

	double conHeig = -(currentStress.volume() - residualPress);
	double scale = -conHeig / (refPressure - residualPress);

	scale = pow(scale, pressDependCoeff);
	refShearModulus *= scale;
	refBulkModulus *= scale;

	double plastModul, size;
	static Vector temp(6);
	for (int i = 1; i <= numOfSurfaces; i++) {
		plastModul = committedSurfaces[i].modulus() * scale;
		size = committedSurfaces[i].size() * conHeig;
		committedSurfaces[i] = MultiYieldSurface(temp, size, plastModul);
	}

}


void PIMYImplex::setTrialStress(T2Vector& stress)
{
	static Vector devia(6);
	//devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
	devia = stress.deviator();
	devia.addVector(1.0, subStrainRate.deviator(), 2. * refShearModulus);

	trialStress.setData(devia, stress.volume());
}


int PIMYImplex::setSubStrainRate(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];

	//if (activeSurfaceNum==numOfSurfaces) return 1;

	//if (strainRate==T2Vector()) return 0;
	if (strainRate.isZero()) return 0;

	double elast_plast_modulus;
	if (activeSurfaceNum == 0)
		elast_plast_modulus = 2 * refShearModulus;
	else {
		double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
		elast_plast_modulus = 2 * refShearModulus * plast_modulus
			/ (2 * refShearModulus + plast_modulus);
	}
	static Vector incre(6);
	//incre = strainRate.deviator()*elast_plast_modulus;
	incre.addVector(0.0, strainRate.deviator(), elast_plast_modulus);

	static T2Vector increStress;
	increStress.setData(incre, 0);
	double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
	double totalCross = 3. * increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross / singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	//incre = strainRate.t2Vector() / numOfSub;
	incre = strainRate.t2Vector();
	incre /= numOfSub;
	subStrainRate.setData(incre);

	return numOfSub;
}


void
PIMYImplex::getContactStress(T2Vector& contactStress)
{
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	static Vector devia(6);
	//devia = trialStress.deviator() - center;
	devia = trialStress.deviator();
	devia -= center;

	double Ms = sqrt(3. / 2. * (devia && devia));
	//devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;
	devia *= theSurfaces[activeSurfaceNum].size() / Ms;
	devia += center;

	// correct implex state variables
	ksi = theSurfaces[activeSurfaceNum].size() / Ms;

	contactStress.setData(devia, trialStress.volume());
}


int PIMYImplex::isLoadReversal(void)
{
	if (activeSurfaceNum == 0) return 0;

	static Vector surfaceNormal(6);
	getSurfaceNormal(currentStress, surfaceNormal);

	//(((trialStress.deviator() - currentStress.deviator()) && surfaceNormal) < 0)
	// return 1;
	static Vector a(6);
	a = trialStress.deviator();
	a -= currentStress.deviator();
	if ((a && surfaceNormal) < 0)
		return 1;

	return 0;
}


void
PIMYImplex::getSurfaceNormal(const T2Vector& stress, Vector& surfaceNormal)
{
	//Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();
	// return Q / sqrt(Q && Q);

	surfaceNormal = stress.deviator();
	surfaceNormal -= theSurfaces[activeSurfaceNum].center();
	surfaceNormal /= sqrt(surfaceNormal && surfaceNormal);
}


double PIMYImplex::getLoadingFunc(const T2Vector& contactStress,
	const Vector& surfaceNormal, int crossedSurface)
{
	double loadingFunc;
	double temp1 = 2. * refShearModulus;
	double temp2 = theSurfaces[activeSurfaceNum].modulus();

	//for crossing first surface
	double temp = temp1 + temp2;
	//loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp;
	static Vector tmp(6);
	tmp = trialStress.deviator();
	tmp -= contactStress.deviator();
	loadingFunc = (surfaceNormal && tmp) / temp;
	//for crossing more than one surface
	if (crossedSurface) {
		double temp3 = theSurfaces[activeSurfaceNum - 1].modulus();
		loadingFunc *= (temp3 - temp2) / temp3;
	}

	return loadingFunc;
}


void PIMYImplex::stressCorrection(int crossedSurface)
{
	static T2Vector contactStress;
	this->getContactStress(contactStress);
	static Vector surfaceNormal(6);
	this->getSurfaceNormal(contactStress, surfaceNormal);
	double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, crossedSurface);
	static Vector devia(6);

	// correct implex state variables
		// lambda
	lambda += loadingFunc;
		// chi
	devia.addVector(0.0, theSurfaces[activeSurfaceNum].center(), 1.0);
	devia.addVector(-1.0, contactStress.deviator(), 1.0);
	chi += sqrt(devia && devia);

	//devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
	devia.addVector(0.0, surfaceNormal, -2.0 * refShearModulus * loadingFunc);
	devia += trialStress.deviator();

	trialStress.setData(devia, trialStress.volume());
	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		stressCorrection(1);  //recursive call
	}
}


void PIMYImplex::updateActiveSurface(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];

	if (activeSurfaceNum == numOfSurfaces) return;

	double A, B, C, X;
	static T2Vector direction;
	static Vector t1(6);
	static Vector t2(6);
	static Vector temp(6);
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector outcenter(6);
	outcenter = theSurfaces[activeSurfaceNum + 1].center();
	double outsize = theSurfaces[activeSurfaceNum + 1].size();

	//t1 = trialStress.deviator() - center;
	//t2 = center - outcenter;
	t1 = trialStress.deviator();
	t1 -= center;
	t2 = center;
	t2 -= outcenter;

	A = t1 && t1;
	B = 2. * (t1 && t2);
	C = (t2 && t2) - 2. / 3. * outsize * outsize;
	X = secondOrderEqn(A, B, C, 0);
	if (fabs(X - 1.) < LOW_LIMIT) X = 1.;
	if (X < 1.) {
		opserr << "FATAL:PIMYImplex::updateActiveSurface(): error in Direction of surface motion."
			<< endln;
		exit(-1);
	}

	//temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	temp = center;
	temp.addVector(1.0, t1, X);
	temp *= (1.0 - size / outsize);
	t2 = center;
	t2.addVector(1.0, outcenter, -size / outsize);
	temp -= t2;

	direction.setData(temp);

	if (direction.deviatorLength() < LOW_LIMIT) return;

	temp = direction.deviator();
	A = temp && temp;
	B = -2 * (t1 && temp);
	if (fabs(B) < LOW_LIMIT) B = 0.;
	C = (t1 && t1) - 2. / 3. * size * size;
	if (fabs(C) < LOW_LIMIT || fabs(C) / (t1 && t1) < LOW_LIMIT) return;

	// Added by Alborz Ghofrani UW to avoid solving quadratic equations with very small coefficients B and C
	if ((fabs(B) < 1.0e-10) && (fabs(C) < 1.0e-10)) {
		return;
	}
	// End Alborz Ghofrani UW

	if (B > 0. || C < 0.) {
		opserr << "FATAL:PIMYImplex::updateActiveSurface(): error in surface motion.\n"
			<< "A= " << A << " B= " << B << " C= " << C << " (t1&&t1)= " << (t1 && t1) << endln;
		exit(-1);
	}
	X = secondOrderEqn(A, B, C, 1);

	//center += temp * X;
	center.addVector(1.0, temp, X);
	theSurfaces[activeSurfaceNum].setCenter(center);
}


void PIMYImplex::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

	static Vector devia(6);
	devia = currentStress.deviator();
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector newcenter(6);

	for (int i = 1; i < activeSurfaceNum; i++) {
		//newcenter = devia - (devia - center) * theSurfaces[i].size() / size;
		newcenter = center;
		newcenter -= devia;
		newcenter *= theSurfaces[i].size() / size;
		newcenter += devia;

		theSurfaces[i].setCenter(newcenter);
	}
}


int PIMYImplex::isCrossingNextSurface(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	if (activeSurfaceNum == numOfSurfaces) return 0;

	if (yieldFunc(trialStress, theSurfaces, activeSurfaceNum + 1) > 0) return 1;

	return 0;
}


int PIMYImplex::implicitSress(void) {

	int numOfSurfaces = numOfSurfacesx[matN];

	// initialize current step implex state variables
	lambda = 0.0;
	chi = 0.0;
	ksi = 0.0;

	for (int i = 1; i <= numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
	activeSurfaceNum = committedActiveSurf;
	subStrainRate = strainRate;
	setTrialStress(currentStress);
	if (isLoadReversal()) {
		updateInnerSurface();
		activeSurfaceNum = 0;
	}
	int numSubIncre = setSubStrainRate();

	for (int i = 0; i < numSubIncre; i++) {
		if (i == 0)
			setTrialStress(currentStress);
		else
			setTrialStress(trialStress);
		if (activeSurfaceNum == 0 && !isCrossingNextSurface()) continue;
		if (activeSurfaceNum == 0) activeSurfaceNum++;
		stressCorrection(0);
		updateActiveSurface();
	}
	//volume stress change
	double volum = refBulkModulus * (strainRate.volume() * 3.);
	volum += currentStress.volume();
	//if (volum > 0) volum = 0.;
	trialStress.setData(trialStress.deviator(), volum);
	return 0;
}