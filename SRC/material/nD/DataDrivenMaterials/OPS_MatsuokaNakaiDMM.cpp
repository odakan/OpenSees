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

// Written: Onur Deniz Akan         (onur.akan@iusspavia.it) 
//          Guido Camata      
//          Enrico Spacone
//          Carlo G. Lai
//
// Created: June 2022
//
// Description: This file contains the implementation for the OPS_MatsuokaNakaiDMM function.

#include "MatsuokaNakaiDMM.h"

#define OPS_Export 

OPS_Export void*
OPS_MatsuokaNakaiDMM(void)
{
	// display kudos
	static int kudos = 0;
	if (++kudos == 1)
		opserr << "Matsuoka-Nakai Multi-Yield Surface Hardening-Softening nDmaterial - Written: OD.Akan, G.Camata, E.Spacone, CG.Lai, IUSS Pavia\n";

	// initialize pointers
	NDMaterial* theMaterial = nullptr;					// Pointer to an nD material to be returned
	DataDrivenNestedSurfaces* theSurfaces = nullptr;	// Pointer to a nested yield surface object

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 1)
	{
		opserr << "Want: nDMaterial MatsuokaNakaiDMM tag? rho? Kref? Gref? Pref? modn? cohesion? TNYS? <HParams?> <-intType type?>\n";
		return theMaterial;
	}

	int tag;
	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0)
	{
		opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid tag? (must be an integer)" << endln;
		return theMaterial;
	}

	if (numArgs < 9)
	{
		opserr << "Want: nDMaterial MatsuokaNakaiDMM " << tag << " rho? Kref? Gref? Pref? modn? cohesion? TNYS? <HParams?> <-intType type?>\n";
		return theMaterial;
	}

	// recieve mass density, rho
	double rho = 0.0;
	numData = 1;
	if (OPS_GetDoubleInput(&numData, &rho) < 0)
	{
		opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid rho? (must be a double)\n";
		return theMaterial;
	}

	// recieve Kref, Gref, Pref, nmod, cohesion
	double params[5]; params[0] = 0.; params[1] = 0.; params[2] = 0.; params[3] = 0.; params[4] = 0.;
	numData = 5;
	if (OPS_GetDoubleInput(&numData, &params[0]) < 0)
	{
		opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid Kref? Gref? Pref? modn? cohesion? (must be doubles)\n";
		return theMaterial;
	}

	// recieve TNYS
	int TNYS = 0;
	int maxTNYS = 1000;
	numData = 1;
	if (OPS_GetIntInput(&numData, &TNYS) < 0)
	{
		opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid TNYS? (must be an integer -" << maxTNYS << " <= TNYS <= " << maxTNYS << " && TNYS != 0)\n";
		return theMaterial;
	}
	if (abs(TNYS) > maxTNYS || TNYS == 0)
	{
		opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid TNYS? (must be an integer -" << maxTNYS << " <= TNYS <= " << maxTNYS << " && TNYS != 0)\n";
		return theMaterial;
	}

	// continue with the optional inputs
		// user defined yield surfaces
	double* HParams = nullptr;
	if (TNYS < 0)
	{
		// recieve hardening parameters
		numData = int(TNYS);
		HParams = new double[++TNYS];
		if (OPS_GetDoubleInput(&numData, HParams) < 0)
		{
			opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid HParams? (must be " << TNYS << " many doubles)\n";
			return theMaterial;
		}
	}
	//other inputs
	int integrationType = 0;                                                        // implicit by default
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* inputstring = OPS_GetString();
		// recieve material integration type
		if (strcmp(inputstring, "-intType") == 0) {									// #1 read type of integration 
			numData = 1;
			if (OPS_GetIntInput(&numData, &integrationType) < 0) {
				opserr << "OPS_MatsuokaNakaiDMM: WARNING: invalid integration type after -intType (must be an integer)\n";
				return theMaterial;
			}
		}
	}
	// all inputs recieved!

	// create a nested yield surfaces object
	theSurfaces = new DataDrivenNestedSurfaces();

	// free the memory
	delete[] HParams;

	// create a MatsuokaNakaiDMM nDmaterial object
	theMaterial = new MatsuokaNakaiDMM();

	if (theMaterial == nullptr) {
		opserr << "FATAL!-- > OPS_MatsuokaNakaiDMM: cannot create MatsuokaNakaiDMM material with tag: " << tag << "\n";
	}

	return theMaterial;
}