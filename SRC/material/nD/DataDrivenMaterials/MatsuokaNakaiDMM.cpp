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

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	June 2022
//
// Description: This file contains the implementation for the MatsuokaNakaiDMM class.

#include "MatsuokaNakaiDMM.h"

// Public methods
	// full constructor
MatsuokaNakaiDMM::MatsuokaNakaiDMM(int tag, double r0,
	double K0, double G0, double P0, double m0,
	int T0, DataDrivenNestedSurfaces* ys,
	int ddtype, int itype, bool verbosity)
	:MultiYieldSurfaceHardeningSoftening(tag, ND_TAG_MatsuokaNakaiDMM, r0,
		K0, G0, P0, m0, T0, ys, ddtype, itype, verbosity)
{

}

// null constructor
MatsuokaNakaiDMM::MatsuokaNakaiDMM()
	:MultiYieldSurfaceHardeningSoftening()
{
}

// destructor
MatsuokaNakaiDMM::~MatsuokaNakaiDMM()
{
}

// return object info
NDMaterial* MatsuokaNakaiDMM::getCopy(void) {
	MatsuokaNakaiDMM* copy = nullptr;
	// set material type and order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}
	copy = new MatsuokaNakaiDMM(*this);
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: MatsuokaNakaiDMM: getCopy: returned NULLPTR! couldn't copy the material...\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

NDMaterial* MatsuokaNakaiDMM::getCopy(const char* type) {
	MatsuokaNakaiDMM* copy = nullptr;
	if (strcmp(type, "MatsuokaNakaiDMM") == 0 || strcmp(type, "ThreeDimensional") == 0) {
		nOrd = 6; // ThreeDimensional
		copy = new MatsuokaNakaiDMM(*this);
	}
	else if (strcmp(type, "PlaneStrain") == 0) {
		nOrd = 3; // PlaneStrain
		copy = new MatsuokaNakaiDMM(*this);
	}
	else {
		opserr << "FATAL: MatsuokaNakaiDMM getCopy: unsupported material type = " << type << "\n";
		opserr << "MatsuokaNakaiDMM materialType can be --> ThreeDimensional (1), PlaneStrain (0)\n";
		exit(-1);
	}
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: MatsuokaNakaiDMM: getCopy: returned NULLPTR! couldn't copy the material.\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

// parallel message passing
int MatsuokaNakaiDMM::sendSelf(int commitTag, Channel& theChannel) {
	int res = 0;

	// place data in a vector
	static Vector data(45);
	data(0) = this->getTag();


	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int MatsuokaNakaiDMM::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	int res = 0;

	// receive data
	static Vector data(45);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}

	// set member variables
	this->setTag((int)data(0));


	return 0;
}

// miscellaneous
void MatsuokaNakaiDMM::Print(OPS_Stream& s, int flag) {
	s << "von Mises Multi-Yield Surface Hardening-Softening nD Material\n";
}

int MatsuokaNakaiDMM::setParameter(const char** argv, int argc, Parameter& param) {
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	// check for material tag
	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "frictionAngle") == 0) {
			return param.addObject(12, this);
		}
		else if (strcmp(argv[0], "cohesion") == 0) {
			return param.addObject(13, this);
		}
	}

	// also go through the base-class method
	return MultiYieldSurfaceHardeningSoftening::setParameter(argv, argc, param);
}

int MatsuokaNakaiDMM::updateParameter(int responseID, Information& info) {
	if (responseID == 12) {
		//frictionAnglex[matN] = info.theDouble;
	}
	else if (responseID == 13) {
		ys.setCohesion(info.theDouble);
	}

	// also go through the base-class method
	return MultiYieldSurfaceHardeningSoftening::updateParameter(responseID, info);
}

// Private methods
	// yield surface operations
double MatsuokaNakaiDMM::yieldFunction(const Vector& stress, const int num_yield_surface, bool yield_stress) {
	// Evaluate and return the yield surface value
	double yield_function = 0;
	return  yield_function;
}

Vector MatsuokaNakaiDMM::get_dF_dS(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t stress
	Vector dfds(6);
	return dfds;
}

Vector MatsuokaNakaiDMM::get_dF_dA(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t alpha(backstress)
	Vector dfda(6);
	return dfda;
}

Vector MatsuokaNakaiDMM::get_dH_dA(const Vector& stress, const int num_yield_surface) {
	// Return the rate of alpha(backstress)
	Vector direction(6);
	return direction;
}

Vector MatsuokaNakaiDMM::get_dP_dS(const Vector& stress, const int num_yield_surface) {
	// Return the plastic flow direction
	Vector mm(6);
	return mm;
}