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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/VonMisesDMM.cpp$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	June 2022
//
// Description: This file contains the implementation for the VonMisesDMM sub-class.

#include "VonMisesDMM.h"

// Public methods
	// full constructor
VonMisesDMM::VonMisesDMM(int tag, double r0,
	double K0, double G0, double P0, double m0,
	int T0, DataDrivenNestedSurfaces* ys,
	int ddtype, int itype)
	:MultiYieldSurfaceHardeningSoftening(tag, ND_TAG_VonMisesDMM, ys, ddtype, itype)
{
	rho = r0;  TNYS = T0; Kref = K0; Gref = G0;
	Pref = P0; Modn = m0; 
	if (theSurfaces->getCohesion() < 0) {
		opserr << "FATAL:VonMisesDMM: cohesion <= 0\n";
		opserr << "Yield strength (cohesion) cannot be zero...\n";
		exit(-1);
	}
}

// null constructor
VonMisesDMM::VonMisesDMM(void)
	:MultiYieldSurfaceHardeningSoftening()
{
}

// destructor
VonMisesDMM::~VonMisesDMM(void)
{
}

// return object info
NDMaterial* VonMisesDMM::getCopy(void) {

	VonMisesDMM* copy = nullptr;
	// set material type and order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}
	copy = new VonMisesDMM(*this);
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theSurfaces->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: VonMisesDMM: getCopy: returned NULLPTR! couldn't copy the material...\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

NDMaterial* VonMisesDMM::getCopy(const char* type) {

	VonMisesDMM* copy = nullptr;
	if (strcmp(type, "VonMisesDMM") == 0 || strcmp(type, "ThreeDimensional") == 0) {
		nOrd = 6; // ThreeDimensional
		copy = new VonMisesDMM(*this);
	}
	else if (strcmp(type, "PlaneStrain") == 0) {
		nOrd = 3; // PlaneStrain
		copy = new VonMisesDMM(*this);
	}
	else {
		opserr << "FATAL: VonMisesDMM getCopy: unsupported material type = " << type << "\n";
		opserr << "VonMisesDMM materialType can be --> ThreeDimensional (1), PlaneStrain (0)\n";
		exit(-1);
	}
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theSurfaces->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: VonMisesDMM: getCopy: returned NULLPTR! couldn't copy the material.\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

// parallel message passing
int VonMisesDMM::sendSelf(int commitTag, Channel& theChannel) {

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

int VonMisesDMM::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {

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
void VonMisesDMM::Print(OPS_Stream& s, int flag) {
	s << "von Mises Multi-Yield Surface Hardening-Softening nD Material\n";
}

int VonMisesDMM::setParameter(const char** argv, int argc, Parameter& param) {

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

int VonMisesDMM::updateParameter(int responseID, Information& info) {

	if (responseID == 1) {
		// update material stage and if nonlinear set up yield surfaces 
		materialStage = info.theInt;
		if (materialStage > 0) {
			// set-up yield surfaces
			theSurfaces = theData->setUpYieldSurfaces(this);
		}
	}
	else if (responseID == 12) {
		//frictionAnglex[matN] = info.theDouble;
	}
	else if (responseID == 13) {
		theSurfaces->updateCohesion(info.theDouble);
	}

	// also go through the base-class method
	return MultiYieldSurfaceHardeningSoftening::updateParameter(responseID, info);
}

// Private methods
	// yield surface operations
double VonMisesDMM::yieldFunction(const Vector& stress, const int num_yield_surface, bool yield_Stress = false) {
	// get material constants
	double Pc = theSurfaces->getCohesion();
	double M = theSurfaces->getHP(num_yield_surface);

	// Evaluate and return the yield surface value
	double yf = sqrt(2./3.) * M * Pc;
	if (!yield_Stress) {
		Vector zeta = getStressDeviator(stress, num_yield_surface);
		yf = sqrt(TensorM::dotdot(zeta, zeta)) - yf;
	}
	return yf;
}

Vector VonMisesDMM::get_dF_dS(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t stress
	Vector zeta = getStressDeviator(stress, num_yield_surface);
	Vector alpha = tools::getColumnVector(num_yield_surface, sv.alpha);
	Vector dfds = zeta / sqrt(TensorM::dotdot(zeta, zeta));
	dfds = dfds - (1. / 3.) * (TensorM::dotdot(dfds, alpha)) * TensorM::I(6);
	return dfds;
}

Vector VonMisesDMM::get_dF_dA(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t alpha(backstress)
	Vector zeta = getStressDeviator(stress, num_yield_surface);
	Vector alpha = tools::getColumnVector(num_yield_surface, sv.alpha);
	Vector dfda = -1 * getMeanStress(stress) * (zeta / sqrt(TensorM::dotdot(zeta, zeta)));
	return dfda;
}

Vector VonMisesDMM::get_dH_dA(const Vector& stress, const int num_yield_surface) {
	// Return the rate of alpha(backstress)
	double radius = getSizeYS(num_yield_surface);

	if (radius == 0) {
		opserr << "\n";
		opserr << "FATAL:VonMisesDMM::get_dR_dA\n";
		opserr << "curr_sz = 0  \n";
		opserr << "N_active_ys: " << num_yield_surface << "\n";
		opserr << "Total theSurfaces  : " << theSurfaces->getTNYS() << "\n";
		opserr << "\n";
		exit(-1);
	}

	double Pavg = getMeanStress(stress);
	double denominator = 0.0;
	Vector dhda(6);
	Vector zeta = getStressDeviator(stress, num_yield_surface);
	Vector Q_prime = zeta / sqrt(TensorM::dotdot(zeta, zeta));

	if (num_yield_surface >= theSurfaces->getTNYS()) {
		double H_prime = theSurfaces->getHref(theSurfaces->getTNYS());
		dhda = (H_prime / Pavg) * Q_prime;
	}
	else {
		double next_radius = getSizeYS(num_yield_surface + 1);
		double H_prime = theSurfaces->getHref(num_yield_surface);
		Vector next_zeta = getStressDeviator(stress, num_yield_surface + 1);

		dhda = ((next_radius / radius) * (zeta)) - (next_zeta);
		denominator = TensorM::dotdot(Q_prime, dhda);

		if (denominator == 0) {
			opserr << "\n";
			opserr << "FATAL:VonMisesDMM::get_dR_dA\n";
			opserr << "denominator = 0  \n";
			opserr << "N_active_ys: " << num_yield_surface << "\n";
			opserr << "Total theSurfaces  : " << theSurfaces->getTNYS() << "\n";
			opserr << "\n";
			exit(-1);
		}

		dhda = (H_prime / (Pavg * denominator)) * dhda;
	}
	return dhda;
}

Vector VonMisesDMM::get_dP_dS(const Vector& stress, const int num_yield_surface) {
	// Return the plastic flow direction 
	return get_dF_dS(stress, num_yield_surface); // Associated flow
}