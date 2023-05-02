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
	DataDrivenNestedSurfaces* ys,
	int ddtype, int itype, bool verbosity)
	:MultiYieldSurfaceHardeningSoftening(tag, ND_TAG_VonMisesDMM, r0,
		K0, G0, P0, m0, ys, ddtype, itype, verbosity)
{
	

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
		theData->checkin(); // do check-in
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
		theData->checkin(); // do check-in
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
	// the get methods
Vector VonMisesDMM::getShiftedDeviator(const Vector& stress, const int num_ys) { return getStressDeviator(stress) - ys.getAlpha(num_ys); }

	// yield surface operations
double VonMisesDMM::yieldFunction(const Vector& stress, const int num_ys, bool yield_stress = false) {
	// if yield_stress is false, return the yield function value. Otherwise, return the yield strength.

	// get the current limit stress
	double strength = sqrt(2.0 / 3.0) * ys.getTau(num_ys);

	// evaluate and return the yield surface value
	if (!yield_stress) {
		Vector zeta = getShiftedDeviator(stress, num_ys);
		strength = sqrt(TensorM::dotdot(zeta, zeta)) - strength;
	}

	// done
	return strength;
}

Vector VonMisesDMM::get_dF_dS(const Vector& stress, const int num_ys) {
	// Return the normal to the yield surface w.r.t stress
	
	// initialize the tensors
	Vector dfds = Vector(6);							// normal tensor
	Vector zeta = getShiftedDeviator(stress, num_ys);	// shifted stress deviator tensor

	// compute the normal tensor
	dfds = zeta / sqrt(TensorM::dotdot(zeta, zeta));

	// done
	return dfds;
}

Vector VonMisesDMM::get_dF_dA(const Vector& stress, const int num_ys) {
	// Return the normal to the yield surface w.r.t alpha(backstress)
	
	// initialize the tensors
	Vector dfda = Vector(6);							// normal tensor
	Vector zeta = getShiftedDeviator(stress, num_ys);	// shifted stress deviator tensor

	// compute the normal tensor
	dfda = -1 * zeta / sqrt(TensorM::dotdot(zeta, zeta));

	// done
	return dfda;
}

Vector VonMisesDMM::get_dH_dA(const Vector& stress, const int num_ys) {
	// Return the rate of alpha(backstress)

	// initialize the tensors
	Vector dhda(6);												// normal tensor

	// compute the normal
	if (num_ys >= ys.getTNYS()) {
		// get material constants
		double H_prime = ys.getEta(ys.getTNYS());

		// initialize tensors
		Vector current_zeta = getShiftedDeviator(stress, num_ys);	// shifted stress deviator tensor (current)

		// compute normal
		Vector Q_prime = current_zeta / sqrt(TensorM::dotdot(current_zeta, current_zeta));
		dhda = H_prime * Q_prime;
	}
	else {
		// get material constants
		double current_strength = getSizeYS(num_ys);
		double next_strength = getSizeYS(num_ys + 1);
		double H_prime = ys.getEta(num_ys);

		// do check
		if (current_strength < ABSOLUTE_TOLERANCE) {
			opserr << "FATAL: VonMisesDMM::get_dH_dA() - current yield surface (no. " << num_ys << ") returned a yield strength of " << current_strength << "!\n";
			exit(-1);
		}
		else if (next_strength < ABSOLUTE_TOLERANCE) {
			opserr << "FATAL: VonMisesDMM::get_dH_dA() - current yield surface (no. " << num_ys + 1 << ") returned a yield strength of " << next_strength << "!\n";
			exit(-1);
		}

		// initialize tensors
		Vector current_zeta = getShiftedDeviator(stress, num_ys);	// shifted stress deviator tensor (current)
		Vector next_zeta = getShiftedDeviator(stress, num_ys + 1);	// shifted stress deviator tensor (next)

		// compute normal
		Vector Q_prime = current_zeta / sqrt(TensorM::dotdot(current_zeta, current_zeta));
		//Vector direction = (next_strength / current_strength) * current_zeta - next_zeta; // Prevost
		Vector direction = current_zeta - (current_strength / next_strength) * next_zeta; // Gu et al.
		double denominator = TensorM::dotdot(Q_prime, direction);
		if (denominator == 0.0) {
			opserr << "\nFATAL: VonMisesDMM::get_dH_dA() - division by 0 while computing the rate of backstress!\n";
			opserr << "Denominator [Q':mu] = "<< denominator << "\n";
			opserr << "Q' = " << Q_prime;
			opserr << "mu = " << direction;
			exit(-1);
		}
		dhda = (H_prime / denominator) * direction;
	}

	// done
	return dhda;
}

Vector VonMisesDMM::get_dP_dS(const Vector& stress, const int num_ys) {
	// Return the plastic flow direction (normal to the plastic potential w.r.t stress)

	// initialize the tensors
	Vector dpds = Vector(6);	// normal tensor

	// compute the normal to the plastic potential
	if (ys.isNonAssociated()) {
		// if dilatancy > 1 -> theta_bar < 0 [contraction] else if if dilatancy < 1 -> theta_bar > 0 [dilation]
		Vector zeta = getShiftedDeviator(stress, num_ys);	// shifted stress deviator tensor
		Vector tau = getStressDeviator(stress);	// shifted stress deviator tensor
		Vector alpha = ys.getAlpha(num_ys);

		// compute P"
		double dilatancy = ys.getTheta(num_ys);				// dilatancy multiplier
		double theta = sqrt(TensorM::dotdot(tau, tau));
		double theta_bar = dilatancy;
		double P_double_prime = ((pow((theta / theta_bar), 2) - 1) / (pow((theta / theta_bar), 2) + 1)) / 3.0;
		//opserr << "theta = " << theta << "\n";
		//opserr << "P_double_prime = " << P_double_prime * 3 << "\n";

		// evaluate the derivative dPdS = Q' + P" * kronecker
		dpds = zeta / sqrt(TensorM::dotdot(zeta, zeta));
		dpds += 1.0 / 3.0 * ( dilatancy) * TensorM::I(6);
	}
	else {
		dpds = get_dF_dS(stress, num_ys); // Associated flow
	}

	return dpds;
}