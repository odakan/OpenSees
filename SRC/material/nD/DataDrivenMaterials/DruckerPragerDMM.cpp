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

#include "DruckerPragerDMM.h"

//	-- TO DO LIST --
// 3) CODE sendSelf and recvSelf
//	- pack and unpack all the variables


// Public methods
	// full constructor
DruckerPragerDMM::DruckerPragerDMM(int tag, double r0,
	double K0, double G0, double P0, double m0,
	DataDrivenNestedSurfaces* ys,
	int ddtype, int itype, bool verbosity)
	:MultiYieldSurfaceHardeningSoftening(tag, ND_TAG_DruckerPragerDMM, r0,
		K0, G0, P0, m0, ys, ddtype, itype, verbosity)
{

}

	// null constructor
DruckerPragerDMM::DruckerPragerDMM()
	:MultiYieldSurfaceHardeningSoftening()
{
}

	// destructor
DruckerPragerDMM :: ~DruckerPragerDMM()
{
}

	// return object info
NDMaterial* DruckerPragerDMM::getCopy(void) {
	DruckerPragerDMM* copy = nullptr;
	// set material type and order
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}
	copy = new DruckerPragerDMM(*this);
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: DruckerPragerDMM: getCopy: returned NULLPTR! couldn't copy the material...\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

NDMaterial* DruckerPragerDMM::getCopy(const char* type) {
	DruckerPragerDMM* copy = nullptr;
	if (strcmp(type, "DruckerPragerDMM") == 0 || strcmp(type, "ThreeDimensional") == 0) {
		nOrd = 6; // ThreeDimensional
		copy = new DruckerPragerDMM(*this);
	}
	else if (strcmp(type, "PlaneStrain") == 0) {
		nOrd = 3; // PlaneStrain
		copy = new DruckerPragerDMM(*this);
	}
	else {
		opserr << "FATAL: DruckerPragerDMM getCopy: unsupported material type = " << type << "\n";
		opserr << "DruckerPragerDMM materialType can be --> ThreeDimensional (1), PlaneStrain (0)\n";
		exit(-1);
	}
	// done

	if (copy != nullptr) {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
		return copy;
	}
	else {
		opserr << "WARNING: DruckerPragerDMM: getCopy: returned NULLPTR! couldn't copy the material.\n";
		opserr << "Moving on...\n";
		opserr << "..\n";
		opserr << ".\n\n";
		return 0;
	}
}

	// parallel message passing
int DruckerPragerDMM::sendSelf(int commitTag, Channel& theChannel) {
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

int DruckerPragerDMM::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
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
void DruckerPragerDMM::Print(OPS_Stream& s, int flag) {
	s << "Drucker-Prager Multi-Yield Surface Hardening-Softening nD Material\n";
}

int DruckerPragerDMM::setParameter(const char** argv, int argc, Parameter& param) {
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

int DruckerPragerDMM::updateParameter(int responseID, Information& info) {
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
Vector DruckerPragerDMM::getShiftedDeviator(const Vector& stress, const int num_ys) {

	// Prevost DYNA1D 5.2.2 eqn. 3
	double P_bar = getMeanStress(stress) - ys.getCohesion();
	return  getStressDeviator(stress) - P_bar * ys.getAlpha(num_ys);
}

	// yield surface operations
double DruckerPragerDMM::yieldFunction(const Vector& stress, const int num_yield_surface,	bool yield_stress) {
	// Evaluate and return the yield surface value
	double yield_function = 0;
	return  yield_function;
}

Vector DruckerPragerDMM::get_dF_dS(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t stress
	Vector dfds(6);
	return dfds;
}

Vector DruckerPragerDMM::get_dF_dA(const Vector& stress, const int num_yield_surface) {
	// Return the normal to the yield surface w.r.t alpha(backstress)
	Vector dfda(6);
	return dfda;
}

Vector DruckerPragerDMM::get_dH_dA(const Vector& stress, const int num_yield_surface) {
	// Return the rate of alpha(backstress)
	Vector direction(6);
	return direction;
}

Vector DruckerPragerDMM::get_dP_dS(const Vector& stress, const int num_yield_surface) {
	// Return the plastic flow direction
	Vector mm(6);
	return mm;
}