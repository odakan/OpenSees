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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MultiYieldSurfaceHardeningSoftening.cpp$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	June 2022
//
// Description: This file contains the implementation for the MultiYieldSurfaceHardeningSoftening class.

#include "MultiYieldSurfaceHardeningSoftening.h"
#include <MaterialResponse.h>


// Public methods
	// full constructor
MultiYieldSurfaceHardeningSoftening::MultiYieldSurfaceHardeningSoftening(int tag, int classTag,
	double r0, double K0, double G0, double P0, double m0,
	DataDrivenNestedSurfaces* data, int ddtype, int itype, bool verbosity)
	:NDMaterial(tag, classTag), rho(r0), Kref(K0), Gref(G0), 
	Pref(P0), Modn(m0), theData(data), beVerbose(verbosity)
{
	// handle solution options
		// handle dimension
	if (OPS_GetNDM() == 2) {		// PlaneStrain
		nOrd = 3;
	}
	else if (OPS_GetNDM() == 3) {	// ThreeDimensional
		nOrd = 6;
	}
	else {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - unknown model dimension...\n";
		exit(-1);
	}

		// material integration
	if (itype == 1) {
		use_implex = true;
	}
	else {
		use_implex = false;
	}

		// material data
	if (ddtype == 1) {						// passive
		use_active_approach = false;
		use_data_driven_surface = true;
	}
	else if (ddtype > 0 && ddtype < 1) {	// active
		use_active_approach = true;
		use_data_driven_surface = true;
	}
	else if (ddtype == 0) {					// automatic surface
		use_active_approach = false;
		use_data_driven_surface = false;
		use_user_custom_surface = false;
	}
	else if (ddtype == -1) {				// user custom surface
		use_active_approach = false;
		use_data_driven_surface = false;
		use_user_custom_surface = true;
	}
	else {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - unknown material data driver type! ddtype is = " << ddtype << "\n";
		exit(-1);
	}

	// do some checks
	if (theData == nullptr) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - yield surfaces are missing!\n";
		opserr << "The pointer that points at the yield surface: ys = " << theData << "\n";
		exit(-1);
	}
	else {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
	}
	if (nOrd != 3 && nOrd != 6) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - dimension error\n";
		opserr << "Material order has to be 3 or 6, but it is nOrd = " << nOrd << "\n";
		exit(-1);
	}
	if (Kref <= 0) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - Kref <= 0\n";
		exit(-1);
	}
	if (Gref <= 0) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - Gref <= 0\n";
		exit(-1);
	}
	if (Pref <= 0) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening() - Pref <= 0\n";
		opserr << "Use Pref = 101.325 kPa?\n";
		exit(-1);
	}
	if (Modn < 0) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening() - modn < 0\n";
		opserr << "modn is assumed to be zero...\n";
		Modn = 0.0;
	}
	if (rho < 0) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening() - mass density < 0\n";
		opserr << "mass density (rho) is assumed to be zero...\n";
		rho = 0.;
	}

	// reset internal parameters
	revertToStart();
}

	// null constructor
MultiYieldSurfaceHardeningSoftening::MultiYieldSurfaceHardeningSoftening(void)
	: NDMaterial()
{

}

	// destructor
MultiYieldSurfaceHardeningSoftening::~MultiYieldSurfaceHardeningSoftening(void)
{
	// free the memory allocated to the data-driven nested yield surface object
	if (theData != nullptr) {
		if (theData->canDelete()) {	// make sure that no other material is using the object.
			delete[] theData;		// if permitted, delete the object
			theData = nullptr;		// and reset pointer to nullptr.
		}
		else {						// if not, inform the object about								
			theData->checkout();	// leave and check this instance out
			theData = nullptr;		// and reset pointer to nullptr.
		}
	}
}

	// iteration control
int MultiYieldSurfaceHardeningSoftening::commitState(void) {

	// do the implicit correction if impl-ex
	if (use_implex) {
		// update material internal variables
		if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::commitState() -> IMPL-EX: implicit correction...\n"; }
		updateInternal(false, false);  // explicit_phase?, do_tangent?
	}

	// commit internal variables
	if (sv.commitState()) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::commitState() - material internal variables failed commit!\n";
		exit(-1);
	}
	
	if (ys.commitState()) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::commitState() - yield surface object failed commit!\n";
		exit(-1);
	}

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::commitState() -> printStats() ->\n"; sv.printStats(false); ys.printStats(false); }

	// done
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::revertToLastCommit(void) {

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::revertToLastCommit()\n"; }

	// restore committed internal variables
	if (sv.revertToLastCommit()) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::revertToLastCommit() - material internal variables failed reverting to last commit!\n";
		exit(-1);
	}

	if (ys.revertToLastCommit()) {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::revertToLastCommit() - yield surface object failed reverting to last commit!\n";
		exit(-1);
	}

	// done
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::revertToStart(void) {
	// reset state variables
	sv = MaterialStateVariables(nOrd);
	ys = YieldSurfacePackage();
	materialStage = 0;

	// reset elastoplastic stiffness
	sv.Ce = Kref * TensorM::IIvol(nOrd) + 2 * Gref * TensorM::IIdev(nOrd);	// compute initial elastic modulus
	sv.Cep = sv.Ce;

	// done
	return 0;
}

	// set material strain
int MultiYieldSurfaceHardeningSoftening::setTrialStrain(const Vector& strain) {
	
	// set the strain state
	sv.eps = strain;

	// implex time step
	if (!sv.dtime_is_user_defined) {
		sv.dtime_n = ops_Dt;
		if (!sv.dtime_first_set) {
			sv.dtime_n_commit = sv.dtime_n;
			sv.dtime_first_set = true;
		}
	}
	
	// trigger material update
	if (use_implex) {
		if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::setTrialStrain() -> IMPL-EX: explicit stage...\n"; }
		updateInternal(true, true);
		sv.sig_implex = sv.sig; // save stress for output
	}
	else {
		if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::setTrialStrain() -> IMPLICIT...\n"; }
		if (use_numerical_tangent) {
			// for the implicit case we use the numerical tangent... seems more stable
			static Vector strain(3);
			static Matrix Cnum(3, 3);
			strain = sv.eps;
			for (int j = 0; j < 3; ++j) {
				sv.eps(j) = strain(j) + SMALL_PERTURBATION;
				updateInternal(false, false);
				for (int i = 0; i < 3; ++i)
					Cnum(i, j) = sv.sig(i);
				sv.eps(j) = strain(j) - SMALL_PERTURBATION;
				updateInternal(false, false);
				for (int i = 0; i < 3; ++i)
					Cnum(i, j) = (Cnum(i, j) - sv.sig(i)) / 2.0 / SMALL_PERTURBATION;
				sv.eps(j) = strain(j);
			}
			updateInternal(false, false);
			sv.Cep = Cnum;
		}
		else {
			updateInternal(false, true);
		}
	}

	// done
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::setTrialStrain(const Vector& strain, const Vector& rate) {
	return setTrialStrain(strain);
}

int MultiYieldSurfaceHardeningSoftening::setTrialStrainIncr(const Vector& strain)
{
	return setTrialStrain(strain);
}

int MultiYieldSurfaceHardeningSoftening::setTrialStrainIncr(const Vector& strain, const Vector& rate) {
	return setTrialStrain(strain);
}

void MultiYieldSurfaceHardeningSoftening::setTheSurfaces(YieldSurfacePackage* theSurfaces) { ys = *theSurfaces; }

	// return material info
NDMaterial* MultiYieldSurfaceHardeningSoftening::getCopy(void) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::getCopy(void) -> subclass responsibility\n";
	exit(-1);
	return 0;
}

NDMaterial* MultiYieldSurfaceHardeningSoftening::getCopy(const char* type) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::getCopy(char*) -> subclass responsibility\n";
	exit(-1);
	return 0;
}

const char* MultiYieldSurfaceHardeningSoftening::getType(void) const {

	if (nOrd == 6) {
		return "ThreeDimensional";
	}
	else if (nOrd == 3) {
		return "PlaneStrain";
	}
	else {
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::getType() - unknown type!\n";
		exit(-1);
		return 0;
	}
}

int MultiYieldSurfaceHardeningSoftening::getDataDriver(void) {

	int preference = 0;		// automatic backbone genration

	if (use_data_driven_surface) {
		preference = 1;		// passive: do not update once generated
		if (use_active_approach) {
			preference = 2;	// active: update on the fly using the data
		}
	}
	else if (use_user_custom_surface) {
		preference = -1;	// user custom surface generation
	}

	return preference;
}

int MultiYieldSurfaceHardeningSoftening::getOrder(void) const {

	if (nOrd != 3 && nOrd != 6)
	{
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::getOrder() - unknown order!\n";
		exit(-1);
		return 0;
	}
	else {
		return nOrd;
	}
}

Vector MultiYieldSurfaceHardeningSoftening::getState(void) {
	Vector mState(nOrd);

	return mState;
}

double MultiYieldSurfaceHardeningSoftening::getGref(void) { return Gref; }

double MultiYieldSurfaceHardeningSoftening::getPref(void) { return Pref; }

double MultiYieldSurfaceHardeningSoftening::getGmod(void) { return sv.Gmod; }

	// return stress & strain
const Vector& MultiYieldSurfaceHardeningSoftening::getStress(void) {

	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stress = gs.p;
	stress = sv.sig;

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getStress() -> " << stress; }

	// done
	return stress;
}

const Vector& MultiYieldSurfaceHardeningSoftening::getStrain(void) {

	auto& gs = tools::getGlobalStorage(nOrd);
	auto& strain = gs.u;
	strain = sv.eps;

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getStrain() -> " << strain; }

	// done
	return strain;
}

	// return the tangent
const Matrix& MultiYieldSurfaceHardeningSoftening::getTangent(void) {
	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stiff = gs.K;
	stiff.Zero();
	stiff = sv.Cep;

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getTangent() -> " << stiff; }

	return stiff;
}

const Matrix& MultiYieldSurfaceHardeningSoftening::getInitialTangent(void) {
	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stiff = gs.K0;
	stiff.Zero();
	stiff = Kref * TensorM::IIvol(nOrd) + 2 * Gref * TensorM::IIdev(nOrd);	// compute initial elastic modulus

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getInitialTangent() -> " << stiff; }

	return stiff;
}

	// return material response
Response* MultiYieldSurfaceHardeningSoftening::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;
	const char* matType = this->getType();

	output.tag("NdMaterialOutput");
	output.attr("matType", this->getClassType());
	output.attr("matTag", this->getTag());

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int MultiYieldSurfaceHardeningSoftening::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;
	case 3:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getState();
		return 0;
	default:
		return -1;
	}
}

	// parallel message passing
int MultiYieldSurfaceHardeningSoftening::sendSelf(int commitTag, Channel& theChannel) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::sendSelf() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::recvSelf() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

	// miscellaneous
void MultiYieldSurfaceHardeningSoftening::Print(OPS_Stream& s, int flag) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::Print() -> subclass responsibility\n";
	exit(-1);
}

int MultiYieldSurfaceHardeningSoftening::setParameter(const char** argv, int argc, Parameter& param) {
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	// check for material tag
	if (theMaterialTag == this->getTag()) {
		if (strcmp(argv[0], "shearModulus") == 0) {
			return param.addObject(10, this);
		}
		else if (strcmp(argv[0], "bulkModulus") == 0) {
			return param.addObject(11, this);
		}
		else if (strcmp(argv[0], "TNYS") == 0) {
			return param.addObject(12, this);
		}
		else if (strcmp(argv[0], "cohesion") == 0) {
			return param.addObject(13, this);
		}
		else if (strcmp(argv[0], "frictionAngle") == 0) {
			return param.addObject(14, this);
		}
		else if (strcmp(argv[0], "dilatancyAngle") == 0) {
			return param.addObject(15, this);
		}
		else if (strcmp(argv[0], "peakStrain") == 0) {
			return param.addObject(16, this);
		}
		else if (strcmp(argv[0], "referencePressure") == 0) {
			return param.addObject(17, this);
		}
		else if (strcmp(argv[0], "tangentType") == 0) {
			return param.addObject(109, this);
		}
		else if (strcmp(argv[0], "materialType") == 0) {
			return param.addObject(110, this);
		}
		else if (strcmp(argv[0], "timeFactor") == 0) {
			return param.addObject(111, this);
		}
	}

	// done
	return -1;
}

int MultiYieldSurfaceHardeningSoftening::updateParameter(int responseID, Information& info) {

	if (responseID == 10) {
		Gref = info.theDouble;
	}
	else if (responseID == 11) {
		Kref = info.theDouble;
	}
	else if (responseID == 12) {
		theData->setTNYS(info.theDouble);
	}
	else if (responseID == 13) {
		theData->setCohesion(info.theDouble);
	}
	else if (responseID == 14) {
		theData->setFrictionAngle(info.theDouble);
	}
	else if (responseID == 15) {
		theData->setDilatancyAngle(info.theDouble);
	}
	else if (responseID == 16) {
		theData->setPeakStrain(info.theDouble);
	}
	else if (responseID == 17) {
		theData->setPref(info.theDouble);
	}
	else if (responseID == 109) {
		if (info.theInt == 0) {
			use_numerical_tangent = false;
		}
		else {
			use_numerical_tangent = true;
		}
	}
	else if (responseID == 110) {
		if (strcmp(info.theString, "ThreeDimensional") == 0) {
			nOrd = 6;
		}
		else if (strcmp(info.theString, "PlaneStrain") == 0) {
			nOrd = 3;
		}
		else {
			opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::updateParameter() - unsupported material type = " << info.theString << "\n";
			opserr << "MultiYieldSurfaceHardeningSoftening materialType can be --> ThreeDimensional (1), PlaneStrain (0)\n";
			opserr << "-----> CONTINUE WITHOUT UPDATING THE MATERIAL TYPE <-----\n";
		}
	}
	else if (responseID == 111) {
		// set user defined current time increment
		// this is useful for rate dependency in implicit mode and for the implex model
		// when using arc length or displacement control methods, where the pseudo time step
		// is actually the load factor.
		// if when this variable is first set or when it is set before the first commit
		// we set the committed variable to the same value
		sv.dtime_n = info.theDouble;
		if (!sv.dtime_first_set) {
			sv.dtime_n_commit = sv.dtime_n;
			sv.dtime_first_set = true;
		}
		sv.dtime_is_user_defined = true;
	}

	// done
	return 0;
}

// Private methods
	// the get methods
double MultiYieldSurfaceHardeningSoftening::getK(void) { return sv.Kmod; }

double MultiYieldSurfaceHardeningSoftening::getG(void) { return sv.Gmod; }

double MultiYieldSurfaceHardeningSoftening::getSizeYS(const int num_ys) { return yieldFunction(Vector(nOrd), num_ys, true); }

const Vector MultiYieldSurfaceHardeningSoftening::getStressVector(void) { return sv.sig_commit; }

const Vector MultiYieldSurfaceHardeningSoftening::getStrainVector(void) { return sv.eps_commit; }

const Vector MultiYieldSurfaceHardeningSoftening::getPlasticStrainVector(void) { return sv.xs_commit; }

double MultiYieldSurfaceHardeningSoftening::getMeanStress(const Vector& stress) { 
	double pressure = 0;
	if (nOrd == 3) {
		pressure = (0.5 * (stress(0) + stress(1)));
	}
	else if (nOrd == 6) {
		pressure = (1. / 3. * (stress(0) + stress(1) + stress(2)));
	}
	return pressure;
}

Vector MultiYieldSurfaceHardeningSoftening::getStressDeviator(const Vector& stress) {
	Vector deviator = Vector(nOrd);
	double pressure = getMeanStress(stress);
	Vector kronecker = TensorM::I(nOrd);
	for (int i = 0; i < nOrd; i++) {
		deviator(i) = stress(i) - (pressure * kronecker(i));
	}
	return deviator;
}

Vector MultiYieldSurfaceHardeningSoftening::getShiftedDeviator(const Vector& stress, const int num_ys) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::yieldFunction() -> subclass responsibility\n";
	exit(-1);
	return Vector(6);
}

	// yield surface operations
double MultiYieldSurfaceHardeningSoftening::yieldFunction(const Vector& stress, const int num_ys, bool yield_stress = false) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::yieldFunction() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dF_dS(const Vector& stress, const int num_ys) {
	Vector zero(nOrd);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dF_dS() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dF_dA(const Vector& stress, const int num_ys) {
	Vector zero(nOrd);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dF_dA() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dH_dA(const Vector& stress, const int num_ys) {
	Vector zero(nOrd);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dH_dA() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dP_dS(const Vector& stress, const int num_ys) {
	Vector zero(nOrd);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dP_dS() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

	// closest point projection methods
double MultiYieldSurfaceHardeningSoftening::linearizedFlow(const double dlambda) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::linearizedFlow() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

double MultiYieldSurfaceHardeningSoftening::yieldf(const Vector& stress, const int num_ys, bool yield_stress) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::yieldf() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MultiYieldSurfaceHardeningSoftening::Qi(const Vector& stress, const int num_ys) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::Qi() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MultiYieldSurfaceHardeningSoftening::Di(const Vector& stress, const int num_ys) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::Di() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

	// material internal operations
void MultiYieldSurfaceHardeningSoftening::updateModuli(const Vector& stress) {
	
	// get parameters
	sv.Gmod = Gref;
	sv.Kmod = Kref;

	// update elastic modulus w.r.t. average pressure
	double Pavg = getMeanStress(stress);	// compute average pressure
	if (Pavg > 0.0) {						// prevent tensile stress update
		Pavg = 0.0;
	}
	sv.Gmod = Gref * pow(abs(Pavg / Pref), Modn);	// update shear modulus
	sv.Kmod = Kref * pow(abs(Pavg / Pref), Modn);	// update bulk modulus
	
	sv.Ce = sv.Kmod * TensorM::IIvol(nOrd) + 2 * sv.Gmod * TensorM::IIdev(nOrd);	// compute elastic modulus
}

void MultiYieldSurfaceHardeningSoftening::updateInternal(bool do_implex, bool do_tangent) {

	// get commited values
	sv.sig = sv.sig_commit;
	sv.xs = sv.xs_commit;
	ys.revertToLastCommit();

	// initialize step plastic multiplier
	sv.lambda = 0.0;
	Vector ST = Vector(nOrd);
	Vector eps_incr = Vector(nOrd);

	// compute the strain increment (independent variable)
	eps_incr = sv.eps - sv.eps_commit;
	
	if (materialStage == 0) {										// linear elastic material
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
		sv.Cep = sv.Ce;
	}
	else if (materialStage == 2) {									// nonlinear elastic material
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
		sv.Cep = sv.Ce;
	}
	else {
		// compute trial stress
		ST = sv.sig + TensorM::inner(sv.Ce, (eps_incr));		    // trial stres
		if (do_implex && use_implex) {
			// EXPLICIT step: do explicit extrapolation
			implExIntegration(ST, do_tangent);
		}
		else {														// elastoplastic material
			// IMPLICIT step: solve constitutive equations
			double curr_yf_value = yieldFunction(ST, ys.now());			// yield function value
			if (curr_yf_value < ABSOLUTE_TOLERANCE) {					// elastic un/re-loading
				sv.sig = ST;
				sv.Cep = sv.Ce;
			}
			else {														// plastic loading
				int converged = 0;
				if (solution_strategy == 0) {
					// do return-mapping in steps
					Vector edev_incr = getStressDeviator(eps_incr);		// strain increment deviator
					double edev_step = sqrt(3.0 * 0.5 * TensorM::dotdot(edev_incr, edev_incr));
					int num_step = fmax(1, ((edev_step / CPLANE_STRAIN_STEP) + 1));
					Vector sig_step = (TensorM::inner(sv.Ce, (eps_incr)) / num_step);
					for (int i = 0; i < num_step; i++) {
						ST = sv.sig + sig_step;
						converged += cuttingPlaneAlgorithm(ST, do_tangent);
					}
					// compute the elastoplastic tangent
					if (do_tangent) { computeElastoplasticTangent(ys.now(), sv.sig); }
					if (converged) { 
						if (beVerbose) { 
							opserr << "MultiYieldSurfaceHardeningSoftening::updateInternal() -> cutting-plane (CP) algorithm converged in "
								<< converged << " iterations!\n";
						}
					} else { opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::updateInternal() - cutting-plane (CP) algorithm failed!\n"; }
				}
				else if (solution_strategy == 1) {
					converged = closestPointProjection(ST, do_tangent);
					if (converged != 0) {
						if (beVerbose) {
							opserr << "MultiYieldSurfaceHardeningSoftening::updateInternal() -> closest-point-projection (CPP) algorithm converged in "
								<< converged << " iterations!\n";
						}
					} else { opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::updateInternal() - closest-point-projection (CPP) algorithm failed!\n"; }
				}
				else {
					opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::updateInternal() - unknown return mapping algorithm!\n";
					exit(-1);
				}
				if (converged == 0) {
					opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::updateInternal() - return-mapping algorithm did not converge!\n";
					exit(-1);
				}
			}
		} // END if-else (implex or implicit stress)
	} // END if-else (material stage: elastic or plastic)
} 

void MultiYieldSurfaceHardeningSoftening::updateFailureSurface(const Vector& stress) {

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::updateFailureSurface()\n"; }

	//DO YOU NEED ZETA OR SIJ BELOW????

	// initialize variables 
	Vector alpha = ys.getAlpha(ys.getTNYS());
	Vector nn = get_dF_dS(stress, ys.getTNYS());
	Vector zeta = getStressDeviator(stress);

	// direction
	Vector unitdir = nn.Normalize();

	// magnitude
	alpha = zeta - sqrt(2. / 3.) * getSizeYS(ys.getTNYS()) * unitdir;

	// update alpha
	ys.setAlpha(alpha, ys.getTNYS());
	//updateInnerYieldSurfaces(ys.getTNYS(), stress);

}

		// compute methods
void MultiYieldSurfaceHardeningSoftening::computeElastoplasticTangent(int num_ys, const Vector& stress) {

	Vector xi = get_dF_dA(stress, num_ys);
	Vector rr = get_dH_dA(stress, num_ys);
	Vector nn = get_dF_dS(stress, num_ys);
	Vector mm = get_dP_dS(stress, num_ys);

	double denominator = TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) - TensorM::dotdot(xi, rr);

	if (denominator < ABSOLUTE_TOLERANCE ){
		opserr << "\FATAL: MultiYieldSurfaceHardeningSoftening::computeElastoplasticTangent() - division by 0 while computing elastoplastic tangent!\n";
		opserr << "Denominator [nn * Ce * mm - xi * rr] = 0. \n";
		opserr << "nn * Ce * mm = " << TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) << "\n";
		opserr << "xi * rr      = " << TensorM::dotdot(xi, rr) << "\n\n";
		exit(-1);
	}

	sv.Cep = sv.Ce - ((TensorM::inner(TensorM::outer(TensorM::inner(sv.Ce, mm), nn), sv.Ce))/denominator);
}

		// solution strategies
int MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm(const Vector& sigma_trial, const bool do_tangent) {

	// convergence status
	int convergence = 0;

	// initialize rate tensors
	Vector nn(nOrd);      Vector mm(nOrd);      Vector xi(nOrd);      Vector rr(nOrd);
	Vector next_nn(nOrd); Vector next_mm(nOrd); Vector next_xi(nOrd); Vector next_rr(nOrd);

	// Algorithm 7.2 Prevost (1985)
	// Step 1 - initialize
	int maxIter = 500;
	int iteration_counter = 0;
	double next_yf_value = 0;
	double curr_yf_value = 0;
	double old_yf_value = ABSOLUTE_TOLERANCE;
	double dlambda = 0.0; double dlambda1 = 0.0; double dlambda2 = 0.0;
	sv.sig = sigma_trial; // trial stress
	while (iteration_counter < maxIter) {
		// Do stress relaxation on to the current yield surface
			// (move the derivates below here to increase efficiency in the expense of non/linear update of of surfaces)
		while (true) {
			if (iteration_counter >= maxIter) {
				opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - return-mapping on to the current surface has failed!\n";
				return 0;
			}
			// Step 2 - check the consistency condition
			curr_yf_value = yieldFunction(sv.sig, ys.now());
			if ((curr_yf_value < RELATIVE_TOLERANCE) || (abs(curr_yf_value / old_yf_value) < ABSOLUTE_TOLERANCE)) {
				old_yf_value = curr_yf_value;
				break; // end the while loop and GO TO STEP 6
			}
				// compute the derivatives 	
			nn = get_dF_dS(sv.sig, ys.now());		// (move these outside the while loop)
			mm = get_dP_dS(sv.sig, ys.now());		// (they are constant for a surface anyway)
			xi = get_dF_dA(sv.sig, ys.now());
			rr = get_dH_dA(sv.sig, ys.now());
			// Step 3 - compute new plastic loading functions
			dlambda = 0.0;
			double denominator = TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) - TensorM::dotdot(xi, rr);
				// do a check
			if (abs(denominator) < ABSOLUTE_TOLERANCE) {
				opserr << "\nWARNING: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - division by 0 while computing the plastic loading function!\n";
				opserr << "Denominator [nn * Ce * mm - xi * rr] = " << denominator;
				opserr << "nn * Ce * mm = " << TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) << "\n";
				opserr << "xi * rr      = " << TensorM::dotdot(xi, rr) << "\n\n";
				return 0;
			}
			dlambda = curr_yf_value / denominator;
			sv.lambda = sv.lambda + dlambda;
			// Step 4 - update palstic strains and state variables
				// update plastic strain
			sv.xs = sv.xs + dlambda * mm;
				// update stress
			sv.sig = sv.sig - dlambda * TensorM::inner(sv.Ce, mm);
				// update current yield surface
			Vector alpha = ys.getAlpha(ys.now());
			alpha = alpha + dlambda * rr;	// update the center
			ys.setAlpha(alpha, ys.now());
				// update inner yield surfaces
			double yield_stress = getSizeYS(ys.now());
			if (yield_stress > ABSOLUTE_TOLERANCE) {
				Vector zeta = getShiftedDeviator(sv.sig, ys.now());
				for (int i = 0; i < ys.getNYS(); ++i) {
					Vector inner_zeta = getShiftedDeviator(sv.sig, i);
					Vector inner_alpha = ys.getAlpha(i);
					double inner_radius = getSizeYS(i);
					inner_alpha = inner_alpha + dlambda * (inner_zeta - inner_radius / yield_stress * zeta);
					ys.setAlpha(inner_alpha, i);
				}
			}
			else {
				opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - division by 0 while updating the internal surfaces!\n";
				return 0;
			}
			// Step 5 - iteration_counter++ and go to step 2
			old_yf_value = curr_yf_value;
			iteration_counter++;
		}
		// Step 6 - check for overshooting of the next yield surface
		next_yf_value = yieldFunction(sv.sig, ys.next());					// next yf_value
		double curr_H_prime = TensorM::dotdot((-1 * xi), rr);				// current H'
		int sign_H_prime = (curr_H_prime > 0) ? 1 : ((curr_H_prime < 0) ? -1 : 0);
		if ((ys.now() >= ys.getTNYS()) || ((next_yf_value * sign_H_prime) < ABSOLUTE_TOLERANCE)) {
			convergence = 1;
			break; // end while loop: algorithm is done...
		}
		// Do stress relaxation on to the next yield surface
			// compute the derivatives
		next_nn = get_dF_dS(sv.sig, ys.next());
		next_mm = get_dP_dS(sv.sig, ys.next());
		next_xi = get_dF_dA(sv.sig, ys.next());
		next_rr = get_dH_dA(sv.sig, ys.next());
		// Step 7 - correct the plastic loading function
		dlambda1 = dlambda2 = 0.0;
			// compute lambda1
		double curr_H0 = TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm);	// current H0
				// do a check
		if (curr_H_prime == 0) {
			opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - division by 0 while correcting the plastic multiplier (lambda 1)!\n";
			exit(-1);
		}
		dlambda1 = next_yf_value / curr_H_prime;
			// compute lambda2
		double next_H0 = TensorM::dotdot(TensorM::inner(next_nn, sv.Ce), next_mm);	// next H0
		double next_H_prime = TensorM::dotdot((-1 * next_xi), next_rr);				// next H'
				// do a check
		if ((next_H0 + next_H_prime) == 0) {
			opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - division by 0 while correcting the plastic multiplier (lambda 2)!\n";
			exit(-1);
		}
		dlambda2 = (next_yf_value / (next_H0 + next_H_prime)) * (1 + (curr_H0 / curr_H_prime));
		sv.lambda = sv.lambda + (dlambda2 - dlambda1);
		// Step 8 - correct plastic strains and state variables
			// correct plastic strain
		sv.xs = sv.xs - dlambda1 * mm + dlambda2 * next_mm;
			// correct stress
		sv.sig = sv.sig + dlambda1 * TensorM::inner(sv.Ce, mm) - dlambda2 * TensorM::inner(sv.Ce, next_mm);
		// Step 9 - increment number of active tield surface, iteration_counter++ and go to step 2
		old_yf_value = curr_yf_value;
		curr_yf_value = next_yf_value;
		ys.increment();
		iteration_counter++;
	}

	// if converged, return number of iterations instead
	if (convergence) { convergence = fmax(iteration_counter, 1); }
	return convergence;
}

int MultiYieldSurfaceHardeningSoftening::closestPointProjection(const Vector& sigma_trial, const bool do_tangent) {
	// Formulation by Gu et al. (2011)

	// convergence status
	int convergence = 0; int maxIter = 500;	int iteration_counter = 0;

	// initialize trial state
	Vector tau_trial, zeta_trial, curr_alpha, next_alpha, zeta_star, tau_star;

	// initialize rate tensors
	Vector nn(nOrd); Vector mm(nOrd);

	// initialize tangent tensors
	Matrix dTtrdE(nOrd, nOrd); Matrix dStrdE(nOrd, nOrd);

	// initialize return-mapping variables
	double next_yf_value = 0;
	double dlambda = 0.0;		double dlambda_bar = 0.0;
	double curr_H_prime = 0;	double old_H_prime = 0;
	double curr_K = 0.0;		double next_K = 0.0;
	double coeff_1 = 0.0;		double coeff_2 = 0.0;

	// set trial stress
	sv.sig = sigma_trial;
	if (do_tangent) { dTtrdE = 2 * sv.Gmod * TensorM::IIdev(nOrd); dStrdE = sv.Kmod * TensorM::IIvol(nOrd); }

	while (iteration_counter < maxIter) {

		// get trial stress deviator
		tau_trial = getStressDeviator(sv.sig);
		zeta_trial = getShiftedDeviator(tau_trial, ys.now());
		opserr << "trial_yf_value = " << yieldFunction(sv.sig, ys.now()) << "\n";
		opserr << "sv.sig     = " << sv.sig;
		opserr << "tau_trial  = " << tau_trial;
		opserr << "zeta_trial = " << zeta_trial;

		// compute contact stress
		curr_alpha = ys.getAlpha(ys.now());
		next_alpha = ys.getAlpha(ys.next());
		double Ki = sqrt((3.0 / 2.0) * TensorM::dotdot(zeta_trial, zeta_trial));	// eqn. 10
		curr_K = fmax((sqrt(3.0 / 2.0) * ys.getTau(ys.now())), SMALL_VALUE);		// current radius
		next_K = fmax((sqrt(3.0 / 2.0) * ys.getTau(ys.next())), SMALL_VALUE);		// next radius
		zeta_star = ((curr_K / Ki) * zeta_trial);							// eqn. 9
		tau_star = zeta_star + curr_alpha;
		opserr << "zeta_star = " << zeta_star;

		// compute the derivatives
		nn = zeta_star / sqrt(TensorM::dotdot(zeta_star, zeta_star));
		mm = zeta_star / sqrt(TensorM::dotdot(zeta_star, zeta_star));

		// compute plastic modulus H'
		sv.Hmod = ys.getEta(ys.now());
		if (abs(2 * sv.Gmod - sv.Hmod) < SMALL_VALUE) {
			curr_H_prime = LARGE_VALUE;
		}
		else {
			curr_H_prime = (2 * sv.Gmod * sv.Hmod) / (2 * sv.Gmod - sv.Hmod);
		}

		// compute iteration plastic multiplier
		coeff_1 = 1.0 / (curr_H_prime + 2 * sv.Gmod);
		if (iteration_counter == 0) {
			coeff_2 = 1.0;
		}
		else {
			double old_Hmod = ys.getEta(ys.prev());
			old_H_prime = (2 * sv.Gmod * old_Hmod) / (2 * sv.Gmod - old_Hmod);
			coeff_2 = (old_H_prime - curr_H_prime) / (old_H_prime);
		}
		dlambda = coeff_1 * coeff_2 * TensorM::dotdot(nn, (tau_trial - tau_star));
		opserr << "coeff_1 = " << coeff_1 << " coeff_2 = " << coeff_2 << "\n";
		opserr << "dlambda = " << dlambda << "\n";

		// update material internal variables
		sv.xs = sv.xs + dlambda * mm;
		sv.sig = sv.sig - dlambda * TensorM::inner(sv.Ce, mm);

		// compute the consistent tangent derivatives
		if (do_tangent) {
			Vector dKdE = (3.0 / (2.0 * Ki)) * TensorM::inner(zeta_trial, dTtrdE); //covariant
			Matrix dTstdE = (curr_K / Ki) * dTtrdE - (curr_K / (Ki * Ki)) * TensorM::outer(zeta_trial, dKdE); //mixed-variant (contra-co)
			Matrix dQdE = 1.0 / sqrt(TensorM::dotdot(zeta_star, zeta_star)) * dTstdE - TensorM::inner(TensorM::outer(zeta_star, zeta_star), dTstdE)
				/ sqrt(pow(TensorM::dotdot(zeta_star, zeta_star), 3)); //mixed-variant (contra-co)
			Vector dDlbdE = coeff_1 * coeff_2 * (TensorM::inner((tau_trial - tau_star), dQdE) + TensorM::inner(nn, (dTtrdE - dTstdE))); //covariant
			Matrix dTpldE = 2 * sv.Gmod * (TensorM::outer(dDlbdE, mm) + dlambda * dQdE);
			// correct stress derivative
			dTtrdE = dTtrdE - TensorM::inner(dTpldE, TensorM::IIdev(nOrd));

		}

		// check if overshhoting the next yield surface
		next_yf_value = yieldFunction(sv.sig, ys.next());			// next yf_value
		if ((ys.now() >= ys.getTNYS()) || (next_yf_value < ABSOLUTE_TOLERANCE)) {
			convergence = 1;
			break; // end while loop: algorithm is done...
		}
		// increment number of yield surface
		ys.increment();
		iteration_counter++;
	}
	
	// compute the direction of translation
	Vector roots(2);
	tau_trial = getStressDeviator(sv.sig);
	zeta_trial = getShiftedDeviator(tau_trial, ys.now());
	double A = TensorM::dotdot(tau_trial, tau_trial);
	double B = 2 * TensorM::dotdot((curr_alpha - next_alpha), (tau_trial - curr_alpha));
	double C = TensorM::dotdot((curr_alpha - next_alpha), (curr_alpha - next_alpha)) - 2.0 / 3.0 * next_K * next_K;
	roots = TensorM::quadratic(A, B, C);
	double ksi = fmax(roots(0), roots(1));
	Vector tau_trans = curr_alpha + ksi * (tau_trial - curr_alpha);
	Vector mu_alpha = (tau_trans - curr_alpha) - (curr_K / next_K) * (tau_trans - next_alpha);
	//mu_alpha = getStressDeviator(mu_alpha);
	//mu_alpha = mu_alpha / sqrt(TensorM::dotdot(mu_alpha, mu_alpha));	// unit vector

	// update rephrased variables (IMPL-EX: implicit correction)
	sv.lambda = sv.lambda + sqrt(3.0 / 2.0) / sqrt(TensorM::dotdot(zeta_star, zeta_star)) * dlambda;
	sv.ksi = sv.ksi + curr_H_prime / (TensorM::dotdot(nn, mu_alpha)) * dlambda;

	// compute translation magnitude
	double D = TensorM::dotdot(mu_alpha, mu_alpha);
	double E = -2 * TensorM::dotdot(mu_alpha, (tau_trial - curr_alpha));
	double F = TensorM::dotdot((tau_trial - curr_alpha), (tau_trial - curr_alpha)) - 2.0 / 3.0 * curr_K * curr_K;
	roots = TensorM::quadratic(D, E, F);
	double ksi_mag = fmin(roots(0), roots(1));

	// update current yield surface
	curr_alpha += ksi_mag * mu_alpha;	// update the center
	ys.setAlpha(curr_alpha, ys.now());

	// update inner yield surfaces
	for (int i = 0; i < ys.getNYS(); i++) {
		double inner_K = sqrt(3.0 / 2.0) * ys.getTau(i);
		Vector inner_alpha = ys.getAlpha(i);
		inner_alpha = tau_trial - (inner_K * zeta_trial) / curr_K;
		ys.setAlpha(inner_alpha, i);
	}

	// update consistent tangent
	if (do_tangent) {
		sv.Cep = dTtrdE + dStrdE;
	}

	opserr << "current_yf_value = " << yieldFunction(sv.sig, ys.now()) << "\n";
	opserr << "next_yf_value    = " << next_yf_value << "\n";

	// if converged, return number of iterations instead
	if (convergence) {convergence = ++iteration_counter; }
	return convergence;
}

int MultiYieldSurfaceHardeningSoftening::implExIntegration(const Vector& sigma_trial, const bool do_tangent) {

	// set trial stress 
	sv.sig = sigma_trial;

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if ((sv.dtime_n_commit > 0.0))
		time_factor = sv.dtime_n / sv.dtime_n_commit;
	// note: the implex method just wants the ratio of the new to the old time step
	// not the real time step, so it is just fine to assume it to be 1.
	// otherwise we have to deal with the problem of the opensees pseudo-time step
	// being the load multiplier in continuation methods...
	time_factor = 1.0;

	// extrapolate material internal variables
	double dlambda = (sv.lambda_commit - sv.lambda_commit_old) * time_factor;
	double dksi = (sv.ksi_commit - sv.ksi_commit_old) * time_factor;
	/* nYs = nYs_commit -> this step is done inside ys.revertToLastCommit() function which
							is called at the beginning of the updateInternal() method...   */
	sv.lambda = sv.lambda_commit + dlambda;
	sv.ksi = sv.ksi_commit + dksi;

	// get trial stress deviator
	Vector tau_trial = getStressDeviator(sv.sig);
	Vector zeta_trial = getShiftedDeviator(tau_trial, ys.now());
	//double dilatancy = ys.getBeta(ys.now());				    // dilatancy parameter

	// compute the derivatives
	Vector nn = zeta_trial / sqrt(TensorM::dotdot(zeta_trial, zeta_trial));
	Vector mm = zeta_trial / sqrt(TensorM::dotdot(zeta_trial, zeta_trial));

	// calculate material state using explicit integration
	sv.xs = sv.xs + dlambda * mm;
	sv.sig = sv.sig - dlambda * TensorM::inner(sv.Ce, mm);

	// compute the translation direction
	double curr_K = fmax((sqrt(3.0 / 2.0) * ys.getTau(ys.now())), SMALL_VALUE);			// current radius
	double next_K = fmax((sqrt(3.0 / 2.0) * ys.getTau(ys.next())), SMALL_VALUE);		// next radius
	Vector curr_alpha = ys.getAlpha(ys.now());											// current backstress
	Vector next_alpha = ys.getAlpha(ys.next());											// next backstress
	Vector dir = next_K / curr_K * (tau_trial - curr_alpha) - (tau_trial - next_alpha);	// direction vector

	// update current yield surface
	curr_alpha = curr_alpha + dksi * dir;	// update the center
	ys.setAlpha(curr_alpha, ys.now());

	// update internal surfaces
	for (int i = 0; i < ys.getNYS(); ++i) {
		double inner_K = sqrt(3.0 / 2.0) * ys.getTau(ys.now());
		Vector inner_alpha = ys.getAlpha(i);
		inner_alpha = tau_trial - (inner_K * zeta_trial) / curr_K;
		ys.setAlpha(inner_alpha, i);
	}

	// compute implex tangent
	//if (do_tangent) {
	//	Vector dKdE = (3.0 / (2.0 * Ki)) * TensorM::inner(zeta_trial, dTtrdE); //covariant
	//	Matrix dTstdE = (curr_K / Ki) * dTtrdE - (curr_K / (Ki * Ki)) * TensorM::outer(zeta_trial, dKdE); //mixed-variant (contra-co)
	//	Matrix dQdE = 1.0 / sqrt(TensorM::dotdot(zeta_star, zeta_star)) * dTstdE - TensorM::inner(TensorM::outer(zeta_star, zeta_star), dTstdE)
	//		/ sqrt(pow(TensorM::dotdot(zeta_star, zeta_star), 3)); //mixed-variant (contra-co)
	//	Vector dDlbdE = coeff_1 * coeff_2 * (TensorM::inner((tau_trial - tau_star), dQdE) + TensorM::inner(nn, (dTtrdE - dTstdE))); //covariant
	//	Matrix dTpldE = 2 * sv.Gmod * (TensorM::outer(dDlbdE, mm) + dlambda * dQdE);
	//	// correct stress derivative
	//	dTtrdE = dTtrdE - TensorM::inner(dTpldE, TensorM::IIdev(nOrd));
	//}

	return 0;
}

		// root search algorithm
double MultiYieldSurfaceHardeningSoftening::zbrentStress(const Vector& start_stress, const Vector& end_stress, 
	const int num_ys, const double x1, const double x2, const double tol)
{
	opserr << "MultiYieldSurfaceHardeningSoftening::zbrentStress::start_stress = " << start_stress << "\n";
	opserr << "MultiYieldSurfaceHardeningSoftening::zbrentStress::end_stress   = " << end_stress << "\n";
	int iter;
	double a = x1;
	double b = x2;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	double min1 = 0.0;
	double min2 = 0.0;
	double fc = 0.0;
	double p = 0.0;
	double q = 0.0;
	double r = 0.0;
	double s = 0.0;
	double tol1 = 0.0;
	double xm = 0.0;

	static Vector sigma_a(nOrd);
	static Vector sigma_b(nOrd);
	sigma_a = start_stress * (1 - a) + end_stress * a;
	sigma_b = start_stress * (1 - b) + end_stress * b;
	double fa = yieldFunction(sigma_a, num_ys);
	double fb = yieldFunction(sigma_b, num_ys);
	if ((fb * fa) > 0.0) {
		if (fabs(fa) < 1E-3) {
			fa = -fabs(fa);
		}
	}

	fc = fb;
	for (iter = 1; iter <= BRENT_MAXITER; iter++) {
		if ((fb * fc) > 0.0) {
			c = a;
			fc = fa;
			e = d = b - a;
		}

		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.0 * MACHINE_EPSILON * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if (fabs(xm) <= tol1 || fb == 0.0) {
			return b;
		}

		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}

			if (p > 0.0) {
				q = -q;
			}

			p = fabs(p);
			min1 = 3.0 * xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			}
			else {
				d = xm;
				e = d;
			}
		}
		else {
			d = xm;
			e = d;
		}

		a = b;
		fa = fb;
		if (fabs(d) > tol1) {
			b += d;
		}
		else {
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		}

		sigma_b = start_stress * (1 - b) + end_stress * b;
		fb = yieldFunction(sigma_b, num_ys);
	}

	//done
	return 0.0;
}