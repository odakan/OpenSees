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

// Constants
constexpr int BRENT_MAXITER = 20;
constexpr double BRENT_TOLERANCE = 1e-6;
constexpr double MACHINE_EPSILON = DBL_EPSILON;
constexpr double ABSOLUTE_TOLERANCE = 1e-4;
constexpr double RELATIVE_TOLERANCE = 1E-6;
constexpr double SMALL_PERTURBATION = 1.0e-9;

//	-- TO DO LIST --
// 
// 2) FIX compute elastoplastic tangent
//
// 3) FIX elastic tangent computation (Voigt notation 
//		compatible with elastoplastic one)
//
// 4) MAKE SURE all the necessary OpenSees functions are there
//
// 5) CODE all OpenSees related public functions
//
// 6) CODE ThreeDimensional, PlaneStrain tangent and stress

// Public methods
	// full constructor
MultiYieldSurfaceHardeningSoftening::MultiYieldSurfaceHardeningSoftening(int tag, int classTag,
	double r0, double K0, double G0, double P0, double m0, int T0,
	DataDrivenNestedSurfaces* data, int ddtype, int itype)
	:NDMaterial(tag, classTag), rho(r0), Kref(K0), Gref(G0), 
	Pref(P0), Modn(m0), TNYS(T0), theData(data)
{
	
	// handle solution options
	if (itype == 1) {
		use_implex = true;
	}
	else {
		use_implex = false;
	}

	if (ddtype == 2) {
		use_online_approach = true;
		use_data_driven_surface = true;
	}
	else if (ddtype == 1) {
		use_online_approach = false;
		use_data_driven_surface = true;
	}
	else {
		use_online_approach = false;
		use_data_driven_surface = false;
	}

	// do some checks
	if (theData == nullptr) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening:: yield surfaces are missing!\n";
		opserr << "The pointer that points at the yield surface: ys = " << theData << "\n";
		exit(-1);
	}
	else {
		// inform the yield surface object about the new instance
		theData->checkin(); // do check-in
	}
	if (nOrd != 3 && nOrd != 6) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening:: dimension error\n";
		opserr << "Material order has to be 3 or 6, but it is nOrd = " << nOrd << "\n";
		exit(-1);
	}
	if (Kref <= 0) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening: Kref <= 0\n";
		exit(-1);
	}
	if (Gref <= 0) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening: Gref <= 0\n";
		exit(-1);
	}
	if (Pref <= 0) {
		opserr << "WARNING:MultiYieldSurfaceHardeningSoftening: Pref <= 0\n";
		opserr << "Use Pref = 101.325 kPa?\n";
		exit(-1);
	}
	if (Modn < 0) {
		opserr << "WARNING:MultiYieldSurfaceHardeningSoftening: modn < 0\n";
		opserr << "modn is assumed to be zero...\n";
		Modn = 0.0;
	}
	if (rho < 0) {
		opserr << "WARNING:VonMisesDMM: mass density < 0\n";
		opserr << "mass density (rho) is assumed to be zero...\n";
		rho = 0.;
	}
	if (theData->getCohesion() < 0) {
		opserr << "FATAL:VonMisesDMM: cohesion <= 0\n";
		opserr << "Yield strength (cohesion) cannot be zero...\n";
		exit(-1);
	}	

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
		updateInternal(false, false);  // explicit_phase?, do_tangent?
	}

	// commit internal variables
	sv.eps_commit_old = sv.eps_commit;
	sv.eps_commit = sv.eps;
	sv.sig_commit_old = sv.sig_commit;
	sv.sig_commit = sv.sig;
	sv.xs_commit_old = sv.xs_commit;
	sv.xs_commit = sv.xs;
	if (ys.commitState()) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening::commitState: Yield surface object failed commit!\n";
		exit(-1);
	}

	if (DEBUG) { opserr << "MultiYieldSurfaceHardeningSoftening::commitState:\n"; sv.printStats(false); ys.printStats(false); }

	// done
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::revertToLastCommit(void) {

	if (DEBUG) { opserr << "MYSHS::revertToLastCommit! \n"; }

	// restore committed internal variables
	sv.eps = sv.eps_commit;
	sv.eps_commit = sv.eps_commit_old;
	sv.sig = sv.sig_commit;
	sv.sig_commit = sv.sig_commit_old;
	sv.xs = sv.xs_commit;
	sv.xs_commit = sv.xs_commit_old;
	if (ys.revertToLastCommit()) {
		opserr << "FATAL:MultiYieldSurfaceHardeningSoftening::revertToLastCommit: Yield surface object failed reverting to last commit!\n";
		exit(-1);
	}

	// done
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::revertToStart(void) {
	// reset state variables
	sv = MaterialStateVariables();
	ys = YieldSurfacePackage();
	materialStage = 0;

	// reset elastoplastic stiffness
	updateModulus(sv.sig, 0);
	sv.Cep = sv.Ce;

	// done
	return 0;
}

	// set material strain
int MultiYieldSurfaceHardeningSoftening::setTrialStrain(const Vector& strain) {
	
	// set the strain state
	if (nOrd == 3) { // PlaneStrain
		static Vector temp(6);
		temp(0) = strain(0);
		temp(1) = strain(1);
		temp(2) = 0.0;
		temp(3) = strain(2);
		temp(4) = 0.0;
		temp(5) = 0.0;
		sv.eps = temp;
	}
	else if (nOrd == 6) { // ThreeDimensional
		sv.eps = strain;
	}
	else {
		opserr << "MultiYieldSurfaceHardeningSoftening:: Material type is: " << getType() << "\n";
		opserr << "But strain vector size is: " << strain.Size() << "\n";
		exit(-1);
	}

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
		updateInternal(true, true);
		sv.sig_implex = sv.sig; // save stress for output
	}
	else {
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
	opserr << "MultiYieldSurfaceHardeningSoftening::getCopy(void) -- subclass responsibility\n";
	exit(-1);
	return 0;
}

NDMaterial* MultiYieldSurfaceHardeningSoftening::getCopy(const char* type) {
	opserr << "MultiYieldSurfaceHardeningSoftening::getCopy(const char* type) -- subclass responsibility\n";
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
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::getType -- unknown type!\n";
		exit(-1);
		return 0;
	}
}


int MultiYieldSurfaceHardeningSoftening::getDataDriver(void) {

	int preference = 0;  // automatic backbone genration

	if (use_data_driven_surface) {
		preference = 1;  // offline: do not update once generated
		if (use_online_approach) {
			preference = 2;  // online: update on the fly using the data
		}
	}

	return preference;
}


int MultiYieldSurfaceHardeningSoftening::getOrder(void) const {

	if (nOrd != 3 && nOrd != 6)
	{
		opserr << "MultiYieldSurfaceHardeningSoftening::getType -- unknown order!\n";
		exit(-1);
		return 0;
	}
	else {
		return nOrd;
	}
}

Vector MultiYieldSurfaceHardeningSoftening::getState(void) {
	Vector mState(6);

	return mState;
}

double MultiYieldSurfaceHardeningSoftening::getGref(void) { return Gref; }
double MultiYieldSurfaceHardeningSoftening::getPref(void) { return Pref; }
double MultiYieldSurfaceHardeningSoftening::getTNYS(void) { return TNYS; }

	// return stress & strain
const Vector& MultiYieldSurfaceHardeningSoftening::getStress(void) {

	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stress = gs.p;
	if (nOrd == 3) { //PlaneStrain
		stress(0) = sv.sig(0); stress(1) = sv.sig(1); stress(3) = sv.sig(3);
	}
	else { //ThreeDimensional
		stress = sv.sig;
	}

	if (DEBUG) { opserr << "MultiYieldSurfaceHardeningSoftening::getStress: " << stress; }

	// done
	return stress;
}

const Vector& MultiYieldSurfaceHardeningSoftening::getStrain(void) {

	auto& gs = tools::getGlobalStorage(nOrd);
	auto& strain = gs.u;
	if (nOrd == 3) { //PlaneStrain
		strain(0) = sv.eps(0); strain(1) = sv.eps(1); strain(3) = sv.eps(3);
	}
	else { //ThreeDimensional
		strain = sv.eps;
	}

	if (DEBUG) { opserr << "MultiYieldSurfaceHardeningSoftening::getStrain: " << strain; }

	// done
	return strain;
}

	// return the tangent
const Matrix& MultiYieldSurfaceHardeningSoftening::getTangent(void) {
	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stiff = gs.K;
	stiff.Zero();
	if (nOrd == 3) { //PlaneStrain
		stiff(0, 0) = sv.Cep(0, 0); stiff(0, 1) = sv.Cep(0, 1); stiff(0, 3) = sv.Cep(0, 3);
		stiff(1, 0) = sv.Cep(1, 0); stiff(1, 1) = sv.Cep(1, 1); stiff(1, 3) = sv.Cep(1, 3);
		stiff(3, 0) = sv.Cep(3, 0); stiff(3, 1) = sv.Cep(3, 1); stiff(3, 3) = sv.Cep(3, 3);
	}
	else { //ThreeDimensional
		stiff = sv.Cep;
	}

	if (DEBUG) { opserr << "MYSHS::getTangent Cep = " << stiff; }

	return stiff;
}

const Matrix& MultiYieldSurfaceHardeningSoftening::getInitialTangent(void) {
	auto& gs = tools::getGlobalStorage(nOrd);
	auto& stiff = gs.K0;
	auto& zeros = gs.p;
	stiff.Zero();
	zeros.Zero();
	sv.Ce = Kref * TensorM::IIvol(6) + 2 * Gref * TensorM::IIdev(6);	// compute initial elastic modulus
	if (nOrd == 3) { //PlaneStrain
		stiff(0, 0) = sv.Ce(0, 0); stiff(0, 1) = sv.Ce(0, 1); stiff(0, 3) = sv.Ce(0, 3);
		stiff(1, 0) = sv.Ce(1, 0); stiff(1, 1) = sv.Ce(1, 1); stiff(1, 3) = sv.Ce(1, 3);
		stiff(3, 0) = sv.Ce(3, 0); stiff(3, 1) = sv.Ce(3, 1); stiff(3, 3) = sv.Ce(3, 3);
	}
	else { //ThreeDimensional
		stiff = sv.Ce;
	}

	if (DEBUG) { opserr << "MYSHS::getInitialTangent Ce = " << stiff; }

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
	opserr << "MultiYieldSurfaceHardeningSoftening::sendSelf -> Subclass responsibility\n";
	exit(-1);
	return 0;
}

int MultiYieldSurfaceHardeningSoftening::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	opserr << "MultiYieldSurfaceHardeningSoftening::recvSelf -> Subclass responsibility\n";
	exit(-1);
	return 0;
}

	// miscellaneous
void MultiYieldSurfaceHardeningSoftening::Print(OPS_Stream& s, int flag) {
	opserr << "MultiYieldSurfaceHardeningSoftening::Print -> Subclass responsibility\n";
	exit(-1);
}

int MultiYieldSurfaceHardeningSoftening::setParameter(const char** argv, int argc, Parameter& param) {
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
		else if (strcmp(argv[0], "tangentType") == 0) {
			return param.addObject(109, this);
		}
		else if (strcmp(argv[0], "materialType") == 0) {
			return param.addObject(110, this);
		}
		else if (strcmp(argv[0], "timeFactor") == 0) {
			return param.addObject(111, this);
		}
		else if (strcmp(argv[0], "resetPlasticInternalVariables") == 0) {
			return param.addObject(112, this);
		}
	}

	// done
	return -1;
}

int MultiYieldSurfaceHardeningSoftening::updateParameter(int responseID, Information& info) {

	if (responseID == 1) {
		// update material stage and if nonlinear set up yield surfaces 
		materialStage = info.theInt;
		if (materialStage > 0) {
			// set-up yield surfaces
			if (theData == nullptr) {
				opserr << "FATAL:MultiYieldSurfaceHardeningSoftening: Could not access yield surface database!";
				exit(-1);
			}
			else {
				ys = theData->generateYieldSurfaces(getTag(), getDataDriver(), Pref, Gref, TNYS);
			}
		}
	}
	else if (responseID == 10) {
		Gref = info.theDouble;
	}
	else if (responseID == 11) {
		Kref = info.theDouble;
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
			opserr << "WARNING: MultiYieldSurfaceHardeningSoftening updateParameter: unsupported material type = " << info.theString << "\n";
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
double MultiYieldSurfaceHardeningSoftening::getSizeYS(const int num_yield_surface) { return yieldFunction(sv.sig, num_yield_surface, true); }
double MultiYieldSurfaceHardeningSoftening::getMeanStress(const Vector& stress) { return 1. / 3. * (stress(0) + stress(1) + stress(2)); }
const Vector MultiYieldSurfaceHardeningSoftening::getStressVector(void) { return sv.sig_commit; }
const Vector MultiYieldSurfaceHardeningSoftening::getStrainVector(void) { return sv.eps_commit; }
const Vector MultiYieldSurfaceHardeningSoftening::getPlasticStrainVector(void) { return sv.xs_commit; }
const Vector MultiYieldSurfaceHardeningSoftening::getStressDeviator(const Vector& stress, int num_yield_surface) {
	// compute the stress deviator and apply Ziegler's Rule (kinematic hardening)
	return stress - (getMeanStress(stress) * (TensorM::I(6) + ys.getAlpha(num_yield_surface, ys.getNYS_commit())));
	
}

	// yield surface operations
double MultiYieldSurfaceHardeningSoftening::yieldFunction(const Vector& stress, const int num_yield_surface, bool yield_stress = false) {
	opserr << "MultiYieldSurfaceHardeningSoftening::yieldFunction -> Subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dF_dS(const Vector& stress, const int num_yield_surface) {
	Vector zero(6);
	opserr << "MultiYieldSurfaceHardeningSoftening::get_dF_dS -> Subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dF_dA(const Vector& stress, const int num_yield_surface) {
	Vector zero(6);
	opserr << "MultiYieldSurfaceHardeningSoftening::get_dF_dS -> Subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dH_dA(const Vector& stress, const int num_yield_surface) {
	Vector zero(6);
	opserr << "MultiYieldSurfaceHardeningSoftening::get_dH_dA -> Subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dP_dS(const Vector& stress, const int num_yield_surface) {
	Vector zero(6);
	opserr << "MultiYieldSurfaceHardeningSoftening::plasticFlowDirect -> Subclass responsibility\n";
	exit(-1);
	return zero;
}

	// material internal operations
		// update methods
void MultiYieldSurfaceHardeningSoftening::updateStress(Vector& stress, const double lambda, const int num_yield_surface) {
	Vector mm = get_dP_dS(stress, num_yield_surface);	// get_dP_dS (plastic flow diraction)
	// Update stress
//updateModulus(stress, num_yield_surface);
	stress = stress - lambda * TensorM::inner(sv.Ce, mm);
}

void MultiYieldSurfaceHardeningSoftening::updateModulus(const Vector& stress, const int num_yield_surface) {
	// get parameters

	double Href = ys.getEta(num_yield_surface, ys.getNYS_commit());

	sv.Gmod = Gref;
	sv.Kmod = Kref;
	sv.Hmod = Href;

	/*
	// update elastic modulus w.r.t. average pressure
	if (materialStage == 0) { // do not update if elastic stage
		sv.Gmod = Gref;
		sv.Kmod = Kref;
		sv.Hmod = Href;
	}
	else {
		double Pavg = getMeanStress(stress);	// compute average pressure
		if (Pavg > 0.0) {						// prevent tensile stress update
			Pavg = 0.0;
		}
		sv.Gmod = Gref * pow(abs(Pavg / Pref), modn);	// update shear modulus
		sv.Kmod = Kref * pow(abs(Pavg / Pref), modn);	// update bulk modulus
		sv.Hmod = Href * pow(abs(Pavg / Pref), modn);	// update plastic shear modulus
	}
	*/
	sv.Ce = sv.Kmod * TensorM::IIvol(6) + 2 * sv.Gmod * TensorM::IIdev(6);	// compute elastic modulus
}

void MultiYieldSurfaceHardeningSoftening::updateInternal(bool do_implex, bool do_tangent) {
	// get commited values
	sv.sig = sv.sig_commit;
	sv.xs = sv.xs_commit;
	ys.revertToLastCommit();
	
	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (do_implex && use_implex && (sv.dtime_n_commit > 0.0))
		time_factor = sv.dtime_n / sv.dtime_n_commit;
	// note: the implex method just wants the ratio of the new to the old time step
	// not the real time step, so it is just fine to assume it to be 1.
	// otherwise we have to deal with the problem of the opensees pseudo-time step
	// being the load multiplier in continuation methods...
	time_factor = 1.0;

	// compute the strain increment (independent variable)
	Vector eps_incr = sv.eps - sv.eps_commit;

	if (materialStage == 0) {														// linear elastic material
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
	}
	else if (materialStage == 2) {													// nonlinear elastic material
		updateModulus(sv.sig, ys.getNYS());												// update elastic modulus
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
	}
	else {					
		if (do_implex && use_implex) {
			// EXPLICIT step: do explicit extrapolation
			//  DO implex ...
			//  ...




		}
		else {																		// elastoplastic material	
			// IMPLICIT step: solve constitutive equations
			updateModulus(sv.sig, ys.getNYS());										// update elastic modulus
			Vector ST = sv.sig + TensorM::inner(sv.Ce, (eps_incr));				// trial stress
			double curr_yf_value = yieldFunction(ST, ys.getNYS());			// yield function value
			if (curr_yf_value < ABSOLUTE_TOLERANCE) {							// elastic un/re-loading
				sv.sig = ST;
			}
			else {																// plastic loading
				int converged = -1;
				if (solution_strategy == 0) {
					converged = cuttingPlaneAlgorithm(ST);				// forward-euler
				}
				else if (solution_strategy == 1) {
					converged = closestPointProjection(ST);				// backward-euler
				}
				else {
					opserr << "FATAL:MultiYieldSurfaceHardeningSoftening::updateInternal: Unknown return mapping algorithm!\n";
					exit(-1);
				}
				if (converged != 0)
				{
					opserr << "FATAL:MultiYieldSurfaceHardeningSoftening::updateInternal: Return mapping algorithm did not converge!\n";
					exit(-1);
				}
			}
		} // END if-else (implex or implicit stress)
	} // END if-else (material stage: elastic or plastic)

	// compute the consistent tangent
	if (do_tangent) {
		if (use_implex) {
			sv.Cep.Zero();

		}
		else {
			if (materialStage != 1 ) {
				sv.Cep = sv.Ce;
			}
			else {
				computeElastoplasticTangent(ys.getNYS(), sv.sig);
			}
		} // END if-else (implex or implicit tangent)
	} // END if (do_tangent)
} 

void MultiYieldSurfaceHardeningSoftening::updatePlasticStrain(Vector& pstrain, const Vector& stress, const double lambda, const int num_yield_surface) {
	// Update plastic strain
	Vector mm = get_dP_dS(stress, num_yield_surface);	// get_dP_dS (plastic flow diraction)
	pstrain = pstrain + lambda * mm;					// add plastic strain increment
}

void MultiYieldSurfaceHardeningSoftening::updateFailureSurface(const Vector& stress) {
	//////////////////////////////
	if (DEBUG) { opserr << "MYSHS::updateFailureSurface in! \n"; }
	//////////////////////////////

	//DO YOU NEED ZETA OR SIJ BELOW????

	// initialize variables 
	Vector alpha = ys.getAlpha(TNYS, ys.getNYS_commit());
	Vector nn = get_dF_dS(stress, TNYS);
	Vector zeta = getStressDeviator(stress, TNYS);

	// direction
	Vector unitdir = nn.Normalize();

	// magnitude
	alpha = zeta - sqrt(2. / 3.) * getSizeYS(TNYS) * unitdir;

	// update alpha
	ys.setAlpha(alpha, TNYS);
	updateInnerYieldSurfaces(TNYS, stress);

}

void MultiYieldSurfaceHardeningSoftening::updateInnerYieldSurfaces(int num_yield_surface, const Vector& stress) {

	//DO YOU NEED ZETA OR SIJ BELOW????

	// Update the inner yield surfaces (alpha)
	Vector curr_alpha = ys.getAlpha(num_yield_surface, ys.getNYS_commit());
	double curr_radius = getSizeYS(num_yield_surface);

	if (curr_radius != 0) {
		Vector zeta = getStressDeviator(stress, num_yield_surface);
		for (int i = 0; i < num_yield_surface; ++i) {
			Vector the_alpha = ys.getAlpha(i, ys.getNYS_commit());
			double the_radius = getSizeYS(i);
			the_alpha = zeta - the_radius / curr_radius * zeta - curr_alpha;
			ys.setAlpha(the_alpha, i);
		}
	}
	else {
		if (curr_radius == 0) {
			opserr << "\n";
			opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::updateInnerYieldSurfaces()\n";
			opserr << "Denominator -- curr_radius = 0...\n\n";
		}
	}
}

void MultiYieldSurfaceHardeningSoftening::updateCurrentYieldSurface(int num_yield_surface, double lambda, const Vector& stress) {
	// Update the current active yield surface (alpha)
	Vector rr = get_dH_dA(stress, num_yield_surface);	// rate of alpha
	Vector alpha = ys.getAlpha(num_yield_surface, ys.getNYS_commit());
	alpha = alpha + lambda * rr;						// update the alpha
	ys.setAlpha(alpha, num_yield_surface);
}

		// compute methods
void MultiYieldSurfaceHardeningSoftening::computeElastoplasticTangent(int num_yield_surface, const Vector& stress) {
	Vector xi = get_dF_dA(stress, num_yield_surface);
	Vector rr = get_dH_dA(stress, num_yield_surface);
	Vector nn = get_dF_dS(stress, num_yield_surface);
	Vector mm = get_dP_dS(stress, num_yield_surface);

//updateModulus(stress, num_yield_surface);

	double denominator = TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) - TensorM::dotdot(xi, rr);

	if (denominator == 0 ){
	    opserr << "DruckerPragerMYSHS::computeElastoplasticTangent" << "\n";
		opserr << "Error denominator                          = 0 " << "\n";
		opserr << "LEFT  = nn(i,j) : sv.Ce(i,j,k,l) : mm(k,l) = " << TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) << "\n";
		opserr << "RIGHT = xi(o,t) : rr(o,t)                  = " << TensorM::dotdot(xi, rr) << "\n";
	}

	sv.Cep = sv.Ce - ((TensorM::inner(TensorM::outer(TensorM::inner(sv.Ce, mm), nn), sv.Ce))/denominator);
}

double MultiYieldSurfaceHardeningSoftening::computePlasticLoadingFunction(const Vector& stress, const double yf_value, const int num_yield_surface) {

	Vector nn = get_dF_dS(stress, num_yield_surface);
	Vector mm = get_dP_dS(stress, num_yield_surface);
	Vector xi = get_dF_dA(stress, num_yield_surface);
	Vector rr = get_dH_dA(stress, num_yield_surface);

	// get lambda
	double denominator = 0;
	denominator = TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) - TensorM::dotdot(xi, rr);
	if (denominator == 0) {
		opserr << "MultiYieldSurfaceHardeningSoftening::computePlasticLoadingFunction():\n";
		opserr << "Error denominator == 0 \n";
		opserr << "LEFT  = curr_nn * Ce * curr_mm == " << TensorM::dotdot(TensorM::inner(nn, sv.Ce), mm) << "\n";
		opserr << "RIGHT = curr_xi * bar_alpha    == " << TensorM::dotdot(xi, rr) << "\n";
	}
	return (yf_value / denominator);
}

		// solution strategies
int MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm(const Vector& sigma_trial) {

	// convergence status
	int converged = -1;

	// get some constants
	double curr_yf_value = 0;

	// Algorithm 7.2 Prevost (1985)
	// Step 1 - initialize
	int maxIter = 500;
	double old_yf_value = ABSOLUTE_TOLERANCE;
	int iteration_counter = 0;
	double lambda = 0.0;
	sv.sig = sigma_trial; // trial stress
	while (iteration_counter < maxIter) {
		// Do stress relaxation on the current yield surface
		while (iteration_counter < maxIter) {
			// Step 2 - update stress and check consistency condition
			curr_yf_value = yieldFunction(sv.sig, ys.getNYS());
			if ((curr_yf_value < RELATIVE_TOLERANCE) || (abs(curr_yf_value / old_yf_value) < ABSOLUTE_TOLERANCE)) {
				old_yf_value = curr_yf_value;
				break; // end the while loop and GO TO STEP 6
			}
			// Step 3 - compute new plastic loading functions
			lambda = computePlasticLoadingFunction(sv.sig, curr_yf_value, ys.getNYS());
			// Step 4 - update palstic strains and state variables
			updatePlasticStrain(sv.xs, sv.sig, lambda, ys.getNYS());
			updateStress(sv.sig, lambda, ys.getNYS());
			updateCurrentYieldSurface(ys.getNYS(), lambda, sv.sig);
			updateInnerYieldSurfaces(ys.getNYS(), sv.sig);
			// Step 5 - iteration_counter++ and go to step 2
			old_yf_value = curr_yf_value;
			iteration_counter++;
		}
		// Step 6 - check for overshooting of the next yield surface
		curr_yf_value = yieldFunction(sv.sig, ys.getNYS() + 1);	// next yf_value
		if ((ys.getNYS() == TNYS) || (curr_yf_value < ABSOLUTE_TOLERANCE)) {
			if (DEBUG) { opserr << "MYSHS::cuttingPlaneAlgorithm: cuting-plane solver converged after " << iteration_counter << " iterations!\n"; }
			converged = 0;
			break; // end while loop: algorithm is done...
		}
		// Step 7 - correct the plastic loading function
		double lambda1, lambda2;
		correctPlasticLoadingFunction(sv.sig, lambda1, lambda2, ys.getNYS());
		// Step 8 - correct plastic strains and state variables
		correctPlasticStrain(sv.xs, sv.sig, lambda1, lambda2, ys.getNYS());
		correctStress(sv.sig, lambda1, lambda2, ys.getNYS());
		// Step 9 - increment number of active tield surface, iteration_counter++ and go to step 2
		ys.incrementNYS();
		old_yf_value = curr_yf_value;
		iteration_counter++;
	}

	if (iteration_counter == maxIter)
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm: return mapping failed!\n";

	return converged;
}

int MultiYieldSurfaceHardeningSoftening::closestPointProjection(const Vector& eps_incr) {

return 0;
}

		// correct methods
void MultiYieldSurfaceHardeningSoftening::correctStress(Vector& stress, const double lambda1, const double lambda2, const int num_yield_surface) {
	//////////////////////////////
	if (DEBUG) { opserr << "MYSHS::correctUpdateStress in! \n"; }
	////////////////////////////

	// After overshooting, Correct the overshooting stress
	Vector curr_mm = get_dF_dS(stress, num_yield_surface);		// curr_nn
	Vector next_mm = get_dF_dS(stress, num_yield_surface + 1);	// next_nn

	// correct stress
//updateModulus(stress, num_yield_surface);
	stress = stress + lambda1 * TensorM::inner(sv.Ce, curr_mm) - lambda2 * TensorM::inner(sv.Ce, next_mm);

	//////////////////////////////
	if (DEBUG) { opserr << "MYSHS::correctUpdateStress out! \n"; }
	//////////////////////////////
}

void MultiYieldSurfaceHardeningSoftening::correctPlasticStrain(Vector& pstrain, const Vector& stress,
	const double lambda1, const double lambda2, const int num_yield_surface) {
	//////////////////////////////
	if (DEBUG) { opserr << "MYSHS::correctPlasticStrain in! \n"; }
	//////////////////////////////
	// After the overshooting, Correct the plastic strain
	Vector curr_mm = get_dP_dS(stress, num_yield_surface);			// curr_mm
	Vector next_mm = get_dP_dS(stress, num_yield_surface + 1);		// next_mm
	pstrain = pstrain - lambda1 * curr_mm + lambda2 * next_mm;		// lambda * plastic_flow
	//////////////////////////////
	if (DEBUG) { opserr << "MYSHS::correctPlasticStrain out! \n"; }
	//////////////////////////////
}

void MultiYieldSurfaceHardeningSoftening::correctPlasticLoadingFunction(Vector& stress, double& lambda1, double& lambda2, const int num_yield_surface) {

	// compute lambda1
	Vector curr_nn = get_dF_dS(stress, num_yield_surface);				// dF_dS
	Vector curr_xi = get_dF_dA(stress, num_yield_surface);				// dF_dA
	Vector curr_rr = get_dH_dA(stress, num_yield_surface);				// alpha_direction
	double curr_H_prime = TensorM::dotdot((-1 * curr_xi), curr_rr);		// curr_H_prime

	if (curr_H_prime == 0) {
		opserr << "\n";
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::correctPlasticLoadingFunction\n ";
		opserr << "Denominator --> curr_H_prime = 0\n ";
		exit(-1);
	}

	double next_yf_val = yieldFunction(stress, num_yield_surface + 1);
	lambda1 = next_yf_val / curr_H_prime;

	// compute lambda2
	Vector next_nn = get_dF_dS(stress, num_yield_surface + 1);			// next dF_dS
	Vector next_xi = get_dF_dA(stress, num_yield_surface + 1);			// next dF_dA
	Vector next_rr = get_dH_dA(stress, num_yield_surface + 1);			// next alpha_direction
	double next_H_prime = TensorM::dotdot((-1 * next_xi), next_rr);		// next_H_prime
	double next_H0 = TensorM::dotdot(TensorM::inner(next_nn, sv.Ce), next_nn);

	if ((next_H0 + next_H_prime) == 0) {
		opserr << "\n";
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::correctPlasticLoadingFunction\n ";
		opserr << "Denominator --> (next_H0 + next_H_prime) = 0\n ";
		exit(-1);
	}

	double numerator = TensorM::dotdot(TensorM::inner(next_nn, sv.Ce), next_nn);
	lambda2 = next_yf_val * (1 + numerator / curr_H_prime) / (next_H0 + next_H_prime);
}

		// root search algorithm
double MultiYieldSurfaceHardeningSoftening::zbrentStress(const Vector& start_stress, const Vector& end_stress, 
	const int num_yield_surface, const double x1, const double x2, const double tol)
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

	static Vector sigma_a(6);
	static Vector sigma_b(6);
	sigma_a = start_stress * (1 - a) + end_stress * a;
	sigma_b = start_stress * (1 - b) + end_stress * b;
	double fa = yieldFunction(sigma_a, num_yield_surface);
	double fb = yieldFunction(sigma_b, num_yield_surface);
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
		fb = yieldFunction(sigma_b, num_yield_surface);
	}

	//done
	return 0.0;
}