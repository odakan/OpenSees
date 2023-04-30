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
	double r0, double K0, double G0, double P0, double m0,
	DataDrivenNestedSurfaces* data, int ddtype, int itype, bool verbosity)
	:NDMaterial(tag, classTag), rho(r0), Kref(K0), Gref(G0), 
	Pref(P0), Modn(m0), theData(data), beVerbose(verbosity)
{
	
	// handle solution options
	if (itype == 1) {
		use_implex = true;
	}
	else {
		use_implex = false;
	}

	if (ddtype == 1) {						// offline
		use_online_approach = false;
		use_data_driven_surface = true;
	}
	else if (ddtype > 0 && ddtype < 1) {	// online
		use_online_approach = true;
		use_data_driven_surface = true;
	}
	else if (ddtype == 0) {					// automatic surface
		use_online_approach = false;
		use_data_driven_surface = false;
		use_user_custom_surface = false;
	}
	else if (ddtype == -1) {				// user custom surface
		use_online_approach = false;
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
		opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::::setTrialStrain() - material type is: " << getType() << "\n";
		opserr << "But the strain vector size is: " << strain.Size() << "\n";
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
		preference = 1;		// offline: do not update once generated
		if (use_online_approach) {
			preference = 2;	// online: update on the fly using the data
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
	Vector mState(6);

	return mState;
}

double MultiYieldSurfaceHardeningSoftening::getGref(void) { return Gref; }
double MultiYieldSurfaceHardeningSoftening::getPref(void) { return Pref; }

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

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getStress() -> " << stress; }

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

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getStrain() -> " << strain; }

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

	if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::getTangent() -> " << stiff; }

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
	}

	// done
	return -1;
}

int MultiYieldSurfaceHardeningSoftening::updateParameter(int responseID, Information& info) {

	if (responseID == 1) {
		// update material stage and if nonlinear set up yield surfaces 
		if (info.theInt > 0) {
			// set-up yield surfaces
			if (theData == nullptr) {
				opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::updateParameter() - nDMaterial " << getTag() <<
					" could not access the yield surface database!";
				exit(-1);
			}
			else {
				if (theData->isAOK(getDataDriver())) {
					ys = theData->generateYieldSurfaces(getTag(), getDataDriver(), Gref, Pref, Modn);
					materialStage = info.theInt;
				}
				else {
					opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::updateParameter() - nDMaterial " << getTag() << 
						" could not update material stage! Keeping the current stage = " << materialStage << " ...\n";
				}
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
double MultiYieldSurfaceHardeningSoftening::getSizeYS(const int num_ys) { return yieldFunction(Vector(6), num_ys, true); }
const Vector MultiYieldSurfaceHardeningSoftening::getStressVector(void) { return sv.sig_commit; }
const Vector MultiYieldSurfaceHardeningSoftening::getStrainVector(void) { return sv.eps_commit; }
const Vector MultiYieldSurfaceHardeningSoftening::getPlasticStrainVector(void) { return sv.xs_commit; }
double MultiYieldSurfaceHardeningSoftening::getMeanStress(const Vector& stress) { return (1. / 3. * (stress(0) + stress(1) + stress(2))); }
Vector MultiYieldSurfaceHardeningSoftening::getStressDeviator(const Vector& stress) {return stress - (getMeanStress(stress) * TensorM::I(6));}

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
	Vector zero(6);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dF_dS() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dF_dA(const Vector& stress, const int num_ys) {
	Vector zero(6);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dF_dA() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dH_dA(const Vector& stress, const int num_ys) {
	Vector zero(6);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dH_dA() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

Vector MultiYieldSurfaceHardeningSoftening::get_dP_dS(const Vector& stress, const int num_ys) {
	Vector zero(6);
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::get_dP_dS() -> subclass responsibility\n";
	exit(-1);
	return zero;
}

	// material internal operations
void MultiYieldSurfaceHardeningSoftening::updateModulus(const Vector& stress, const int num_ys) {
	
	// get parameters
	double Href = ys.getEta(num_ys);
	sv.Gmod = Gref;
	sv.Kmod = Kref;
	sv.Hmod = Href;

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
		sv.Gmod = Gref * pow(abs(Pavg / Pref), Modn);	// update shear modulus
		sv.Kmod = Kref * pow(abs(Pavg / Pref), Modn);	// update bulk modulus
		sv.Hmod = Href * pow(abs(Pavg / Pref), Modn);	// update plastic shear modulus
	}
	
	sv.Ce = sv.Kmod * TensorM::IIvol(6) + 2 * sv.Gmod * TensorM::IIdev(6);	// compute elastic modulus
}

void MultiYieldSurfaceHardeningSoftening::updateInternal(bool do_implex, bool do_tangent) {

	// get commited values
	sv.sig = sv.sig_commit;
	sv.xs = sv.xs_commit;
	ys.revertToLastCommit();

	// initialize step plastic multiplier
	sv.lambda = 0.0;
	
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

	if (materialStage == 0) {										// linear elastic material
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
		sv.Cep = sv.Ce;
	}
	else if (materialStage == 2) {									// nonlinear elastic material
		updateModulus(sv.sig, ys.now());								// update elastic modulus
		sv.sig = sv.sig + TensorM::inner(sv.Ce, (eps_incr));
		sv.Cep = sv.Ce;
	}
	else {					
		if (do_implex && use_implex) {
			// EXPLICIT step: do explicit extrapolation
			//  DO implex ...
			//  ...

			/* ys.num = ys.num_commit -> this step is done inside ys.revertToLastCommit() function which
			                             is called at the beginning of the updateInternal() method...   */
			sv.lambda = sv.lambda_commit + (sv.lambda_commit - sv.lambda_commit_old) * time_factor;



		}
		else {														// elastoplastic material	
			// IMPLICIT step: solve constitutive equations
			updateModulus(sv.sig, ys.now());							// update elastic modulus
			Vector ST = sv.sig + TensorM::inner(sv.Ce, (eps_incr));		// trial stress
			double curr_yf_value = yieldFunction(ST, ys.now());			// yield function value
			if (curr_yf_value < ABSOLUTE_TOLERANCE) {					// elastic un/re-loading
				sv.sig = ST;
			}
			else {														// plastic loading
				int converged = 0;
				if (solution_strategy == 0) {
					Vector edev_incr = getStressDeviator(eps_incr);		// strain increment deviator
					double edev_step = sqrt(3.0 * 0.5 * TensorM::dotdot(edev_incr, edev_incr));
					int num_step = fmax(1, ((edev_step / CPLANE_STRAIN_STEP) + 1));
					Vector sig_step = (TensorM::inner(sv.Ce, (eps_incr)) / num_step);
					for (int i = 0; i < num_step; i++) {
						ST = sv.sig + sig_step;
						converged = cuttingPlaneAlgorithm(ST, do_tangent);
					}
					// compute the elastoplastic tangent
					if (do_tangent) { computeElastoplasticTangent(ys.now(), sv.sig); }
				}
				else if (solution_strategy == 1) {
					converged = piecewiseLinearSolution(ST, do_tangent);
				}
				else {
					opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::updateInternal() - unknown return mapping algorithm!\n";
					exit(-1);
				}
				if (!converged) {
					opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::updateInternal() - return mapping algorithm did not converge!\n";
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
	Vector nn(6);      Vector mm(6);      Vector xi(6);      Vector rr(6);
	Vector next_nn(6); Vector next_mm(6); Vector next_xi(6); Vector next_rr(6);

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
		//opserr << "curr_H_prime = " << curr_H_prime << "\n";
		//opserr << "sign_H_prime = " << sign_H_prime << "\n";
		if ((ys.now() >= ys.getTNYS()) || ((next_yf_value * sign_H_prime) < ABSOLUTE_TOLERANCE)) {
			if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() -> return-mapping converged after " << iteration_counter << " iterations!\n"; }
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

	if (iteration_counter >= maxIter) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::cuttingPlaneAlgorithm() - return-mapping failed!\n";
		return 0;
	}

	return convergence;
}

int MultiYieldSurfaceHardeningSoftening::piecewiseLinearSolution(const Vector& sigma_trial, const bool do_tangent) {

	// convergence status
	int convergence = 0;

	// initialize rate tensors
	Vector nn(6); Vector mm(6); Vector rr(6);

	// Formulation by Gu et al. (2011)
	// Step 1 - initialize
	int maxIter = 500;
	int iteration_counter = 0;
	double next_yf_value = 0;
	double curr_H_prime = 0;
	double old_H_prime = 0;
	double dlambda = 0.0;
	double curr_K = 0.0;
	double next_K = 0.0;

	// set trial stress
	sv.sig = sigma_trial;

	while (iteration_counter < maxIter) {

		// get trial stress deviator
		Vector tau_trial = getStressDeviator(sv.sig);
		Vector zeta_trial = getShiftedDeviator(tau_trial, ys.now());

		// Differentiation of the elastic trial deviatoric stress with respect to the current deviatoric strain
		if (do_tangent) { Matrix dTtr_dE = 2 * sv.Gmod * TensorM::II4(6); }

		// compute contact stress
		double Ki = sqrt((3.0 / 2.0) * TensorM::dotdot(zeta_trial, zeta_trial));	// eqn. 10
		curr_K = ys.getTau(ys.now());												// current radius
		//////////////////////////////////////////////////////////////////////////////
		// Needs to be written generic for compatibility with other type of surfaces// 
		Vector alpha = ys.getAlpha(ys.now());
		Vector tau_star = ((curr_K / Ki) * zeta_trial) + alpha;	// eqn. 9
		//////////////////////////////////////////////////////////////////////////////
		Vector sig_star = tau_star + getMeanStress(sv.sig) * TensorM::I(6);

		// Differentiation of the contact stress with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// compute the derivatives
		nn = get_dF_dS(sig_star, ys.now());
		mm = get_dP_dS(sig_star, ys.now());
		rr = get_dH_dA(sig_star, ys.now());

		// Differentiation of the unit vector normal to the yield surface with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// compute lambda <L>/H' in eqn. 5
		sv.Hmod = ys.getEta(ys.now());
		old_H_prime = curr_H_prime;
		curr_H_prime = 1 / ((1 / sv.Hmod) - (0.5 * sv.Gmod));
		dlambda = fmax(TensorM::dotdot(nn, (tau_trial - tau_star)), 0) / curr_H_prime;

		// update plastic strain and plastic multiplier
		sv.lambda += dlambda;
		sv.xs += dlambda * mm;

		// Differentiation of the plastic multiplier with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// Differentiation of the plastic strain increment tensor with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// compute plastic stress eqn. 13
		Vector  SP = Vector(6);
		if (ys.now() == 0) {
			SP = ((2 * sv.Gmod * curr_H_prime) / (curr_H_prime + 2 * sv.Gmod)) * sv.xs;
		}
		else {
			SP = ((2 * sv.Gmod * curr_H_prime) / (curr_H_prime + 2 * sv.Gmod)) *
				((old_H_prime - curr_H_prime) / (old_H_prime)) * sv.xs;
		}

		// Differentiation of the plastic stress correction tensor with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// update stress
		sv.sig = sv.sig - SP;

		// Differentiation of the new trial stress after plastic correction with respect to the current deviatoric strain
		if (do_tangent) {

		}

		// check if overshhoting the next yield surface
		next_yf_value = yieldFunction(sv.sig, ys.next());			// next yf_value
		if ((ys.now() >= ys.getTNYS()) || (next_yf_value < ABSOLUTE_TOLERANCE)) {
			if (beVerbose) { opserr << "MultiYieldSurfaceHardeningSoftening::piecewiseLinearSolution() -> return-mapping converged after " << iteration_counter << " iterations!\n"; }
			convergence = 1;
			break; // end while loop: algorithm is done...
		}
		// increment number of yield surface
		ys.increment();
		iteration_counter++;
	}

	if (iteration_counter == maxIter) {
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::piecewiseLinearSolution() - return-mapping failed!\n";
		return 0;
	}

	// After convergence compute the direction of translation
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
		opserr << "WARNING: MultiYieldSurfaceHardeningSoftening::piecewiseLinearSolution() - division by 0 while updating the internal surfaces!\n";
		return 0;
	}

	// Compute the consistent tangent operator
	if (do_tangent) {

	}

	return convergence;
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

	static Vector sigma_a(6);
	static Vector sigma_b(6);
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