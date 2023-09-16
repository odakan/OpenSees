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

namespace global {
	class MaterialList {
	public:
		int size = 0;
		int next_tag = 0;
		std::vector<MatsuokaNakaiDMM*> mptr;

	public:
		MaterialList(void) = default;
		~MaterialList(void) = default;

		void append(MatsuokaNakaiDMM* matptr) {
			matptr->setSubTag(next_tag);
			next_tag++;
			mptr.push_back(matptr);
			size = mptr.size();
		}

		void remove(MatsuokaNakaiDMM* matptr) {
			auto it = std::find(mptr.begin(), mptr.end(), matptr);
			if (it != mptr.end()) {
				mptr.erase(it);
			}
			// clean up if the list is empty or if not update the size
			if (mptr.empty()) { cleanup(); }
			else { size = mptr.size(); }
		}

		void cleanup(void) {
			size = 0;
			next_tag = 0;
			mptr.clear();
		}
	};
}

// global material list
auto MatsuokaNakailist = global::MaterialList();

// Public methods
	// full constructor
MatsuokaNakaiDMM::MatsuokaNakaiDMM(int tag, double r0,
	double K0, double G0, double P0, double m0,
	std::shared_ptr<DataDrivenNestedSurfaces> ys,
	int ddtype, int itype, bool verbosity)
	:MultiYieldSurfaceHardeningSoftening(tag, ND_TAG_MatsuokaNakaiDMM, r0,
		K0, G0, P0, m0, ys, ddtype, itype, verbosity)
{
	// append material to the list
	MatsuokaNakailist.append(this);
}

// null constructor
MatsuokaNakaiDMM::MatsuokaNakaiDMM()
	:MultiYieldSurfaceHardeningSoftening()
{
}

// destructor
MatsuokaNakaiDMM::~MatsuokaNakaiDMM()
{
	MatsuokaNakailist.remove(this);
}

// return object info
NDMaterial* MatsuokaNakaiDMM::getCopy(void) {

	MatsuokaNakaiDMM* copy = nullptr;
	copy = new MatsuokaNakaiDMM(*this);

	if (copy != nullptr) {
		// inform the list about the new instance
		MatsuokaNakailist.append(copy);
		return copy;
	}
	else {
		opserr << "FATAL: MatsuokaNakaiDMM::getCopy() - returned NULLPTR! couldn't copy the material\n";
		exit(-1);
	}
	// done
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
		opserr << "FATAL: MatsuokaNakaiDMM::getCopy() - unsupported material type = " << type << "\n";
		opserr << "MatsuokaNakaiDMM materialType can be --> ThreeDimensional (1), PlaneStrain (0)\n";
		exit(-1);
	}
	// done

	if (copy != nullptr) {
		// inform the list about the new instance
		MatsuokaNakailist.append(copy);
		return copy;
	}
	else {
		opserr << "FATAL: MatsuokaNakaiDMM::getCopy() - returned NULLPTR! couldn't copy the material\n";
		exit(-1);
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
	s << "Matsuoka-Nakai Multi-Yield Surface Hardening-Softening nD Material " << getTag() << "::" << getSubTag() << "\n";
}

int MatsuokaNakaiDMM::setParameter(const char** argv, int argc, Parameter& param) {

	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	// check for material tag
	if (theMaterialTag == this->getTag()) {
		if (strcmp(argv[0], "updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
	}

	// also go through the base-class method
	return MultiYieldSurfaceHardeningSoftening::setParameter(argv, argc, param);
}

int MatsuokaNakaiDMM::updateParameter(int responseID, Information& info) {

	if (responseID == 1) {
		// update material stage and if nonlinear set up yield surfaces 
		if (info.theInt > 0) {
			// set-up yield surfaces
			for each (auto ptr in MatsuokaNakailist.mptr) {
				if (ptr->getSubTag() == 0) { continue; }	// Skip the first material since its not used (the master copy)
				if (ptr->theData == nullptr) {
					opserr << "FATAL: MatsuokaNakaiDMM::updateParameter() - nD Material " << ptr->getTag() << "."
						<< ptr->getSubTag() << " cannot access the yield surface library!";
					exit(-1);
				}
				else {
					// if the desired type of surface is available and the material stage has not been already updated
					if ((ptr->theData->isAOK(ptr->getDataDriver())) && (ptr->materialStage != info.theInt)) {
						if (beVerbose) {
							opserr << "WARNING: MatsuokaNakaiDMM::updateParameter() - nD Material " << ptr->getTag() <<
								"." << ptr->getSubTag() << " -> material stage has updated to: " << info.theInt << "\n";
						}
						ptr->updateModuli(ptr->sv.sig);
						ptr->ys = YieldSurfacePackage(ptr->getTag(), ptr->getSubTag(), ptr->getOrder(), ptr->getDataDriver(),
							ptr->theData, ptr->getGmod(), ptr->getPref(), ptr->sv.sig, ptr->sv.eps, ptr->beVerbose);
						ptr->materialStage = info.theInt;
					}
					else {
						opserr << "WARNING: MatsuokaNakaiDMM::updateParameter() - nD Material " << ptr->getTag() << "." << ptr->getSubTag() <<
							" cannot update material stage! Keeping the current stage: " << ptr->materialStage << " ...\n";
					}
				}
			}
		}
	}

	// also go through the base-class method
	return MultiYieldSurfaceHardeningSoftening::updateParameter(responseID, info);
}

// Private methods
	// the get methods
Vector MatsuokaNakaiDMM::getShiftedDeviator(const Vector& stress, const int num_ys) {
	Vector zeta = Vector(nOrd);
	Vector deviator = getStressDeviator(stress);
	Vector backpressure = ys.getAlpha(num_ys);
	for (int i = 0; i < nOrd; i++) {
		zeta(i) = deviator(i) - backpressure(i);
	}
	return zeta;
}

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

double MatsuokaNakaiDMM::linearizedFlow(const double dlambda) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::linearizedFlow() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

double MatsuokaNakaiDMM::yieldf(const Vector& stress, const int num_ys, bool yield_stress) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::yieldf() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MatsuokaNakaiDMM::Qi(const Vector& stress, const int num_ys) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::Qi() -> subclass responsibility\n";
	exit(-1);
	return 0;
}

Vector MatsuokaNakaiDMM::Di(const Vector& stress, const int num_ys) {
	opserr << "FATAL: MultiYieldSurfaceHardeningSoftening::Di() -> subclass responsibility\n";
	exit(-1);
	return 0;
}