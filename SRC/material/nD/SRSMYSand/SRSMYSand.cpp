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

// $Revision: 0.0 $
// $Date: 2023-11-01 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/SRSMYSand/SRSMYSand.cpp $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata
//				Enrico Spacone
//				Carlo G. Lai
//              Claudio Tamagnini
//
// Created in:	November 2023
//
// Description: This file contains the implementation for the SRSMYSand class.

#include "SRSMYSand.h"
#include <MaterialResponse.h>

using tc = CTensor::Constants;

namespace tools {
	// keep track of all material instances local to the current process
	class MaterialList {
	public:
		size_t size = 0;
		int next_tag = 0;
		std::vector<std::unique_ptr<SRSMYSand>> mptr;

	public:
		MaterialList(void) = default;
		~MaterialList(void) = default;

		void append(std::unique_ptr<SRSMYSand> matptr) {
			matptr->setSubTag(next_tag);
			next_tag++;
			mptr.push_back(std::move(matptr));
			size = mptr.size();
		}

		void remove(std::unique_ptr<SRSMYSand> matptr) {
			auto it = std::find(mptr.begin(), mptr.end(), matptr);
			if (it != mptr.end()) {
				mptr.erase(it);
			}
			// clean up, if the list is empty or, if not, update size
			if (mptr.empty()) { cleanup(); }
			else { size = mptr.size(); }
		}

		void cleanup(void) {
			size = 0;
			next_tag = 0;
			mptr.clear();
		}
	};

	// allocate static global storage for all instances
	class GlobalStorage {
	public:

	public:

	};
}

// Public methods
	// full constructor
SRSMYSand::SRSMYSand(int tag, int classTag,
	double r0, double K0, double G0, double P0, double m0,
	double ddtype, int itype, bool verbosity)
	:NDMaterial(tag, classTag)
{
	
}

// null constructor
SRSMYSand::SRSMYSand(void)
	: NDMaterial()
{

}

// destructor
SRSMYSand::~SRSMYSand(void)
{

}

// iteration control
int SRSMYSand::commitState(void) {

	// done
	return 0;
}

int SRSMYSand::revertToLastCommit(void) {


	// done
	return 0;
}

int SRSMYSand::revertToStart(void) {

	// done
	return 0;
}

// set material strain
int SRSMYSand::setTrialStrain(const Vector& strain) {

	// done
	return 0;
}

int SRSMYSand::setTrialStrain(const Vector& strain, const Vector& rate) {
	return setTrialStrain(strain);
}

int SRSMYSand::setTrialStrainIncr(const Vector& strain)
{
	return setTrialStrain(strain);
}

int SRSMYSand::setTrialStrainIncr(const Vector& strain, const Vector& rate) {
	return setTrialStrain(strain);
}


// return material info
NDMaterial* SRSMYSand::getCopy(void) {

	return 0;
}

NDMaterial* SRSMYSand::getCopy(const char* type) {

	return 0;
}

const char* SRSMYSand::getType(void) const {


	return 0;

}

int SRSMYSand::getOrder(void) const {


	return 0;
}

// return stress & strain
const Vector& SRSMYSand::getStress(void) {

	// done
	return Vector(6);
}

const Vector& SRSMYSand::getStrain(void) {

	// done
	return Vector(6);
}

// return the tangent
const Matrix& SRSMYSand::getTangent(void) {

	return Matrix(6 ,6);
}

const Matrix& SRSMYSand::getInitialTangent(void) {

	return Matrix(6, 6);
}

// return material response
Response* SRSMYSand::setResponse(const char** argv, int argc, OPS_Stream& output)
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
	else
		return 0;
}

int SRSMYSand::getResponse(int responseID, Information& matInfo)
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
	default:
		return -1;
	}
}

// parallel message passing
int SRSMYSand::sendSelf(int commitTag, Channel& theChannel) {

	return 0;
}

int SRSMYSand::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {

	return 0;
}

// miscellaneous
void SRSMYSand::Print(OPS_Stream& s, int flag) {

}

int SRSMYSand::setParameter(const char** argv, int argc, Parameter& param) {
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);


	// done
	return -1;
}

int SRSMYSand::updateParameter(int responseID, Information& info) {



	// done
	return 0;
}

