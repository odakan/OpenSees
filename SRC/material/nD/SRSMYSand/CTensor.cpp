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

// $Revision: 1.0 $
// $Date: 2023/11/01 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/SRSMYSand/CTensor.cpp,v $

// Written: Onur Deniz Akan
// Created: 11/23
// Based on fmk's Matrix class and UWMaterials' symmetric tensor operations suite
//
// Description: This file contains the class implementation for CTensor (Compressed Tensor).
//

#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "ID.h"
#include "CTensor.h"

using std::nothrow;

constexpr double SMALL_VALUE = 1e-8;

// constructors
CTensor::CTensor(void)
{
	// initialized empty
}

CTensor::CTensor(const CTensor& other)
{
	*this = other;
}

CTensor::CTensor(const CTensor& deviatoric, const double volumetric, const bool linearized)
{
	setData(deviatoric, volumetric, linearized);
}

CTensor::CTensor(const CTensor& deviatoric, const CTensor& volumetric, const bool linearized)
{
	setData(deviatoric, volumetric, linearized);
}

// second-order tensor
CTensor::CTensor(int nRows, int rep)
	:ct(nRows, 1), repr(rep)
{
	if (rep < 0 || rep > 2) {
		opserr << "FATAL! CTensor::CTensor() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	setOrder(2);
	matrix_dim(nRows);
	// initialized with zeros
}

CTensor::CTensor(double* data, int nRows, int rep)
	:ct(data, nRows, 1), repr(rep)
{
	if (rep < 0 || rep > 2) {
		opserr << "FATAL! CTensor::CTensor() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	setOrder(2);
	matrix_dim(nRows);
	// initialized with data
}

CTensor::CTensor(const Vector& V, int rep)
	: repr(rep)
{
	if (rep < 0 || rep > 2) {
		opserr << "FATAL! CTensor::CTensor() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	// initialize from Vector
	ct = Matrix(V.Size(), 1);
	setOrder(2);
	matrix_dim(ct.noRows());
	for (int i = 0; i < ct.noRows(); i++)
		ct(i, 0) = V(i);
}

// fourth-order tensor
CTensor::CTensor(int nRows, int nCols, int rep)
	:ct(nRows, nCols), repr(rep)
{
	if (rep < 0 || rep > 4) {
		opserr << "FATAL! CTensor::CTensor() - unsupported matrix representation!\n";
		exit(-1);
	}
	setOrder(4);
	matrix_dim(nRows);
	// initialized with zeros
}

CTensor::CTensor(double* data, int nRows, int nCols, int rep)
	:ct(data, nRows, nCols), repr(rep)
{
	if (rep < 0 || rep > 4) {
		opserr << "FATAL! CTensor::CTensor() - unsupported matrix representation!\n";
		exit(-1);
	}
	setOrder(4);
	matrix_dim(nRows);
	// initialized with data
}

CTensor::CTensor(const Matrix& M, int rep)
	:ct(M), repr(rep)
{
	if (rep < 0 || rep > 4) {
		opserr << "FATAL! CTensor::CTensor() - unsupported matrix representation!\n";
		exit(-1);
	}
	setOrder(4);
	matrix_dim(ct.noRows());
	// initialized from Matrix
}

// destructor
CTensor::~CTensor()
{
	// no pointers tp clean
}

// utility methods
	// general
void CTensor::Zero(void) { ct.Zero(); }

int CTensor::setOrder(int ord) { 
	if (ord == 2 || ord == 4) {
		order = ord;
	}
	else {
		opserr << "FATAL! CTensor::setRepresentation() - unsupported matrix representation!\n";
		exit(-1);
	}
	return 0;
}

int CTensor::makeRep(int rep) {
	if (order == 2 && rep > 2) {
		opserr << "FATAL! CTensor::makeRep() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	if (rep < 0 || rep > 4) {
		opserr << "FATAL! CTensor::makeRep() - unsupported matrix representation!\n";
		exit(-1);
	}
	if (repr == rep) {
		// check for a quick return
		return 0;
	}
	else if (repr == -1) {
		// if unset, simply set
		repr = rep;
	}
	else {
		// if already set, switch to new representation
		if (rep == 0) {
			toFull();
		}
		else if (rep == 1) {
			toCov();
		}
		else if (rep == 2) {
			toContr();
		}
		else if (rep == 3) {
			toCovContr();
		}
		else {
			toContrCov();
		}
	}
	return 0;
}

int CTensor::getRep(void) const { return repr; }
int CTensor::getOrder(void) const { return order; }
int CTensor::length(void) const { return (ct.noRows() * ct.noCols()); }
int CTensor::noRows(void) const { return ct.noRows(); }
int CTensor::noCols(void) const { return ct.noCols(); }

double CTensor::trace(void) const {
	double sum = 0;
	if (order == 2) {
		sum = ct(0, 0) + ct(1, 0);
		if (dim == 3) { sum += ct(2, 0); }
	}
	else {
		opserr << "FATAL! CTensor::trace() - 4th order trace operator is not implemented yet!\n";
		exit(-1);
	}
	return sum;
}

CTensor CTensor::volumetric(const bool linearized) {
	CTensor vol(*this);
	int rep = 0;
	if (linearized) {
		if (order == 2) {
			if (repr == 1) {
				rep = 3; // CovContr
			}
			else if (repr == 2) {
				rep = 3; // ContrCov
			}
			else { rep = 0; }
			// vol = Constants::IIvol(dim, rep) ^ *this;
			addDoubleDotProduct(1.0, Constants::IIvol(dim, rep), 1.0, true);
		}
		else {
			opserr << "FATAL! CTensor::volumetric() - 4th order volumetric return is not implemented yet!\n";
			exit(-1);
		}
	}
	else {
		// multiplicative decomposition

	}

	return vol;
}

CTensor CTensor::deviator(const bool linearized) {
	CTensor dev(*this);
	int rep = 0;
	if (linearized) {
		if (order == 2) {
			if (repr == 1) {
				rep = 3; // CovContr
			}
			else if (repr == 2) {
				rep = 3; // ContrCov
			}
			else { rep = 0; }
			// dev = Constants::IIdev(dim, rep) ^ *this;
			addDoubleDotProduct(1.0, Constants::IIdev(dim, rep), 1.0, true);
		}
		else {
			opserr << "FATAL! CTensor::volumetric() - 4th order deviator return is not implemented yet!\n";
			exit(-1);
		}
	}
	else {
		// multiplicative decomposition

	}

	return dev;
}

Vector CTensor::makeVector(void) {
	if (order != 2) {
		opserr << "WARNING! CTensor::makeVector() - cannot make vector from a 4th order tensor!\n";
		return Vector();
	}
	Vector result(ct.noRows());
	for (int i = 0; i < ct.noRows(); i++)
		result(i) = this->ct(i, 0);
	return result;
}

Matrix CTensor::makeMatrix(void) {
	if (order != 4) {
		opserr << "WARNING! CTensor::makeMatrix() - cannot make matrix from a 2nd order tensor!\n";
		return Matrix();
	}
	Matrix result(this->ct);
	return result; }

// second-order tensor
int CTensor::setData(double* newData, int nRows, int rep) {
	if (rep < 0 || rep > 2) {
		opserr << "FATAL! CTensor::setData() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	ct = Matrix(newData, nRows, 1);
	repr = rep;
	setOrder(2);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(int nRows) {
	ct.resize(nRows, 1);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(const Vector& V) {
	int nRows = V.Size();
	ct.resize(nRows, 1);
	matrix_dim(nRows);
	return 0;
}

// fourth-order tensor
int CTensor::setData(double* newData, int nRows, int nCols, int rep) {
	if (rep < 0 || rep > 4) {
		opserr << "FATAL! CTensor::setData() - unsupported matrix representation!\n";
		exit(-1);
	}
	ct = Matrix(newData, nRows, nCols);
	repr = rep;
	setOrder(4);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(int nRows, int nCols) {
	ct.resize(nRows, nCols);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(const Matrix& M) {
	int nRows = M.noRows();
	int nCols = M.noCols();
	ct.resize(nRows, nCols);
	matrix_dim(nRows);
	return 0;
}

int CTensor::setData(const CTensor& deviatoric, const double volumetric, const bool linearized) {
	// move the input deviatoric CTensor
	*this = deviatoric;
	if (linearized) {
		// linearly add the volumetric part (T = Tvol + Tdev)
		if (order == 2) {
			addTensor(1.0, Constants::I(dim, repr) * volumetric, 1.0);
		}
		else {
			addTensor(1.0, Constants::IIvol(dim, repr) * volumetric, 1.0);
		}
	}
	else {
		// multiplicative composition (T = TvolTdev)
		if (order == 2) {
			*this =  dot(Constants::I(dim, repr) * volumetric, true);
		}
		else {
			*this = Constants::I(dim, repr) * volumetric ^ *this;
		}
	}

	if (ct.noCols() > 1 && order == 2) {
		opserr << "WARNING! CTensor::CTensor() - deviatoric had the wrong order specified! ctensor order is corrected to 4...\n";
		order = 4;
	}
	if (ct.noCols() == 1 && order == 4) {
		opserr << "WARNING! CTensor::CTensor() - deviatoric had the wrong order specified! ctensor order is corrected to 2...\n";
		order = 2;
	}

	return 0;
}

int CTensor::setData(const CTensor& deviatoric, const CTensor& volumetric, const bool linearized) {
	// move the input deviatoric CTensor
	*this = deviatoric;
	if (linearized) {
		// linearly add the volumetric part (T = Tvol + Tdev)
		addTensor(1.0, volumetric, 1.0);
	}
	else {
		// multiplicative composition (T = TvolTdev)
		if (order == 2) {
			addDotProduct(1.0, volumetric, 1.0, true);
		}
		else {
			addDoubleDotProduct(1.0, volumetric, 1.0, true);
		}
	}

	if (ct.noCols() > 1 && order == 2) {
		opserr << "WARNING! CTensor::CTensor() - deviatoric had the wrong order specified! ctensor order is corrected to 4...\n";
		order = 4;
	}
	if (ct.noCols() == 1 && order == 4) {
		opserr << "WARNING! CTensor::CTensor() - deviatoric had the wrong order specified! ctensor order is corrected to 2...\n";
		order = 2;
	}

	return 0;
}

// Symmetric Tensor Operations
int CTensor::Normalize(void) {
	// normalize self
	double self_norm = norm();
	if (self_norm > SMALL_VALUE) {
		ct /= self_norm;
		return 0;
	}
	return -1;
}

double CTensor::det(void) const {
	// compute the determinant of ctensor
	opserr << "FATAL! CTensor::det() - function not implemented yet!\n"; exit(-1);
	double deteminant = 0;;

	return deteminant;
}

double CTensor::norm(void) const {
	// compute the norm of ctensor
	double result = 0.0;
	if (order == 2) {		// 2nd order tensor
		if (repr == 0) {		// Full
			for (int i = 0; i < ct.noRows(); i++) {
				result += ct(i, 0) * ct(i, 0);
			}
		}
		else if (repr == 1) {	// Cov
			for (int i = 0; i < ct.noRows(); i++) {
				result += ct(i, 0) * ct(i, 0) - (i >= dim) * 0.5 * ct(i, 0) * ct(i, 0);
			}
		}
		else if (repr == 2) {	// Contr
			for (int i = 0; i < ct.noRows(); i++) {
				result += ct(i, 0) * ct(i, 0) + (i >= dim) * ct(i, 0) * ct(i, 0);
			}
		}
	}
	else if (order == 4) {	// 4th order tensor
		if (repr == 0) {		// Full
			for (int i = 0; i < ct.noRows(); i++) {
				for (int j = 0; j < ct.noCols(); j++) {
					result += ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 1) {	// Cov
			for (int i = 0; i < ct.noRows(); i++) {
				for (int j = 0; j < ct.noCols(); j++) {
					result += ct(i, j) * ct(i, j) - (i >= dim) * 0.5 * ct(i, j) * ct(i, j) - (j >= dim) * 0.5 * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 2) {	// Contr
			for (int i = 0; i < ct.noRows(); i++) {
				for (int j = 0; j < ct.noCols(); j++) {
					result += ct(i, j) * ct(i, j) + (i >= dim) * ct(i, j) * ct(i, j) + (j >= dim) * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 3) {	// CovContr
			for (int i = 0; i < ct.noRows(); i++) {
				for (int j = 0; j < ct.noCols(); j++) {
					result += ct(i, j) * ct(i, j) - (i >= dim) * 0.5 * ct(i, j) * ct(i, j) + (j >= dim) * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 4) {	// ContrCov
			for (int i = 0; i < ct.noRows(); i++) {
				for (int j = 0; j < ct.noCols(); j++) {
					result += ct(i, j) * ct(i, j) + (i >= dim) * ct(i, j) * ct(i, j) - (j >= dim) * 0.5 * ct(i, j) * ct(i, j);
				}
			}
		}
	}
	return sqrt(result);
}

double CTensor::operator%(const CTensor& other) const {
	// double dot operation between two 2nd order CTensors
	if (dim != other.dim) {
		// make trouble
		opserr << "FATAL! CTensor::operator%() - ctensor dimensions do not match!\n";
		exit(-1);
	}
	if (order != 2 || other.order !=2) {
		// make trouble
		opserr << "FATAL! CTensor::operator%() - both ctensors must be 2nd order! Use operator^() for double dot operations involving 4th order tensors...\n";
		exit(-1);
	}
	// compute double dot
	double result = 0.0;
	if (repr == 0 && other.repr == 0) {		// Full - Full
		for (int i = 0; i < ct.noRows(); i++) {
			result += ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 1 && other.repr == 1) {	// Cov - Cov
		for (int i = 0; i < ct.noRows(); i++) {
			result += ct(i, 0) * other.ct(i, 0) - (i >= dim) * 0.5 * ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 2 && other.repr == 2) {	// Contr - Contr
		for (int i = 0; i < ct.noRows(); i++) {
			result += ct(i, 0) * other.ct(i, 0) + (i >= dim) * ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 1 && other.repr == 2) {	// Cov - Contr
		for (int i = 0; i < ct.noRows(); i++) {
			result += ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 2 && other.repr == 1) {	// Contr - Cov
		for (int i = 0; i < ct.noRows(); i++) {
			result += ct(i, 0) * other.ct(i, 0);
		}
	}
	else {
		// make trouble
		opserr << "FATAL! CTensor::operator%() - double dot operation between the requested combination of matrix representations is not supported!\n";
		exit(-1);
	}
	return result;
}

CTensor CTensor::operator^(const CTensor& other) const {
	// double dot operation between two [2-4] or [4-4] CTensors
	if (dim != other.dim) { // dimensions must match [2D&2D or 3D&3D]
		// make trouble
		opserr << "FATAL! CTensor::operator^() - ctensor dimensions do not match!\n";
		exit(-1);
	}
	// decide output properties
	int rep = -1;
	double* theData = nullptr;
	CTensor result;
	if (order == 2 && other.order == 4) {		// double dot between 2nd and 4th order tensors
		if (ct.noRows() != other.ct.noRows()) {
			opserr << "FATAL! CTensor::operator^() - 2nd-4th order ctensor size do not match! [ct.noRows() != other.ct.noRows()]\n";
			exit(-1);
		}
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
		}
		else if (repr == 1 && other.repr == 2) {	// Cov(1) x Contr(2) = Contr(2)
			rep = 2;
		}
		else if (repr == 2 && other.repr == 1) {	// Contr(2) x Cov(1) = Cov(1) 
			rep = 1;
		}
		else if (repr == 1 && other.repr == 3) {	// Cov(1) x ContrCov(3) = Cov(1)
			rep = 1;
		}
		else if (repr == 2 && other.repr == 4) {	// Contr(2) x CovContr(4) = Contr(2)
			rep = 2;
		}
		else {
			opserr << "FATAL! CTensor::operator^() - matrices with unsupported combination of representations. cannot double dot!\n";
			exit(-1);
		}
		// compute the double dot operation
		result.order = 2;
		result.dim = dim;
		result.repr = rep;
		result.resize(other.ct.noCols(), 1);
		result.Zero();
		for (int i = 0; i < result.ct.noRows(); i++) {
			for (int j = 0; j < ct.noRows(); j++) {
				result.ct(i, 0) += other.ct(j, 0) * ct(j, i);
			}
		}
	}
	else if (order == 4 && other.order == 2) {	// double dot between 4th and 2nd order tensors
		if (ct.noCols() != other.ct.noRows()) {
			opserr << "FATAL! CTensor::operator^() - 4th-2nd order ctensor size do not match! [ct.noCols() != other.ct.noRows()]\n";
			exit(-1);
		}
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
		}
		else if (repr == 1 && other.repr == 2) {	// Cov(1) x Contr(2) = Cov(1)
			rep = 1;
		}
		else if (repr == 2 && other.repr == 1) {	// Contr(2) x Cov(1) = Contr(2)
			rep = 2;
		}
		else if (repr == 4 && other.repr == 2) {	// ContrCov(4) x Contr(2) = Contr(2)
			rep = 2;
		}
		else if (repr == 3 && other.repr == 1) {	// CovContr(3) x Cov(1) = Cov(1)
			rep = 1;
		}
		else {
			opserr << "FATAL! CTensor::operator^() - matrices with unsupported combination of representations. cannot double dot!\n";
			exit(-1);
		}
		// compute the double dot operation
		result.order = 2;
		result.dim = dim;
		result.repr = rep;
		result.resize(ct.noRows(), 1);
		result.Zero();
		for (int i = 0; i < result.ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				result.ct(i, 0) += ct(i, j) * other.ct(j, 0);
			}
		}
	}
	else if (order == 4 && other.order == 4) {	// double dot between two 4th order tensors
		if (ct.noRows() != other.ct.noCols()) {
			opserr << "FATAL! CTensor::operator^() - 4th order ctensor size do not match! [ct.noRows() != other.ct.noCols()]\n";
			exit(-1);
		}
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
		}
		else if (repr == 1 && other.repr == 2) {	// Cov(1) x Contr(2) = CovContr(3)
			rep = 3;
		}
		else if (repr == 2 && other.repr == 1) {	// Contr(2) x Cov(1) = ContrCov(4)
			rep = 4;
		}
		else if (repr == 4 && other.repr == 2) {	// ContrCov(4) x Contr(2) = Contr(2)
			rep = 2;
		}
		else if (repr == 3 && other.repr == 1) {	// CovContr(3) x Cov(1) = Cov(1)
			rep = 1;
		}
		else if (repr == 4 && other.repr == 4) {	// ContrCov(4) x ContrCov(4) = ContrCov(4)
			rep = 4;
		}
		else if (repr == 3 && other.repr == 3) {	// CovContr(3) x CovContr(3) = CovContr(3)
			rep = 3;
		}
		else {
			opserr << "FATAL! CTensor::operator^() - matrices with unsupported combination of representations. cannot double dot!\n";
			exit(-1);
		}
		// compute the double dot operation
		result.order = 4;
		result.dim = dim;
		result.repr = rep;
		result.resize(ct.noRows(), other.ct.noCols());
		result.Zero();
		for (int i = 0; i < result.ct.noRows(); i++) {
			for (int j = 0; j < result.ct.noCols(); j++) {
				for (int k = 0; k < result.ct.noRows(); k++) {
					result.ct(i, j) += ct(i, k) * other.ct(k, j);
				}
			}
		}
	}
	else {
		// make trouble
		opserr << "FATAL! CTensor::operator^() - matrices with unsupported combination of orders. cannot double dot!\n";
		exit(-1);
	}
	return result;
}

CTensor CTensor::operator*(const CTensor& other) const {
	// dyadic product operation between two 2nd order CTensors
	if (dim != other.dim) {
		// make trouble
		opserr << "FATAL! CTensor::operator*() - ctensor dimensions do not match!\n";
		exit(-1);
	}
	if (order != 2 || other.order != 2) {
		// make trouble
		opserr << "FATAL! CTensor::operator*() - both ctensors must be 2nd order. cannot commpute the dyadic product!\n";
		exit(-1);
	}
	// decide output matrix size
	int m = 0;
	if (repr == 0 && other.repr == 0) {			// Full - Full
		m = int((dim == 3) * (9) + (dim == 2) * (4));
	}
	else {
		m = int((dim == 3) * (6) + (dim == 2) * (3));
	}
	// decide output repsentation
	int rep = 0;
	if (repr == 0 && other.repr == 0) {		// Full - Full
		rep = 0;
	}
	else if (repr == 1 && other.repr == 1) {	// Cov - Cov
		rep = 1;
	}
	else if (repr == 2 && other.repr == 2) {	// Contr - Contr
		rep = 2;
	}
	else if (repr == 1 && other.repr == 2) {	// Cov - Contr
		rep = 3;
	}
	else if (repr == 2 && other.repr == 1) {	// Contr - Cov
		rep = 4;
	}
	else {
		// make trouble
		opserr << "FATAL! CTensor::operator*() - dyadic product between the requested combination of matrix representations is not supported!\n";
		exit(-1);
	}
	// compute dyadic product
	CTensor result(m, m, rep);
	for (int i = 0; i < ct.noRows(); i++) {
		for (int j = 0; j < other.ct.noRows(); j++)
			result.ct(i, j) += ct(i, 0) * other.ct(j, 0);
	}
	return result;
}

CTensor CTensor::dot(const CTensor& other, bool pre = false) const {
	// single dot operation between two CTensors
	opserr << "FATAL! CTensor::dot() - function not implemented yet!\n"; exit(-1);
	CTensor result;

	if (order == 2 && other.order == 2) {
		if (repr == 1 && other.repr == 1) {

		}
		else if (repr == 2 && other.repr == 2) {

		}
		else if (repr == 1 && other.repr == 2) {

		}
		else if (repr == 2 && other.repr == 1) {

		}

		if (dim == 2) {
			result(0) = ct(0, 0) * other.ct(0, 0) + ct(2, 0) * other.ct(2, 0);
			result(1) = ct(2, 0) * other.ct(2, 0) + ct(1, 0) * other.ct(1, 0);
			result(2) = 0.5 * (ct(0, 0) * other.ct(2, 0) + ct(2, 0) * other.ct(0, 0) + ct(3, 0) * other.ct(1, 0) + ct(1, 0) * other.ct(3, 0) + ct(5, 0) * other.ct(4, 0) + ct(4, 0) * other.ct(5, 0));
		}
		else {
			result(0) = ct(0, 0) * other.ct(0, 0) + ct(3, 0) * other.ct(3, 0) + ct(5, 0) * other.ct(5, 0);
			result(1) = ct(3, 0) * other.ct(3, 0) + ct(1, 0) * other.ct(1, 0) + ct(4, 0) * other.ct(4, 0);
			result(2) = ct(5, 0) * other.ct(5, 0) + ct(4, 0) * other.ct(4, 0) + ct(2, 0) * other.ct(2, 0);
			result(3) = 0.5 * (ct(0, 0) * other.ct(3, 0) + ct(3, 0) * other.ct(0, 0) + ct(3, 0) * other.ct(1, 0) + ct(1, 0) * other.ct(3, 0) + ct(5, 0) * other.ct(4, 0) + ct(4, 0) * other.ct(5, 0));
			result(4) = 0.5 * (ct(3, 0) * other.ct(5, 0) + ct(5, 0) * other.ct(3, 0) + ct(1, 0) * other.ct(4, 0) + ct(4, 0) * other.ct(1, 0) + ct(4, 0) * other.ct(2, 0) + ct(2, 0) * other.ct(4, 0));
			result(5) = 0.5 * (ct(0, 0) * other.ct(5, 0) + ct(5, 0) * other.ct(0, 0) + ct(3, 0) * other.ct(4, 0) + ct(4, 0) * other.ct(3, 0) + ct(5, 0) * other.ct(2, 0) + ct(2, 0) * other.ct(5, 0));
		}

	}
	else if (order == 2 && other.order == 4) {

	}
	else if (order == 4 && other.order == 2) {

	}
	else if (order == 4 && other.order == 4) {

	}
	else {

	}

	return result;
}

CTensor& CTensor::invert(void) {
	// return the inverse of ctensor
	opserr << "FATAL! CTensor::invert() - function not implemented yet!\n"; exit(-1);
	CTensor result;

	return result;
}

// special norms
double CTensor::J2(void) const {
	return 0.5 * (*this % *this);
}

double CTensor::octahedral(void) const {
	return sqrt(2.0 / 3.0 * J2());
}

int CTensor::addTensor(double factThis, const CTensor& other, double factOther) {
	if (ct.noRows() != other.ct.noRows() && ct.noCols() != other.ct.noCols()) {
		opserr << "WARNING! CTensor::addTensor() - ctensor dimension mismatch!\n";
		return -1;;
	}
	return ct.addMatrix(factThis, other.ct, factOther);
}

int CTensor::addTensorTranspose(double factThis, const CTensor& other, double factOther) {
	if (order != 4 && other.order != 4) {
		opserr << "WARNING! CTensor::addTensorTranspose() - ctensor order must be 4!\n";
		return -1;;
	}
	if (ct.noRows() != other.ct.noRows() && ct.noCols() != other.ct.noCols()) {
		opserr << "WARNING! CTensor::addTensorTranspose() - ctensor dimension mismatch!\n";
		return -1;;
	}
	return ct.addMatrixTranspose(factThis, other.ct, factOther);
}

int CTensor::addDotProduct(double factThis, const CTensor& other, double factOther, bool premultiply = false) {

}

int CTensor::addDoubleDotProduct(double factThis, const CTensor& other, double factOther, bool premultiply = false) {

}


int CTensor::addTensorProduct(double factThis, const CTensor& other, double factOther, bool premultiply = false) {

}

// overloaded operators
bool CTensor::operator==(const CTensor& other) const {
	if (order != other.order) return false;
	if (dim != other.dim) return false;
	if (ct.noRows() != other.ct.noRows()) return false;
	if (ct.noCols() != other.ct.noCols()) return false;
	if (repr != other.repr) return false;
	for (int i = 1; i < ct.noRows(); i++)
	{
		for (int j = 0; j < ct.noCols(); j++)
		{
			if (ct(i, j) != other.ct(i, j))
				return false;
		}
	}
	return true;
}

bool CTensor::operator!=(const CTensor& other) const {
	if (order != other.order) return true;
	if (dim != other.dim) return true;
	if (ct.noRows() != other.ct.noRows()) return true;
	if (ct.noCols() != other.ct.noCols()) return true;
	if (repr != other.repr) return true;
	for (int i = 1; i < ct.noRows(); i++)
	{
		for (int j = 0; j < ct.noCols(); j++)
		{
			if (ct(i, j) != other.ct(i, j))
				return true;
		}
	}
	return false;
}

double& CTensor::operator()(int row) { 
	if (order == 2) { 
		return ct(row, 0);
	} 
	else { 
		opserr << "FATAL! CTensor::operator() - this is a 4th order tensor, but (row) data was requested!\n"; 
		exit(-1); 
	} 
}

double CTensor::operator()(int row) const { 
	if (order == 2) { 
		return ct(row, 0); 
	} 
	else { 
		opserr << "FATAL! CTensor::operator() - this is a 4th order tensor, but (row) data was requested!\n"; 
		exit(-1); 
	} 
}

double& CTensor::operator()(int row, int col) { 
	if (order == 4) { 
		return ct(row, col); 
	} else { 
		opserr << "FATAL! CTensor::operator() - this is a 2nd order tensor, but (row, col) data was requested!\n"; 
		exit(-1); 
	} 
}

double CTensor::operator()(int row, int col) const { 
	if (order == 4) { 
		return ct(row, col); 
	} 
	else { 
		opserr << "FATAL! CTensor::operator() - this is a 2nd order tensor, but (row, col) data was requested!\n"; 
		exit(-1); 
	} 
}

CTensor CTensor::operator()(const ID& rows, const ID& cols) const {
	int nRows = rows.Size();
	int nCols = cols.Size();
	CTensor result(nRows, nCols, repr);
	result.ct = this->ct(rows, cols);
	return result;
}

CTensor& CTensor::operator=(const CTensor& other) {
	// first check we are not trying other = other
	if (this == &other)
		return *this;

	// assignment operation
	order = other.order;
	dim = other.dim;
	repr = other.repr;
	ct.resize(other.noRows(), other.noCols());
	ct = other.ct;
	if (ct.noCols() > 1 && order == 2) {
		opserr << "WARNING! CTensor::operator=() - other had the wrong order specified! ctensor order is corrected to 4...\n";
		order = 4;
	}
	if (ct.noCols() == 1 && order == 4) {
		opserr << "WARNING! CTensor::operator=() - other had the wrong order specified! ctensor order is corrected to 2...\n";
		order = 2;
	}
	return *this;
}

	// ctensor-scalar operations
CTensor CTensor::operator+(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result.ct += fact;
	return result;
}

CTensor CTensor::operator-(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result.ct -= fact;
	return result;
}

CTensor CTensor::operator*(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result.ct *= fact;
	return result;
}

CTensor CTensor::operator/(double fact) const {
	if (fact == 0.0) {
		opserr << "FATAL! CTensor::operator/() - error divide-by-zero\n";
		exit(0);
	}
	CTensor result(*this);
	result.ct /= fact;
	return result;
}

CTensor& CTensor::operator+=(double fact) {
	// check if quick return
	if (fact == 0.0)
		return *this;

	this->ct += fact;
	return *this;
}

CTensor& CTensor::operator-=(double fact) {
	// check if quick return
	if (fact == 0.0)
		return *this;

	this->ct += fact;
	return *this;
}

CTensor& CTensor::operator*=(double fact) {
	// check if quick return
	if (fact == 0.0)
		return *this;

	this->ct += fact;
	return *this;
}

CTensor& CTensor::operator/=(double fact) {
	if (fact == 0.0) {
		opserr << "FATAL! CTensor::operator/=() - error divide-by-zero\n";
		exit(0);
	}
	this->ct /= fact;
	return *this;
}

	// other ctensor-ctensor operations
CTensor CTensor::operator+(const CTensor& other) const {
	if (other.repr != repr) {
		opserr << "WARNING! CTensor::operator+() - ctensor matrix representations does not match!\n";
		exit(-1);
	}
	if ((ct.noCols() != other.ct.noCols()) || (ct.noRows() != other.ct.noRows())) {
		opserr << "WARNING! CTensor::operator+() - ctensor dimensions do not match!\n";
		exit(-1);
	}
	CTensor result(*this);
	result += other;
	return result;
}

CTensor CTensor::operator-(const CTensor& other) const {
	if (other.repr != repr) {
		opserr << "WARNING! CTensor::operator-() - ctensor matrix representations does not match!\n";
		exit(-1);
	}
	if ((ct.noCols() != other.ct.noCols()) || (ct.noRows() != other.ct.noRows())) {
		opserr << "WARNING! CTensor::operator-() - ctensor dimensions do not match!\n";
		exit(-1);
	}
	CTensor result(*this);
	result -= other;
	return result;
}

CTensor& CTensor::operator+=(const CTensor& other) {
	if (other.repr != repr) {
		opserr << "WARNING! CTensor::operator+=() - ctensor matrix representations does not match!\n";
		exit(-1);
	}
	if ((ct.noCols() != other.ct.noCols()) || (ct.noRows() != other.ct.noRows())) {
		opserr << "WARNING! CTensor::operator+=() - ctensor dimensions do not match!\n";
		exit(-1);
	}

	this->ct += other.ct;
	return *this;;
}

CTensor& CTensor::operator-=(const CTensor& other) {
	if (other.repr != repr) {
		opserr << "WARNING! CTensor::operator-=() - ctensor matrix representations does not match!\n";
		exit(-1);
	}
	if ((ct.noCols() != other.ct.noCols()) || (ct.noRows() != other.ct.noRows())) {
		opserr << "WARNING! CTensor::operator-=() - ctensor dimensions do not match!\n";
		exit(-1);
	}

	this->ct -= other.ct;
	return *this;;
}

OPS_Stream& operator<<(OPS_Stream& s, const CTensor& tensor) {
	char* rep = "";
	if (tensor.repr == 0) {
		rep = "-> rep: full\n";
	}
	else if (tensor.repr == 1) {
		rep = "-> rep: covariant\n";
	}
	else if (tensor.repr == 2) {
		rep = "-> rep: contravariant\n";
	}
	else if (tensor.repr == 3) {
		rep = "-> rep: mixed (covcontr)\n";
	}
	else if (tensor.repr == 4) {
		rep = "-> rep: mixed (contrcov)\n";
	}

	s << endln;
	tensor.ct.Output(s);
	s << rep;
	s << endln;
	return s;
}

int CTensor::toCov(void) {
	if (repr == 0) {		// Full to Covariant
		opserr << "FATAL! CTensor::toCov() - function not implemented yet!\n"; exit(-1);
	}
	else if (repr == 1) {	// Covariant to Covariant
		return 0;
	}
	else if (repr ==2) {	// Contravariant to Covariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 2.0;
				}
				if (j >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 1;
		return 0;
	}
	else if (repr == 3) {	// CovContr to Covariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (j >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 1;
		return 0;
	}
	else if (repr == 4) {	// ContrCov to Covariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 1;
		return 0;
	}
	else {
		opserr << "WARNING! CTensor::toCov() - cannot convert ctensor to Cov representation\n";
		return -1;
	}
}

int CTensor::toFull(void) {
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (repr == 0) {		// Full to Full
		return 0;
	}
	else if (repr == 1) {	// Covariant to Full

	}
	else if (repr == 2) {	// Contravariant to Full

	}
	else if (repr == 3) {	// CovContr to Full

	}
	else if (repr == 4) {	// ContrCov to Full

	}
	else {
		opserr << "WARNING! CTensor::toFull() - cannot convert ctensor to Full representation\n";
		return -1;
	}
}

int CTensor::toContr(void) {
	if (repr == 0) {	// Full to Contravariant
		opserr << "FATAL! CTensor::toContr() - function not implemented yet!\n"; exit(-1);
	}
	else if (repr == 1) {	// Covariant to Contravariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 0.5;
				}
				if (j >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 2;
		return 0;
	}
	else if (repr == 2) {	// Contravariant to Contravariant
		return 0;
	}
	else if (repr == 3) {	// CovContr to Contravariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 2;
		return 0;
	}
	else if (repr == 4) {	// ContrCov to Contravariant
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (j >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 2;
		return 0;
	}
	else {
		opserr << "WARNING! CTensor::toContr() - cannot convert ctensor to Contr representation\n";
		return -1;
	}
}

int CTensor::toCovContr(void) {
	if (order != 4) {
		opserr << "WARNING! CTensor::toCovContr() - cannot convert 2nd order ctensor to mixed representation\n";
		return -1;
	}
	if (repr == 0) {		// Full to CovContr
		opserr << "FATAL! CTensor::toCovContr() - function not implemented yet!\n"; exit(-1);
	}
	else if (repr == 1) {	// Covariant to CovContr
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (j >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 3;
		return 0;
	}
	else if (repr == 2) {	// Contravariant to CovContr
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 3;
		return 0;
	}
	else if (repr == 3) {	// CovContr to CovContr
		return 0;
	}
	else if (repr == 4) {	// ContrCov to CovContr
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 2.0;
				}
				if (j >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 3;
		return 0;
	}
	else {
		opserr << "WARNING! CTensor::toCovContr() - cannot convert ctensor to CovContr representation\n";
		return -1;
	}
}

int CTensor::toContrCov(void) {
	if (order != 4) {
		opserr << "WARNING! CTensor::toContrCov() - cannot convert 2nd order ctensor to mixed representation\n";
		return -1;
	}
	if (repr == 0) {		// Full to ContrCov
		opserr << "FATAL! CTensor::toContrCov() - function not implemented yet!\n"; exit(-1);
	}
	else if (repr == 1) {	// Covariant to ContrCov
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 0.5;
				}
			}
		}
		repr = 4;
		return 0;
	}
	else if (repr == 2) {	// Contravariant to ContrCov
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (j >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 4;
		return 0;
	}
	else if (repr == 3) {	// CovContr to ContrCov
		for (int i = 0; i < ct.noRows(); i++) {
			for (int j = 0; j < ct.noCols(); j++) {
				if (i >= dim) {
					ct(i, j) *= 0.5;
				}
				if (j >= dim) {
					ct(i, j) *= 2.0;
				}
			}
		}
		repr = 4;
		return 0;
	}
	else if (repr == 4) {	// ContrCov to ContrCov
		return 0;
	}
	else {
		opserr << "WARNING! CTensor::toContrCov() - cannot convert ctensor to ContrCov representation\n";
		return -1;
	}
}

void CTensor::matrix_dim(int nRows) { dim = int((nRows < 6) * 2 + (nRows >= 6) * 3); }

// tensor valued contants
const CTensor CTensor::Constants::I(const int nD, int rep = 2) {
	// return identity matrix (independent of representation)
	int m;
	if (rep == 0) { // Full
		m = int((nD == 3) * (9) + (nD == 2) * (4));
	}
	else {
		m = int((nD == 3) * (6) + (nD == 2) * (3));
	}
	CTensor T(m, rep);
	for (int i = 0; i < nD; ++i) {
		T(i) = 1.0;
	}
	return T;
}

const CTensor CTensor::Constants::IIsymm(const int nD, int rep = 2) {
	// return 4th order symmetric operator
	// 
	// e.g., (stress_deviator) = 2G[IIsymm - 1/3*p*IIvol](strain) 
	//							-> [IIsymm - 1/3*p*IIvol] = IIdev -> Contr representation
	// 
	// e.g., (stress_deviator) = 2G[IIsymm - 1/3*p*IIvol](stress) 
	//							-> [IIsymm - 1/3*p*IIvol] = IIdev -> Mixed representation

	int m;
	CTensor T;
	if (rep == 0) {				// Full
		m = int((nD == 3) * (9) + (nD == 2) * (4));
		T = CTensor(m, m, rep);
		for (int i = 0; i < m; i++)
			T(i, i) = 1.0 - (i >= nD) * 0.5;
	}
	else {
		m = int((nD == 3) * (6) + (nD == 2) * (3));
		T = CTensor(m, m, rep);
		if (rep == 1) {			// Cov
			for (int i = 0; i < m; i++)
				T(i, i) = 1.0 + (i >= nD) * 1.0;
		}
		else if (rep == 2) {	// Contr
			for (int i = 0; i < m; i++)
				T(i, i) = 1.0 - (i >= nD) * 0.5;
		}
		else if (rep > 2 && rep < 5) {	// CovContr and ContrCov (both mixed representations are the same)
			for (int i = 0; i < m; i++)
				T(i, i) = 1.0;
		}
		else {
			// make trouble
			opserr << "FATAL! CTensor::Constants::IIsymm() - unsupported matrix representation!\n";
			exit(-1);
		}
	}
	return T;
}

const CTensor CTensor::Constants::IIvol(const int nD, int rep = 2) {
	// return 4th order volumetric operator (independent of representation)
	// IIvol = I1 tensor I1
	int m;
	if (rep == 0) { // Full
		m = int((nD == 3) * (9) + (nD == 2) * (4));
	}
	else {
		m = int((nD == 3) * (6) + (nD == 2) * (3));
	}
	CTensor T(m, m, rep);
	if (nD == 3) {
		T(0, 0) = 1;
		T(0, 1) = 1;
		T(0, 2) = 1;
		T(1, 0) = 1;
		T(1, 1) = 1;
		T(1, 2) = 1;
		T(2, 0) = 1;
		T(2, 1) = 1;
		T(2, 2) = 1;
	}
	else if (nD == 2) {
		T(0, 0) = 0.5;
		T(0, 1) = 0.5;
		T(1, 0) = 0.5;
		T(1, 1) = 0.5;
	}
	return T;
}

const CTensor CTensor::Constants::IIdev(const int nD, int rep = 2) {
	//return 4th order deviatoric operator
	// 
	// e.g., (stress_deviator) = 2G[IIdev](strain) 
	//							 -> IIdev -> Contr representation
	// 
	// e.g., (stress_deviator) = 2G[IIdev](stress) 
	//							 -> IIdev -> Mixed representation

	int m;
	CTensor T;
	if (rep == 0) {				// Full
		m = int((nD == 3) * (9) + (nD == 2) * (4));
		T = CTensor(m, m, rep);
		if (nD == 3) {
			T(0, 0) = 2.0 / 3.0;
			T(0, 1) = -1.0 / 3.0;
			T(0, 2) = -1.0 / 3.0;
			T(1, 0) = -1.0 / 3.0;
			T(1, 1) = 2.0 / 3.0;
			T(1, 2) = -1.0 / 3.0;
			T(2, 0) = -1.0 / 3.0;
			T(2, 1) = -1.0 / 3.0;
			T(2, 2) = 2.0 / 3.0;
			T(3, 3) = 0.5;
			T(4, 4) = 0.5;
			T(5, 5) = 0.5;
			T(6, 6) = 0.5;
			T(7, 7) = 0.5;
			T(8, 8) = 0.5;
		}
		else if (nD == 2) {
			T(0, 0) = 1.0 / 3.0;
			T(0, 1) = -0.5 / 3.0;
			T(1, 0) = -0.5 / 3.0;
			T(1, 1) = 1.0 / 3.0;
			T(2, 2) = 0.25;
			T(3, 3) = 0.25;
		}
	}
	else {
		m = int((nD == 3) * (6) + (nD == 2) * (3));
		T = CTensor(m, m, rep);
		if (rep == 1) {			// Cov
			if (nD == 3) {
				T(0, 0) = 2.0 / 3.0;
				T(0, 1) = -1.0 / 3.0;
				T(0, 2) = -1.0 / 3.0;
				T(1, 0) = -1.0 / 3.0;
				T(1, 1) = 2.0 / 3.0;
				T(1, 2) = -1.0 / 3.0;
				T(2, 0) = -1.0 / 3.0;
				T(2, 1) = -1.0 / 3.0;
				T(2, 2) = 2.0 / 3.0;
				T(3, 3) = 2.0;
				T(4, 4) = 2.0;
				T(5, 5) = 2.0;
			}
			else if (nD == 2) {
				T(0, 0) = 1.0 / 3.0;
				T(0, 1) = -0.5 / 3.0;
				T(1, 0) = -0.5 / 3.0;
				T(1, 1) = 1.0 / 3.0;
				T(2, 2) = 1.0;
			}
		}
		else if (rep == 2) {	// Contr
			if (nD == 3) {
				T(0, 0) = 2.0 / 3.0;
				T(0, 1) = -1.0 / 3.0;
				T(0, 2) = -1.0 / 3.0;
				T(1, 0) = -1.0 / 3.0;
				T(1, 1) = 2.0 / 3.0;
				T(1, 2) = -1.0 / 3.0;
				T(2, 0) = -1.0 / 3.0;
				T(2, 1) = -1.0 / 3.0;
				T(2, 2) = 2.0 / 3.0;
				T(3, 3) = 0.5;
				T(4, 4) = 0.5;
				T(5, 5) = 0.5;
			}
			else if (nD == 2) {
				T(0, 0) = 1.0 / 3.0;
				T(0, 1) = -0.5 / 3.0;
				T(1, 0) = -0.5 / 3.0;
				T(1, 1) = 1.0 / 3.0;
				T(2, 2) = 0.25;
			}
		}
		else if (rep > 2 && rep < 5) {	// CovContr and ContrCov (both mixed representations are the same)
			if (nD == 3) {
				T(0, 0) = 2.0 / 3.0;
				T(0, 1) = -1.0 / 3.0;
				T(0, 2) = -1.0 / 3.0;
				T(1, 0) = -1.0 / 3.0;
				T(1, 1) = 2.0 / 3.0;
				T(1, 2) = -1.0 / 3.0;
				T(2, 0) = -1.0 / 3.0;
				T(2, 1) = -1.0 / 3.0;
				T(2, 2) = 2.0 / 3.0;
				T(3, 3) = 1.0;
				T(4, 4) = 1.0;
				T(5, 5) = 1.0;
			}
			else if (nD == 2) {
				T(0, 0) = 1.0 / 3.0;
				T(0, 1) = -0.5 / 3.0;
				T(1, 0) = -0.5 / 3.0;
				T(1, 1) = 1.0 / 3.0;
				T(2, 2) = 0.5;
			}
		}
		else {
			// make trouble
			opserr << "FATAL! CTensor::Constants::IIdev() - unsupported matrix representation!\n";
			exit(-1);
		}
	}
	return T;
}