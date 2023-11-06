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
	// move operation
	order = other.order;
	dim = other.dim;
	numRows = other.numRows;
	numCols = other.numCols;
	repr = other.repr;
	ct = other.ct;
}

// second-order tensor
CTensor::CTensor(int nRows, int rep)
	:ct(nRows, 1), numRows(nRows), numCols(1), repr(rep)
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
	:ct(data, nRows, 1), numRows(nRows), numCols(1), repr(rep)
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
	:numCols(1), repr(rep)
{
	if (rep < 0 || rep > 2) {
		opserr << "FATAL! CTensor::CTensor() - second order tensor does not have a mixed representation!\n";
		exit(-1);
	}
	// initialize from Vector
	numRows = V.Size();
	ct = Matrix(numRows, 1);
	setOrder(2);
	matrix_dim(numRows);
	for (int i = 0; i < numRows; i++)
		ct(i, 0) = V(i);
}

// fourth-order tensor
CTensor::CTensor(int nRows, int nCols, int rep)
	:ct(nRows, nCols), numRows(nRows), numCols(nCols), repr(rep)
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
	:ct(data, nRows, nCols), numRows(nRows), numCols(nCols), repr(rep)
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
	numRows = M.noRows();
	numCols = M.noCols();
	setOrder(4);
	matrix_dim(numRows);
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

int CTensor::getRep(void) { return repr; }
int CTensor::getOrder(void) { return order; }
int CTensor::length(void) const { return (numRows * numCols); }
int CTensor::noRows(void) const { return numRows; }
int CTensor::noCols(void) const { return numCols; }

double CTensor::trace(void) {
	double sum = 0;
	if (order == 2) {
		sum = ct(0, 0) + ct(1, 0);
		if (dim == 3) { sum += ct(2, 0); }
	}
	else {
		sum = operator%(Constants::IIvol(dim, repr));
	}
	return sum;
}

CTensor CTensor::deviator(void) {
	CTensor dev(ct, repr);
	if (order == 2) {
		dev -= Constants::I(dim, repr) * (trace() / double(dim));
	}
	else {
		dev -= Constants::IIvol(dim, repr) * (trace() / double(dim));
	}

	return dev;
}

Vector CTensor::makeVector(void) {
	if (order != 2) {
		opserr << "WARNING! CTensor::makeVector() - cannot make vector from a 4th order tensor!\n";
		return Vector();
	}
	Vector result(numRows);
	for (int i = 0; i < numRows; i++)
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
	numRows = nRows;
	numCols = 1;
	repr = rep;
	setOrder(2);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(int nRows) {
	ct.resize(nRows, 1);
	numRows = nRows;
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(const Vector& V) {
	int nRows = V.Size();
	ct.resize(nRows, 1);
	numRows = nRows;
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
	numRows = nRows;
	numCols = nCols;
	repr = rep;
	setOrder(4);
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(int nRows, int nCols) {
	ct.resize(nRows, nCols);
	numRows = nRows;
	numCols = nCols;
	matrix_dim(nRows);
	return 0;
}

int CTensor::resize(const Matrix& M) {
	int nRows = M.noRows();
	int nCols = M.noCols();
	ct.resize(nRows, nCols);
	numRows = nRows;
	numCols = nCols;
	matrix_dim(nRows);
	return 0;
}

int CTensor::setData(const CTensor& deviatoric, const double volumetric) {
	// move the input deviatoric CTensor
	*this = deviatoric;
	// add the volumetric part
	if (order == 2) {
		addTensor(1.0, Constants::I(dim, repr) * volumetric, 1.0);
	}
	else {
		addTensor(1.0, Constants::IIvol(dim, repr) * volumetric, 1.0);
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
	opserr << "Function not implemented yet!\n"; exit(-1);
	double deteminant = 0;;

	return deteminant;
}

double CTensor::norm(void) const {
	// compute the norm of ctensor
	double result = 0.0;
	if (order == 2) {		// 2nd order tensor
		if (repr == 0) {		// Full
			for (int i = 0; i < numRows; i++) {
				result += ct(i, 0) * ct(i, 0);
			}
		}
		else if (repr == 1) {	// Cov
			for (int i = 0; i < numRows; i++) {
				result += ct(i, 0) * ct(i, 0) - (i >= dim) * 0.5 * ct(i, 0) * ct(i, 0);
			}
		}
		else if (repr == 2) {	// Contr
			for (int i = 0; i < numRows; i++) {
				result += ct(i, 0) * ct(i, 0) + (i >= dim) * ct(i, 0) * ct(i, 0);
			}
		}
	}
	else if (order == 4) {	// 4th order tensor
		if (repr == 0) {		// Full
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					result += ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 1) {	// Cov
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					result += ct(i, j) * ct(i, j) - (i >= dim) * 0.5 * ct(i, j) * ct(i, j) - (j >= dim) * 0.5 * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 2) {	// Contr
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					result += ct(i, j) * ct(i, j) + (i >= dim) * ct(i, j) * ct(i, j) + (j >= dim) * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 3) {	// CovContr
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					result += ct(i, j) * ct(i, j) - (i >= dim) * 0.5 * ct(i, j) * ct(i, j) + (j >= dim) * ct(i, j) * ct(i, j);
				}
			}
		}
		else if (repr == 4) {	// ContrCov
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
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
		for (int i = 0; i < numRows; i++) {
			result += ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 1 && other.repr == 1) {	// Cov - Cov
		for (int i = 0; i < numRows; i++) {
			result += ct(i, 0) * other.ct(i, 0) - (i >= dim) * 0.5 * ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 2 && other.repr == 2) {	// Contr - Contr
		for (int i = 0; i < numRows; i++) {
			result += ct(i, 0) * other.ct(i, 0) + (i >= dim) * ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 1 && other.repr == 2) {	// Cov - Contr
		for (int i = 0; i < numRows; i++) {
			result += ct(i, 0) * other.ct(i, 0);
		}
	}
	else if (repr == 2 && other.repr == 1) {	// Contr - Cov
		for (int i = 0; i < numRows; i++) {
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
	int r = (dim == 2) * 3 + (dim == 3) * 6;
	int c = 1;
	double* theData = nullptr;
	CTensor result;
	if (order == 2 && other.order == 4) {		// double dot between 2nd and 4th order tensors
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
			r = (dim == 2) * 4 + (dim == 3) * 9;
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
		theData = new (nothrow) double[r];
		if (theData == nullptr) {
			opserr << "FATAL! CTensor::operator^() - memory allocation failed!\n";
			exit(-1);
		}
		for (int i = 0; i < r; i++)
			theData[i] = this->ct(i, 0);
		Vector temp(theData, r);
		delete[] theData;
		theData = nullptr;
		temp = other.ct ^ temp;
		result = CTensor(temp, rep);

	}
	else if (order == 4 && other.order == 2) {	// double dot between 4th and 2nd order tensors
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
			r = (dim == 2) * 4 + (dim == 3) * 9;
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
		theData = new (nothrow) double[r];
		if (theData == nullptr) {
			opserr << "FATAL! CTensor::operator^() - memory allocation failed!\n";
			exit(-1);
		}
		for (int i = 0; i < r; i++)
			theData[i] = other.ct(i, 0);
		Vector temp(theData, r);
		delete[] theData;
		theData = nullptr;
		temp = this->ct * temp;
		result = CTensor(temp, rep);
	}
	else if (order == 4 && other.order == 2) {	// double dot between two 4th order tensors
		if (repr == 0 && other.repr == 0) {			// Full(0) x Full(0) = Full(0)
			rep = 0;
			r = (dim == 2) * 4 + (dim == 3) * 9;
			c = r;
		}
		else if (repr == 1 && other.repr == 2) {	// Cov(1) x Contr(2) = CovContr(3)
			rep = 3;
			c = r;
		}
		else if (repr == 2 && other.repr == 1) {	// Contr(2) x Cov(1) = ContrCov(4)
			rep = 4;
			c = r;
		}
		else if (repr == 4 && other.repr == 2) {	// ContrCov(4) x Contr(2) = Contr(2)
			rep = 2;
			c = r;
		}
		else if (repr == 3 && other.repr == 1) {	// CovContr(3) x Cov(1) = Cov(1)
			rep = 1;
			c = r;
		}
		else if (repr == 4 && other.repr == 4) {	// ContrCov(4) x ContrCov(4) = ContrCov(4)
			rep = 4;
			c = r;
		}
		else if (repr == 3 && other.repr == 3) {	// CovContr(3) x CovContr(3) = CovContr(3)
			rep = 3;
			c = r;
		}
		else {
			opserr << "FATAL! CTensor::operator^() - matrices with unsupported combination of representations. cannot double dot!\n";
			exit(-1);
		}
		// compute the double dot operation
		result = CTensor(r, c, rep);
		result.ct = this->ct * other.ct;
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
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < other.numRows; j++)
			result.ct(i, j) = ct(i, 0) * other.ct(j, 0);
	}
	return result;
}

CTensor& CTensor::dot(const CTensor& C) {
	// single dot operation between two CTensors
	opserr << "Function not implemented yet!\n"; exit(-1);
	CTensor result;

	return result;
}

CTensor& CTensor::invert(void) {
	// return the inverse of ctensor
	opserr << "Function not implemented yet!\n"; exit(-1);
	CTensor result;

	return result;
}

// special norms
double CTensor::J2(void) {
	return 0.5 * deviator().norm();
}

double CTensor::octahedral(void) {
	return sqrt(1.0 / 3.0) * deviator().norm();
}

int CTensor::addTensor(double factThis, const CTensor& other, double factOther) {
	if (numRows != other.numRows && numCols != other.numCols) {
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
	if (numRows != other.numRows && numCols != other.numCols) {
		opserr << "WARNING! CTensor::addTensorTranspose() - ctensor dimension mismatch!\n";
		return -1;;
	}
	return ct.addMatrixTranspose(factThis, other.ct, factOther);
}

// overloaded operators
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
	this->order = other.order;
	this->dim = other.dim;
	this->numRows = other.numRows;
	this->numCols = other.numCols;
	this->repr = other.repr;
	this->ct = other.ct;

	return *this;
}

	// ctensor-scalar operations
CTensor CTensor::operator+(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result += fact;
	return result;
}

CTensor CTensor::operator-(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result -= fact;
	return result;
}

CTensor CTensor::operator*(double fact) const {
	// check if quick return
	if (fact == 0.0)
		return *this;

	CTensor result(*this);
	result *= fact;
	return result;
}

CTensor CTensor::operator/(double fact) const {
	if (fact == 0.0) {
		opserr << "FATAL! CTensor::operator/() - error divide-by-zero\n";
		exit(0);
	}
	CTensor result(*this);
	result /= fact;
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
	if ((numCols != other.numCols) || (numRows != other.numRows)) {
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
	if ((numCols != other.numCols) || (numRows != other.numRows)) {
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
	if ((numCols != other.numCols) || (numRows != other.numRows)) {
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
	if ((numCols != other.numCols) || (numRows != other.numRows)) {
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
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (repr == 0) {

	}
	else if (repr == 1) {
		return 0;
	}
	else if (repr ==2) {

	}
	else if (repr == 3) {

	}
	else if (repr == 4) {

	}
	else {
		opserr << "WARNING! CTensor::toCov() - cannot convert ctensor to Cov representation\n";
		return -1;
	}
}

int CTensor::toFull(void) {
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (repr == 0) {
		return 0;
	}
	else if (repr == 1) {

	}
	else if (repr == 2) {

	}
	else if (repr == 3) {

	}
	else if (repr == 4) {

	}
	else {
		opserr << "WARNING! CTensor::toFull() - cannot convert ctensor to Full representation\n";
		return -1;
	}
}

int CTensor::toContr(void) {
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (repr == 0) {

	}
	else if (repr == 1) {

	}
	else if (repr == 2) {
		return 0;
	}
	else if (repr == 3) {

	}
	else if (repr == 4) {

	}
	else {
		opserr << "WARNING! CTensor::toContr() - cannot convert ctensor to Contr representation\n";
		return -1;
	}
}

int CTensor::toCovContr(void) {
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (order != 4) {
		opserr << "WARNING! CTensor::toCovContr() - cannot convert 2nd order ctensor to mixed representation\n";
		return -1;
	}
	if (repr == 0) {

	}
	else if (repr == 1) {

	}
	else if (repr == 2) {

	}
	else if (repr == 3) {
		return 0;
	}
	else if (repr == 4) {

	}
	else {
		opserr << "WARNING! CTensor::toCovContr() - cannot convert ctensor to CovContr representation\n";
		return -1;
	}
}

int CTensor::toContrCov(void) {
	opserr << "Function not implemented yet!\n"; exit(-1);
	if (order != 4) {
		opserr << "WARNING! CTensor::toContrCov() - cannot convert 2nd order ctensor to mixed representation\n";
		return -1;
	}
	if (repr == 0) {

	}
	else if (repr == 1) {

	}
	else if (repr == 2) {

	}
	else if (repr == 3) {

	}
	else if (repr == 4) {
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
	CTensor T(m, m, rep);
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