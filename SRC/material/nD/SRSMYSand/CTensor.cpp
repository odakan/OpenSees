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
// Based on: fmk's Vector and Matrix classes
//
// Description: This file contains the class implementation for CTensor (Compressed Tensor).
//

#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "ID.h"
#include "CTensor.h"

using std::nothrow;

#define CTENSOR_WORK_AREA 400
#define INT_WORK_AREA 20

int CTensor::sizeDoubleWork = CTENSOR_WORK_AREA;
int CTensor::sizeIntWork = INT_WORK_AREA;
double CTensor::CTENSOR_NOT_VALID_ENTRY = 0.0;
double* CTensor::ctensorWork = 0;
int* CTensor::intWork = 0;

// constructors
CTensor::CTensor(void)
{
	// allocate work areas if the first ctensor
	if (ctensorWork == 0) {
		ctensorWork = new (nothrow) double[sizeDoubleWork];
		intWork = new (nothrow) int[sizeIntWork];
		if (ctensorWork == 0 || intWork == 0) {
			opserr << "WARNING: CTensor(void) - out of memory creating work area's\n";
			exit(-1);
		}
	}
}

CTensor::CTensor(const CTensor& other)
{
	// allocate work areas if the first ctensor
	if (ctensorWork == 0) {
		ctensorWork = new (nothrow) double[sizeDoubleWork];
		intWork = new (nothrow) int[sizeIntWork];
		if (ctensorWork == 0 || intWork == 0) {
			opserr << "WARNING: CTensor(CTensor& other) - out of memory creating work area's\n";
			exit(-1);
		}
	}

	order = other.order;
	numRows = other.numRows;
	numCols = other.numCols;
	dataSize = other.dataSize;
	representation = other.representation;
	fromFree = other.fromFree;

	if (other.order == 2) {
		if (dataSize != 0) {
			theData = new (nothrow) double[other.dataSize];

			if (theData == 0) {
				opserr << "WARNING: CTensor(CTensor& other) - ran out of memory on init of size " << dataSize << "\n";
			}
		}
		// copy the component data
		for (int i = 0; i < dataSize; i++)
			theData[i] = other.theData[i];
	}
	else if (other.order == 4) {
		if (dataSize != 0) {
			theData = new (nothrow) double[dataSize];
			if (theData == 0) {
				opserr << "WARNING: CTensor(CTensor& other) - ran out of memory on init of size " << dataSize << "\n";
				numRows = 0; numCols = 0; dataSize = 0;
			}
			else {
				// copy the data
				double* dataPtr = theData;
				double* otherDataPtr = other.theData;
				for (int i = 0; i < dataSize; i++)
					*dataPtr++ = *otherDataPtr++;
			}
		}
	}
	else {
		opserr << "\n";
		exit(-1);
	}
}

// second-order tensor
CTensor::CTensor(int size, char* rep)
	:numRows(size), dataSize(size), representation(rep), order(2)
{
	// allocate work areas if the first ctensor
	if (ctensorWork == 0) {
		ctensorWork = new (nothrow) double[sizeDoubleWork];
		intWork = new (nothrow) int[sizeIntWork];
		if (ctensorWork == 0 || intWork == 0) {
			opserr << "WARNING: CTensor(int size, char* rep) - out of memory creating work area's\n";
			exit(-1);
}
	}

#ifdef _G3DEBUG
	if (sz < 0) {
		opserr << "WARNING: CTensor(int size, char* rep) - size " << size << " specified < 0\n";
		sz = 1;
	}
#endif

	if (size > 0) {
		theData = new (nothrow) double[size];
		if (theData == 0) {
			opserr << "WARNING: CTensor(int size, char* rep) - ran out of memory on init of size " << size << "\n";
			dataSize = 0;
		}
		// zero the components
		for (int i = 0; i < dataSize; i++)
			theData[i] = 0.0;
	}
}



// fourth-order tensor
CTensor::CTensor(int nRows, int nCols, char* rep)
	:numRows(nRows), numCols(nCols), dataSize(nRows * nCols), representation(rep), order(4)
{
	// allocate work areas if the first ctensor
	if (ctensorWork == 0) {
		ctensorWork = new (nothrow) double[sizeDoubleWork];
		intWork = new (nothrow) int[sizeIntWork];
		if (ctensorWork == 0 || intWork == 0) {
			opserr << "WARNING: CTensor(int nRows, int nCols, char* rep) - out of memory creating work area's\n";
			exit(-1);
		}
	}

#ifdef _G3DEBUG
	if (nRows < 0) {
		opserr << "WARNING: CTensor(int nRows, int nCols, char* rep) - tried to init ctensor ";
		opserr << "with num rows: " << nRows << " <0\n";
		numRows = 0; numCols = 0; dataSize = 0; data = 0;
	}
	if (nCols < 0) {
		opserr << "WARNING: CTensor(int nRows, int nCols, char* rep) - tried to init ctensor";
		opserr << "with num cols: " << nCols << " <0\n";
		numRows = 0; numCols = 0; dataSize = 0; data = 0;
	}
#endif
	dataSize = numRows * numCols;
	theData = 0;

	if (dataSize > 0) {
		theData = new (nothrow) double[dataSize];
		//data = (double *)malloc(dataSize*sizeof(double));
		if (theData == 0) {
			opserr << "WARNING: CTensor(int nRows, int nCols, char* rep) - ran out of memory on init of size " << dataSize << "\n";
			numRows = 0; numCols = 0; dataSize = 0;
		}
		else {
			// zero the data
			double* dataPtr = theData;
			for (int i = 0; i < dataSize; i++)
				*dataPtr++ = 0.0;
		}
	}
}


// destructor
CTensor::~CTensor()
{
	if (theData != 0) {
		if (fromFree == 0 && dataSize > 0) {
			delete[] theData;
			theData = 0;
		}
	}
}





























// tensor valued contants
const CTensor CTensor::Constants::I(const int N) {
	// return identity matrix in Voigt notation
	CTensor kD(N, "covariant");
	int Dim = 0;
	kD.Zero();
	if (N == 3) {
		Dim = 2;
	}
	else if (N == 6) {
		Dim = 3;
	}
	else {
		opserr << "FATAL: CTensor::I() - invalid material dimension!!\n";
		exit(-1);
	}
	for (int i = 0; i < Dim; ++i) {
		kD(i) = 1.0;
	}
	return kD;
}

const CTensor CTensor::Constants::I4(const int N) {
	// return 4th order identity tennsor in Voigt notation
	CTensor I4(N, N, "covariant");
	if (N == 6) {
		I4(0, 0) = 1;
		I4(1, 1) = 1;
		I4(2, 2) = 1;
		I4(3, 3) = 0.5;
		I4(4, 4) = 0.5;
		I4(5, 5) = 0.5;
	}
	else if (N == 3) {
		I4(0, 0) = 0.5;
		I4(1, 1) = 0.5;
		I4(2, 2) = 0.5;
	}
	else {
		opserr << "FATAL: CTensor::I4() - invalid material dimension!!\n";
		exit(-1);
	}
	return I4;
}

const CTensor CTensor::Constants::IIvol(const int N) {
	// return 4th order Volumetric Tensor in Voigt notation
	// IIvol = I1 tensor I1
	CTensor IIvol(N, N, "covariant");
	if (N == 6) {
		IIvol(0, 0) = 1;
		IIvol(0, 1) = 1;
		IIvol(0, 2) = 1;
		IIvol(1, 0) = 1;
		IIvol(1, 1) = 1;
		IIvol(1, 2) = 1;
		IIvol(2, 0) = 1;
		IIvol(2, 1) = 1;
		IIvol(2, 2) = 1;
	}
	else if (N == 3) {
		IIvol(0, 0) = 0.5;
		IIvol(0, 1) = 0.5;
		IIvol(1, 0) = 0.5;
		IIvol(1, 1) = 0.5;
	}
	else {
		opserr << "FATAL: CTensor::Constants::IIvol() - invalid material dimension!!\n";
		exit(-1);
	}
	return IIvol;
}

const CTensor CTensor::Constants::IIdev(const int N) {
	//return 4th order Deviatoric Tensor in Voigt notation
	CTensor IIdev(N, N, "covariant");
	if (N == 6) {
		IIdev(0, 0) = 2.0 / 3.0;
		IIdev(0, 1) = -1.0 / 3.0;
		IIdev(0, 2) = -1.0 / 3.0;
		IIdev(1, 0) = -1.0 / 3.0;
		IIdev(1, 1) = 2.0 / 3.0;
		IIdev(1, 2) = -1.0 / 3.0;
		IIdev(2, 0) = -1.0 / 3.0;
		IIdev(2, 1) = -1.0 / 3.0;
		IIdev(2, 2) = 2.0 / 3.0;
		IIdev(3, 3) = 0.5;
		IIdev(4, 4) = 0.5;
		IIdev(5, 5) = 0.5;
	}
	else if (N == 3) {
		IIdev(0, 0) = 1.0 / 3.0;
		IIdev(0, 1) = -0.5 / 3.0;
		IIdev(1, 0) = -0.5 / 3.0;
		IIdev(1, 1) = 1.0 / 3.0;
		IIdev(2, 2) = 0.25;
	}
	else {
		opserr << "FATAL: CTensor::Constants::IIdev() - invalid material dimension!!\n";
		exit(-1);
	}
	return IIdev;
}


// inlined ctensor functions 
inline int CTensor::length(void) const { return dataSize; }
inline int CTensor::noRows(void) const { return numRows; }
inline int CTensor::noCols(void) const { return numCols; }
inline void CTensor::Zero(void) { for (int i = 0; i < dataSize; i++) theData[i] = 0.0; }

inline double& CTensor::operator()(int row)
{
#ifdef _G3DEBUG
	if ((row < 0) || (row >= numRows)) {
		opserr << "CTensor::operator() - row " << row << " our of range [0, " << numRows - 1 << "\n";
		return data[0];
	}
	else if ((col < 0) || (col >= numCols)) {
		opserr << "CTensor::operator() - row " << col << " our of range [0, " << numCols - 1 << "\n";
		return MATRIX_NOT_VALID_ENTRY;
	}
#endif
	return theData[row];
}


inline double CTensor::operator()(int row) const
{
#ifdef _G3DEBUG
	if ((row < 0) || (row >= numRows)) {
		opserr << "CTensor::operator() - row " << row << " our of range [0, " << numRows - 1 << "\n";
		return data[0];
	}
	else if ((col < 0) || (col >= numCols)) {
		opserr << "CTensor::operator() - row " << col << " our of range [0, " << numCols - 1 << "\n";
		return MATRIX_NOT_VALID_ENTRY;
	}
#endif
	return theData[row];
}

inline double& CTensor::operator()(int row, int col)
{
#ifdef _G3DEBUG
	if ((row < 0) || (row >= numRows)) {
		opserr << "CTensor::operator() - row " << row << " our of range [0, " << numRows - 1 << "\n";
		return data[0];
	}
	else if ((col < 0) || (col >= numCols)) {
		opserr << "CTensor::operator() - row " << col << " our of range [0, " << numCols - 1 << "\n";
		return MATRIX_NOT_VALID_ENTRY;
	}
#endif
	return theData[col * numRows + row];
}


inline double CTensor::operator()(int row, int col) const
{
#ifdef _G3DEBUG
	if ((row < 0) || (row >= numRows)) {
		opserr << "CTensor::operator() - row " << row << " our of range [0, " << numRows - 1 << "\n";
		return data[0];
	}
	else if ((col < 0) || (col >= numCols)) {
		opserr << "CTensor::operator() - row " << col << " our of range [0, " << numCols - 1 << "\n";
		return MATRIX_NOT_VALID_ENTRY;
	}
#endif
	return theData[col * numRows + row];
}