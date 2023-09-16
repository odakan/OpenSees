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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialResponseDatabase.h$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

#ifndef MaterialResponseDatabase_h
#define MaterialResponseDatabase_h

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//

#include <math.h>
#include <Vector.h>
#include "MaterialDataPoint.h"
#include "MaterialFileOperations.h"

class Database {
public:
	// constructors
	Database(void) = default;					// null constructor
	Database(const Database& other);			// copy constructor
	Database(const int tag, const char* dir,	// full constructor 
		const char* main);	

	// destructor
	~Database();

	// operator overloading
	Database& operator= (const Database& other);
	friend OPS_Stream& operator<<(OPS_Stream& s, const Database& obj);

	// operational functions
	int size(void) const;
	int getTag(void) const;
	int getDim(void) const;

	// get methods
	double getVoidRatio(const int index) const;
	double getMeanPressure(const int index) const;
	double getVolumetricStrain(const int index) const;
	double getOctahedralStress(const int index) const;
	double getOctahedralStrain(const int index) const;

	// permitted queries
	int seek(const double pavg, const double gamma) const;
	int seek(const char* token, const Vector& tensor) const;
	int seek(const double I1, const double J2, const double J3) const;
	int seek(const double pavg, const char* token, const Vector& tensor) const;

private:
	// variables
	int nDim = 0;
	int matID = 0;
	int length = 0;
	int nTests = 0;

	// path to main file
	std::string fullpath;

	// min-max data
	double Goct_max = SMALL_VALUE;
	double Goct_min = LARGE_VALUE;
	double Toct_max = SMALL_VALUE;
	double Toct_min = LARGE_VALUE;
	double Pavg_max = -LARGE_VALUE;
	double Pavg_min = LARGE_VALUE;
	double Evol_max = SMALL_VALUE;
	double Evol_min = LARGE_VALUE;

	// data points
	std::vector<std::unique_ptr<DataPoint>> datapts;
};
#endif