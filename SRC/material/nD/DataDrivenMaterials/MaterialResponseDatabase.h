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
	Database(void) = default;							// null constructor
	Database(const Database&) = default;				// copy constructor
	Database(const int tag, const char* dir, const char* main);		// full constructor

	// destructor
	~Database();

	// operator overloading
	Database& operator= (const Database&) = default;
	friend OPS_Stream& operator<<(OPS_Stream& s, const Database& obj);

	// operational functions
	int size(void);
	int getTag(void);
	int getDim(void);

	// get methods
	double getVoidRatio(const int index);
	double getDilatancy(const int index);
	double getShearModulus(const int index);
	double getVolumetricStrain(const int index);
	double getOctahedralStress(const int index);
	double getOctahedralStrain(const int index);

	// permitted queries
	int seek(const char* token, const Vector& tensor);
	int seek(const double I1, const double J2, const double J3);
	int seek(const double pavg, const char* token, const Vector& tensor);

private:
	// variables
	int nDim = 0;
	int matID = 0;
	int length = 0;
	int nTests = 0;

	// path to main file
	std::string fullpath;

	// min-max data
	double Goct_max = 0.0;
	double Goct_min = 0.0;
	double Toct_max = 0.0;
	double Toct_min = 0.0;
	double Pavg_max = 0.0;
	double Pavg_min = 0.0;
	double Evol_max = 0.0;
	double Evol_min = 0.0;

	// data points
	std::vector<DataPoint*> datapts;
};
#endif