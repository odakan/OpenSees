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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialResponseDatabase.cpp$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//
// Description: This file contains the implementation for the Database class.

#include "MaterialResponseDatabase.h"

Database::Database(const int tag, const char* dir, const char* fname)
{
	// if database path is present
	if (dir != nullptr && fname != nullptr) {		
		// recieve associated material ID and path
		matID = tag;
		std::string directory(dir);
		std::string mainfile(fname);
		fullpath = directory + mainfile;

		// initialize
		std::string filename;
		std::string line;
		Repository* main = nullptr;
		Repository* file = nullptr;
		Vector Info;
		Vector Data;

		//make sure the datapts vector is empty
		datapts.clear();
		length = 0;

		// load the database into the RAM
			// inform user about process
		opserr << "WARNING: Database() - please wait while the database is being loaded...\n";
		opserr << "          -> material ID: " << tag << "\n";
		opserr << "          -> database at: " << fullpath.c_str() << "\n";
		// open main file
		main = new Repository(directory + mainfile);
		// read file names
		while (main->readLine(filename)) {
			// open file with read name
			file = new Repository(directory + filename);
			Info = main->splitString(filename, '_');
			while (file->readLine(line)) {
				// split line and extract the data
				Data = file->splitString(line, ' ');
				// create a data container
				std::unique_ptr<DataPoint> point = std::make_unique<DataPoint>(Data, Info);
				//DataPoint* point = new DataPoint(Data, Info);
				// check dimension
				if (nDim == 0) nDim = point->getDim();  // assume databse has the dimension of the first point
				if (nDim != point->getDim()) {
					opserr << "WARNING: Database() - encountered with an incompatible data point! Skipping the point...";
					continue;
				}
				// compute some statistics
				if (point->Toct > Toct_max) Toct_max = point->Toct;
				if (point->Toct < Toct_min) Toct_min = point->Toct;
				if (point->Goct > Goct_max) Goct_max = point->Goct;
				if (point->Goct < Goct_min) Goct_min = point->Goct;
				if (point->Pavg > Pavg_max) Pavg_max = point->Pavg;
				if (point->Pavg < Pavg_min) Pavg_min = point->Pavg;
				if (point->Evol > Evol_max) Evol_max = point->Evol;
				if (point->Evol < Evol_min) Evol_min = point->Evol;
				length++;
				// store the pointer to the container
				datapts.push_back(std::move(point));
				//datapts.push_back(point);
			}
			nTests++;
		}

		// Inform user about the sucessful outcome
		opserr << "WARNING: Database() - database successfully loaded! -> number of points: " << length << "\n\n";

		// free memory before continue
		delete main;
		delete file;
	}
	else {
		// if the path is not present, do nothing. i.e. use the default values...
	}
}

Database::Database(const Database& other):
	nDim(other.nDim), matID(other.matID), length(other.length),
	nTests(other.nTests), fullpath(other.fullpath),
	Goct_max(other.Goct_max), Goct_min(other.Goct_min), Toct_max(other.Toct_max),
	Toct_min(other.Toct_min), Pavg_max(other.Pavg_max), Pavg_min(other.Pavg_min),
	Evol_max(other.Evol_max), Evol_min(other.Evol_min)
{
	// Copy the data points
	datapts.reserve(other.datapts.size());
	for (const auto& ptr : other.datapts) {
		datapts.push_back(std::make_unique<DataPoint>(*ptr));
	}

	opserr << "WARNING: Database() - copy constructer was called!\n";
}

Database::~Database()
{
	// smart pointer "unique_ptr" deletes datapts objects when out of scope
	
	//if (!datapts.empty()) {
	//	for each (auto ptr in datapts)
	//	{
	//		delete ptr;
	//	}
	//	datapts.clear();
	//}
}

Database& Database::operator=(const Database& other) {
	if (this != &other) {
		// Copy non-pointer members
		nDim = other.nDim;
		matID = other.matID;
		length = other.length;
		nTests = other.nTests;
		fullpath = other.fullpath;
		Goct_max = other.Goct_max;
		Goct_min = other.Goct_min;
		Toct_max = other.Toct_max;
		Toct_min = other.Toct_min;
		Pavg_max = other.Pavg_max;
		Pavg_min = other.Pavg_min;
		Evol_max = other.Evol_max;
		Evol_min = other.Evol_min;

		// Clear the current data in datapts
		datapts.clear();

		// Copy or move the data from the other object
		datapts.reserve(other.datapts.size());
		for (const auto& ptr : other.datapts) {
			datapts.push_back(std::make_unique<DataPoint>(*ptr));
		}
	}
	opserr << "WARNING: Database::operator=() - copy assignment operator was called!\n";
	return *this;
}

OPS_Stream& operator<<(OPS_Stream& s, const Database& obj) {

	s << "\nDatabase " << obj.matID << " summary\n";
	s << "---------------------------------------------\n";
	s << "Database nD          : " << obj.nDim << "\n";
	s << "Number of tests      : " << obj.nTests << "\n";
	s << "Number of points     : " << obj.length << "\n";
	s << "Max oct. shear stress: " << obj.Toct_max << "\n";
	s << "Min oct. shear stress: " << obj.Toct_min << "\n";
	s << "Max oct. shear strain: " << obj.Goct_max << "\n";
	s << "Min oct. shear strain: " << obj.Goct_min << "\n";
	s << "Max average pressure : " << obj.Pavg_max << "\n";
	s << "Min average pressure : " << obj.Pavg_min << "\n";
	s << "Max volumetric strain: " << obj.Evol_max << "\n";
	s << "Min volumetric strain: " << obj.Evol_min << "\n";
	s << "Database main file   : " << obj.fullpath.c_str() << "\n\n";

	return s << endln;
}

// operational functions
int Database::size(void) const { return length; }
int Database::getTag(void) const { return matID; }
int Database::getDim(void) const { return nDim; }

// get methods
double Database::getVoidRatio(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->getVR();
	}
	else {
		opserr << "FATAL - Database::getVoidRatio() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

double Database::getDilatancy(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->Bsec;
	}
	else {
		opserr << "FATAL - Database::getDilatancy() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

double Database::getShearModulus(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->Gsec;
	}
	else {
		opserr << "FATAL - Database::getShearModulus() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

double Database::getVolumetricStrain(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->Evol;
	}
	else {
		opserr << "FATAL - Database::getVolumetricStrain() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

double Database::getOctahedralStress(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->Toct;
	}
	else {
		opserr << "FATAL - Database::getOctahedralStress() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

double Database::getOctahedralStrain(const int index) const {
	double value = 0.0;
	if (index >= 0 && index < length) {
		value = datapts[index]->Goct;
	}
	else {
		opserr << "FATAL - Database::getOctahedralStrain() - requested index exceeds the length of the database!\n";
		exit(-1);
	}
	return value;
}

// permitted queries
int Database::seek(const char* token, const Vector& tensor) const {
	// searh for the point with the closest desired tensor
	
	// initialize
	int index = -1;
	int size = tensor.Size();
	double TOL = 1-4;
	double closest_distance = 1;
	Vector& data = Vector(size);

	// check compatibility
	if (size != 3 * (nDim - 1)) { 
		opserr << "FATAL - Database::seek() - material dimension does not match the database!\n";
		exit(-1);
	}

	// compute norms
	for (int i = 0; i < length; i++) {
		double distance = 0.0;
		// fetch data tensor
		if (strcmp(token, "stress") == 0) {
			TOL = 1;
			closest_distance = 1e4;
			data = datapts[i]->getStress();
		}
		else if (strcmp(token, "strain") == 0) {
			TOL = 1e-4;
			closest_distance = 1;
			data = datapts[i]->getStrain();
		}
		else if (strcmp(token, "fabric") == 0) {
			TOL = 1e-4;
			closest_distance = 1;
			data = datapts[i]->getFabric();
		}
		else {
			opserr << "FATAL - Database::seek() - the desired tensor is not available in the database!\n";
			exit(-1);
		}
		// compute the tensor-to-tensor distance (in Voigt notation)
		for (int j = 0; j < size; j++) {
			if (j < nDim) {
				distance += ((data(j) * data(j)) - (tensor(j) * tensor(j)));
			}
			else {
				distance += 2 * ((data(j) * data(j)) - (tensor(j) * tensor(j)));
			}
		}
		distance = sqrt(distance);
		// keep the minimum distance
		if (distance < closest_distance) {
			closest_distance = distance;
			// save index
			index = i;
		}
	}

	// return index
	if (index < 0 || index >= length) {
		opserr << "FATAL - Database::seek() - search algorithm failed!\n";
		exit(-1);
	}
	else if (closest_distance > TOL) {
		opserr << "WARNING - Database::seek() - failed to find a sufficiently close data point!\n";
	}

	return index;
}

int Database::seek(const double I1, const double J2, const double J3) const {
	// search for the point with the closest strain tensor invariants

	// initialize
	int index = -1;
	double TOL = 1 - 4;
	double distance = 0.0;
	double closest_distance = 1.0;

	// compute norms
	for (int i = 0; i < length; i++) {
		distance = sqrt(((datapts[index]->I1 * datapts[index]->I1) - (I1 * I1)) +
						((datapts[index]->J2 * datapts[index]->J2) - (J2 * J2)) +
						((datapts[index]->J3 * datapts[index]->J3) - (J3 * J3)));
		// find the minimum distance
		if (distance < closest_distance) {
			closest_distance = distance;
			// save index
			index = i;
		}
	}

	// return index
	if (index < 0 || index >= length) {
		opserr << "FATAL - Database::seek() - search algorithm failed!\n";
		exit(-1);
	}
	else if (closest_distance > TOL) {
		opserr << "WARNING - Database::seek() - failed to find a sufficiently close data point!\n";
	}

	return index;
}

int Database::seek(const double pavg, const char* token, const Vector& tensor) const {
	// search for the point with the closest first stress invariant and the strain deviator tensor
	
	// initialize
	int index = -1;
	int size = tensor.Size();
	double Evol = 0.0;
	double TOL = 0.1;
	double closest_distance = 1;
	Vector data_dev(size);
	Vector tensor_dev(size);
	Vector& data = Vector(size);

	// check compatibility
	if (size != 3 * (nDim - 1)) {
		opserr << "FATAL - Database::seek() - material dimension does not match the database!\n";
		exit(-1);
	}

	// compute the strain deviator tensor
	if (size == 3) {
		Evol = tensor(0) + tensor(1);
		tensor_dev = tensor; tensor_dev(0) -= Evol * 0.5; tensor_dev(1) -= Evol * 0.5;
	}
	else {
		Evol = tensor(0) + tensor(1) + tensor(2);
		tensor_dev = tensor; tensor_dev(0) -= Evol * 1.0 / 3.0; tensor_dev(1) -= Evol * 1.0 / 3.0; tensor_dev(2) -= Evol * 1.0 / 3.0;
	}
	
	// compute norms
	for (int i = 0; i < length; i++) {
		// only check points with similar average pressure
		if (abs(datapts[i]->Pavg - pavg) <= TOL) {
			double distance = 0.0;
			// fetch data tensor
			if (strcmp(token, "stress") == 0) {
				TOL = 1;
				closest_distance = 1e4;
				data = datapts[i]->getStress();
			}
			else if (strcmp(token, "strain") == 0) {
				TOL = 1e-4;
				closest_distance = 1;
				data = datapts[i]->getStrain();
				opserr << "Database::seek() - works until here!\n";
			}
			else if (strcmp(token, "fabric") == 0) {
				TOL = 1e-4;
				closest_distance = 1;
				data = datapts[i]->getFabric();
			}
			else {
				opserr << "FATAL - Database::seek() - the desired tensor is not available in the database!\n";
				exit(-1);
			}
			// compute the strain deviator tensor
			if (size == 3) {
				Evol = data(0) + data(1);
				data_dev = data; data_dev(0) -= Evol * 0.5; data_dev(1) -= Evol * 0.5;
			}
			else {
				Evol = data(0) + data(1) + data(2);
				data_dev = data; data_dev(0) -= Evol * 1.0 / 3.0; data_dev(1) -= Evol * 1.0 / 3.0; data_dev(2) -= Evol * 1.0 / 3.0;
			}
			// compute the tensor-to-tensor distance (in Voigt notation)
			for (int j = 0; j < size; j++) {
				if (j < nDim) {
					distance += ((data_dev(j) * data_dev(j)) - (tensor_dev(j) * tensor_dev(j)));
				}
				else {
					distance += 2 * ((data_dev(j) * data_dev(j)) - (tensor_dev(j) * tensor_dev(j)));
				}
			}
			distance = sqrt(distance);
			// keep the minimum distance
			if (distance < closest_distance) {
				closest_distance = distance;
				// save index
				index = i;
			}
		}

	}

	// return index
	if (index < 0 || index >= length) {
		opserr << "FATAL - Database::seek() - search algorithm failed!\n";
		exit(-1);
	}
	else if (closest_distance > TOL) {
		opserr << "WARNING - Database::seek() - failed to find a sufficiently close data point!\n";
	}

	return index;
}
