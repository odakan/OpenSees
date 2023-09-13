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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialFileOperations.cpp$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//
// Description: This file contains the implementation for the Repository class.

#include "MaterialFileOperations.h"

Repository::Repository(const std::string pathtofile)
{
	// Create an ifstream object to read the file
	file_stream = std::ifstream(pathtofile);

	// Check if the file opened successfully
	if (!file_stream.is_open()) {
		const char* charPtrFromStr = pathtofile.c_str();
		opserr << "WARNING: Repository() - failed to open file at: " << charPtrFromStr << "\n";
	}
}

Repository::~Repository()
{
	// Close the file when done
	file_stream.close();
}


bool Repository::readLine(std::string& line) {
	// initialize
	bool END = false;

	// Read and print the contents of the file line by line
	if (std::getline(file_stream, line)) {
		END = true;
	}
	return END;
}


Vector Repository::splitString(const std::string& input, char delimiter) {
	// split input string with given delimiter and convert items to double
	// finally stack converted items in a vector and return as Vector
	std::vector<std::string> tokens;
	std::istringstream tokenStream(input);
	std::string token;

	while (std::getline(tokenStream, token, delimiter)) {
		tokens.push_back(token);
	}

	// Get the number of elements in the vector
	size_t numElements = tokens.size();

	// Allocate Vector
	Vector data(numElements);
	int count = 0;

	for (int i = 0; i < numElements; i++)
	{
		try {
			data(count) = std::stod(tokens[i]);
			count++;
		}
		catch (const std::invalid_argument& e) {
			//opserr << "WARNING: Repository::splitString() - invalid argument: " << e.what() << "\n";
			continue;
		}
		catch (const std::out_of_range& e) {
			//opserr << "WARNING: Repository::splitString() - out of range: " << e.what() << "\n";
			continue;
		}
	}

	// only return successfully converted elements
	if (count != numElements) {
		Vector temp(count);
		for (int i = 0; i < count; i++) { temp(i) = data(i); }
		data = temp;
	}

	return data;
}