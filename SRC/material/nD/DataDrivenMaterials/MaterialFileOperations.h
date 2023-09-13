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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MaterialFileOperations.h$
// $Revision: 1.0 $
// $Date: 2023-XX-XX XX:XX:XX $

#ifndef MaterialFileOperations_h
#define MaterialFileOperations_h

// Written by:	Onur Deniz Akan		(onur.akan@iusspavia.it)
//				Guido Camata      
//				Enrico Spacone
//				Carlo G. Lai
//
// Created in:	September 2023
//

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <Vector.h>


class Repository {
public:
	// null constructor
	Repository(void) = default;

	// full constructor
	Repository(const std::string pathtofile);

	// destructor
	~Repository();

	// file operations
	bool readLine(std::string& line);
	Vector splitString(const std::string& input, char delimiter);

private:
	std::ifstream file_stream;
	int current_line_number = 0;
	char* current_line = nullptr;


};
#endif