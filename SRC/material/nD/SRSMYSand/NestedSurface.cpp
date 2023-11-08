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

#include "NestedSurface.h"

#include <math.h>
#include <stdlib.h>


// YieldSurface class methods
NestedSurface::NestedSurface()
{
    if (theCenter.length() != 6) {
        opserr << "FATAL:NestedSurface::NestedSurface(): center size not equal 6" << endln;
        exit(-1);
    }
    if (theCenter.getRep() != 2) {
        opserr << "FATAL! NestedSurface::NestedSurface() - center with incompatible matrix reprentation recieved! center must be contravariant\n";
        exit(-1);
    }
    if (theCenter.getOrder() != 2) {
        opserr << "FATAL! NestedSurface::NestedSurface() - center with incompatible matrix order recieved! order is corrected to 2...\n";
        theCenter.setOrder(2);
    }
}

NestedSurface::NestedSurface(const CTensor& theCenter_init, double theSize_init, double plas_modul)
    :theSize(theSize_init), theCenter(theCenter_init), plastShearModulus(plas_modul)
{
    if (theCenter.length() != 6) {
        opserr << "FATAL:NestedSurface::NestedSurface(): center size not equal 6" << endln;
        exit(-1);
    }
    if (theCenter.getRep() != 2) {
        opserr << "FATAL! NestedSurface::NestedSurface() - center with incompatible matrix reprentation recieved! center must be contravariant\n";
        exit(-1);
    }
    if (theCenter.getOrder() != 2) {
        opserr << "FATAL! NestedSurface::NestedSurface() - center with incompatible matrix order recieved! order is corrected to 2...\n";
        theCenter.setOrder(2);
    }
}

NestedSurface::~NestedSurface()
{

}

void NestedSurface::setData(const CTensor& theCenter_init, double theSize_init, double plas_modul)
{
    if (theCenter_init.length() != 6) {
        opserr << "FATAL:NestedSurface::setData(): center size not equal 6" << endln;
        exit(-1);
    }
    if (theCenter_init.getRep() != 2) {
        opserr << "FATAL! NestedSurface::setData() - center with incompatible matrix reprentation recieved! center must be contravariant\n";
        exit(-1);
    }
    theSize = theSize_init;
    theCenter = theCenter_init;
    plastShearModulus = plas_modul;
    if (theCenter.getOrder() != 2) {
        opserr << "FATAL! NestedSurface::setData() - center with incompatible matrix order recieved! order is corrected to 2...\n";
        theCenter.setOrder(2);
    }
}

void NestedSurface::setCenter(const CTensor& newCenter)
{
    if (newCenter.length() != 6) {
        opserr << "FATAL:NestedSurface::setCenter(Vector &): vector size not equal 6" << endln;
        exit(-1);
    }
    if (newCenter.getRep() != 2) {
        opserr << "FATAL! NestedSurface::setCenter() - center with incompatible matrix reprentation recieved! center must be contravariant\n";
        exit(-1);
    }
    theCenter = newCenter;
    if (theCenter.getOrder() != 2) {
        opserr << "FATAL! NestedSurface::setCenter() - center with incompatible matrix order recieved! order is corrected to 2...\n";
        theCenter.setOrder(2);
    }
}

OPS_Stream& operator<<(OPS_Stream& s, const NestedSurface& ns)
{
  s << "  theSize = " << ns.theSize << endln
     << "  theCenter = " << ns.theCenter << endln
     << "  plastShearModulus = " << ns.plastShearModulus << endln;
  return s;
}
