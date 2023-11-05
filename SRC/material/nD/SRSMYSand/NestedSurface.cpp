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
NestedSurface::NestedSurface() :
    theSize(0.0), theCenter(6, 2), plastShearModulus(0.0)
{

}

NestedSurface::NestedSurface(const CTensor& theCenter_init,
    double theSize_init, double plas_modul) :
    theSize(theSize_init), theCenter(theCenter_init), plastShearModulus(plas_modul)
{

}

NestedSurface::~NestedSurface()
{

}

void NestedSurface::setData(const CTensor& theCenter_init,
    double theSize_init, double plas_modul)
{
    theSize = theSize_init;
    theCenter = theCenter_init;
    plastShearModulus = plas_modul;
}

void NestedSurface::setCenter(const CTensor& newCenter)
{
    if (newCenter.length() != 6) {
        opserr << "FATAL:NestedSurface::setCenter(Vector &): vector size not equal 6" << endln;
        exit(-1);
    }

    theCenter = newCenter;
}

OPS_Stream& operator<<(OPS_Stream& s, const NestedSurface& ns)
{
  s << "  theSize = " << ns.theSize << endln
     << "  theCenter = " << ns.theCenter << endln
     << "  plastShearModulus = " << ns.plastShearModulus << endln;
  return s;
}
