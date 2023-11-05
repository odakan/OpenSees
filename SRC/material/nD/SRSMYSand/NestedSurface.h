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

// Written: Onur Deniz Akan, Guido Camata, Carlo G. Lai, Enrico Spacone and Claudio Tamagnini
// Created: 11/23
// Based on the MultiYieldSurface object of PDMY02
//
// A highly stable stress-ratio formulated, strain softening capable, multi-yield-surface model for sand

#ifndef NestedSurface_h
#define NestedSurface_h

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                             SRSMYSand nD material                              |
 +                                                                                +
 |--------------------------------------------------------------------------------|
 |                                                                                |
 +         Authors: Onur Deniz Akan (IUSS),                                       +
 |                  Guido Camata, Enrico Spacone (UNICH),                         |
 |                  Carlo G. Lai (UNIPV) and                                      |
 +                  Claudio Tamagnini (UNIPG)                                     +
 |                                                                                |
 |         Istituto Universitario di Studi Superiori di Pavia (IUSS)              |
 +		   Universita degli Studi Chieti - Pescara	          (UNICH)             +
 |         Universita degli Studi di Pavia                    (UNIPV)             |
 |         Università degli Studi di Perugia                  (UNIPG)             |
 +			                                                                      +
 |                                                                                |
 |         Email: onur.akan@iusspavia.it                                          |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- November 2023                                                |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

#include <float.h>
#include "CTensor.h"

#define UP_LIMIT    1.0e+30
#define LOW_LIMIT   20.*DBL_EPSILON

class NestedSurface {
public:
    //constructors
    NestedSurface();
    NestedSurface(const CTensor& center_init, double size_init,
        double plas_modul);
    ~NestedSurface();
    void setData(const CTensor& center_init, double size_init,
        double plas_modul);
    const CTensor& center() const { return theCenter; }
    double size() const { return theSize; }
    double modulus() const { return plastShearModulus; }
    void  setCenter(const CTensor& newCenter);

    friend OPS_Stream& operator<<(OPS_Stream& s, const NestedSurface& C);

protected:

private:
    double theSize;
    CTensor theCenter;
    double plastShearModulus;

};
#endif