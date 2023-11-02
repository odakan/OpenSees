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

// Onur Deniz Akan, Guido Camata, Carlo G. Lai, Enrico Spacone and Claudio Tamagnini - IUSS Pavia 
//
// A highly stable stress-ratio formulated, strain softening capable, multi-yield-surface model for sand
//

#ifndef SRSMYSand_h
#define SRSMYSand_h

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

 /*
  Description: This file contains the class definition for SRSMYSand.
  (1) A ZeroLengthContactASDimplex element is defined by two nodes in R^2 (x,y) or R^3 (x,y,z).
  (2) Penalty (Kn,Kt) method is used to constrain displacements in the normal direction,
	  and in the tangential direction along with a Mohr-Coulomb frictional yield surface.
  (3) In 2D and 3D, 2, 4, 6 or 12 DOFs with the possibility of having  different DOFsets at
	  each of the two nodes is supported automatically. (Example: Solid <--> Shell contact in 3D)
  (4) Optionally, the normal vector of the contact plane, orienting the contact axis can be
	  specified by the user. Otherwise, the global X vector is chosen by default.
  (5) Optionally, the user can prefer using the IMPL-EX scheme over the standard Backward Euler
	  return mapping scheme for the implicit-explicit integration of the contact material residual
	  and tangent modulus.
  (6) Currently, only small deformations are considered, however an update on large deformations
	  may be implemented in the near future.
  (7) In IMPL-EX scheme, the material tangent and residual at the current step (n+1) are computed
	  by lineary extrapolating the material parameters committed at steps (n) and (n-1). Then, after
	  the convergence criteria is achieved, extrapolated material parameters are corrected with a
	  step of traditional implicit computation from the last committed step (n). Finally, corrected
	  parameters are saved as the current step (n+1) committed variables.
  (8) In a Newton-Raphson scheme, IMPL-EX integration translates to a step-wise linear solution of the material
	  non-linearity, and a symmetric, semi-positive definite tangent tensor regardless of what the analytical
	  tangent might be, hence improving the robustness of the solution. The order of accuracy is the same with
	  the implicit (backward) euler integration scheme, however the committed error at a time step is larger.
	  An automatic time-stepping algorithm should be employed in order to control (a-priori) this error. (Oliver et al, 2008)
   References:
   Oliver, J., Huespe, A.E. and Cante, J.C. "An Implicit/Explicit Integration Scheme to
		Increase Computability of Non-linear Material and Contact/Friction Problems", Comp.
		Methods Appl. Mech. Eng., 197, 1865-1889 (2008)
 */

#include <ID.h>
#include <math.h> 
#include <float.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include <Channel.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <NDMaterial.h>
#include <Information.h>
#include <FEM_ObjectBroker.h>

#include "CTensor.h"
#include "NestedSurface.h"
#include "MaterialStateVariables.h"

class SRSMYSand : public NDMaterial {
public:
	// Full Constructor
	SRSMYSand(int tag, int classTag, double r0,
		double K0, double G0, double P0, double m0, 
		double dtype, int itype, bool verbosity);


	// Null Constructor
	SRSMYSand(void);

	// Destructor
	~SRSMYSand(void);

	// iteration control
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// set material strain
	int setTrialStrain(const Vector& strain);
	int setTrialStrain(const Vector& strain, const Vector& rate);
	int setTrialStrainIncr(const Vector& strain);
	int setTrialStrainIncr(const Vector& strain, const Vector& rate);

	// return object info
	virtual NDMaterial* getCopy(void);
	virtual NDMaterial* getCopy(const char* type);
	const char* getType(void) const;
	int getOrder(void) const;
	// used by the yield surface objects
	double getGref(void);			// return reference shear modulus
	double getPref(void);			// return reference pressure
	double getGmod(void);			// return updated shear modulus

	// return stress & strain
	const Vector& getStress(void);
	const Vector& getStrain(void);

	// return tangent 
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// return material response
	virtual Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	virtual int getResponse(int responseID, Information& matInformation);

	// parallel message passing
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// miscellaneous
	virtual void Print(OPS_Stream& s, int flag = 0);
	virtual int setParameter(const char** argv, int argc, Parameter& param);
	virtual int updateParameter(int responseID, Information& eleInformation);

private:


};
#endif