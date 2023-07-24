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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/VonMisesDMM.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef VonMisesDMM_h
#define VonMisesDMM_h

/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                           VonMisesDMM nD material                            |
 +                                                                                +
 |--------------------------------------------------------------------------------|
 |                                                                                |
 +         Authors: Onur Deniz Akan (IUSS),                                       +
 |                  Guido Camata, Enrico Spacone (UNICH)                          |
 |                  and Carlo G. Lai (UNIPV)                                      |
 |                                                                                |
 +      Istituto Universitario di Studi Superiori di Pavia          (IUSS)        +
 |		Universita degli Studi 'G. d'Annunzio' Chieti - Pescara	    (UNICH)       |
 |      Universita degli Studi di Pavia                             (UNIPV)       |
 +			                                                                      +
 |                                                                                |
 |           Email: onur.akan@iusspavia.it                                        |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- June 2022                                                    |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

 /* command
  nDMaterial VonMisesDMM $matID $rho $Kref $Gref $Pref $modN $Pc $TNYS <$HModuli> <$HParams> <-intType $type>
  where:
   $matID	: unique material tag
   $rho		: mass density
   $Kref	: reference Bulk modulus
   $Gref    : reference Shear modulus
   $Pref	: reference pressure
   $modN	: modulus update power
   $Pc		: reference yield stress (cohesion)
   $TNYS	: total number of yield surfaces (> 0 for automatic |
				< 0 for user defined surface generation)
   $Hmoduli	: if TNYS < 0, vector of secant plastic shear moduli
   $HParams	: if TNYS < 0, vector of hardening parameters (shear strain)
   -intType : type of integration scheme (optional)
				rFlag = 0 implicit, backward-euler integration scheme (default)
				rFlag = 1 impl-ex, implicit/explicit integration scheme
  Description: This file contains the class definition for VonMisesDMM.
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

#include "MultiYieldSurfaceHardeningSoftening.h"

// Associated flow von Mises material definition
// IMPL-EX VonMisesDMM internal integration with Cutting plane (forward-euler) algorithm
// Stress-space defintion

class VonMisesDMM : public MultiYieldSurfaceHardeningSoftening
{
public:
	//full constructor
	VonMisesDMM(int tag, double rho,
		double Kref, double Gref, double Pref, double modn,
		DataDrivenNestedSurfaces* theData,
		int dataDriverType, int integrationType, bool verbosity);

	//null constructor
	VonMisesDMM();

	//destructor
	~VonMisesDMM();

	// return object info
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* type);

	// parallel message passing
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// miscellaneous
	void Print(OPS_Stream& s, int flag = 0);
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

private:
	// the get methods
	Vector getShiftedDeviator(const Vector& stress, const int num_ys);	// shifted deviatoric stress tensor (Ziegler's rule)

	// yield surface operations
	double yieldFunction(const Vector& stress, const int num_ys, bool yield_stress);	// Return yield function (F) value
	Vector get_dF_dS(const Vector& stress, const int num_ys);							// Return normal to the yield surface w.r.t stress
	Vector get_dF_dA(const Vector& stress, const int num_ys);							// Return normal to the yield surface w.r.t alpha (backstress)
	Vector get_dH_dA(const Vector& stress, const int num_ys);							// Return normal to the hardening potential w.r.t alpha (backstress)
	Vector get_dP_dS(const Vector& stress, const int num_ys);							// Return normal to the plastic potential w.r.t stress

	// closest point projection methods
	double linearizedFlow(const double dlambda);
	double yieldf(const Vector& stress, const int num_ys, bool yield_stress);
	Vector Qi(const Vector& stress, const int num_ys);
	Vector Di(const Vector& stress, const int num_ys);

};
#endif