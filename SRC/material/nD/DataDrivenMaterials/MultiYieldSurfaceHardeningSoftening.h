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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DataDrivenMaterials/MultiYieldSurfaceHardeningSoftening.h$
// $Revision: 1.0 $
// $Date: 2022-XX-XX XX:XX:XX $

#ifndef MultiYieldSurfaceHardeningSoftening_h
#define MultiYieldSurfaceHardeningSoftening_h

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                MultiYieldSurfaceHardeningSoftening nD material                 |
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

 /* 
  Description: This file contains the class definition for ZeroLengthContactASDimplex.
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
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <NDMaterial.h>
#include <Information.h>
#include <FEM_ObjectBroker.h>
#include "YieldSurfacePackage.h"
#include "MaterialStateVariables.h"
#include "MaterialTensorOperations.h"
#include "DataDrivenNestedSurfaces.h"


class MultiYieldSurfaceHardeningSoftening : public NDMaterial {
public:
	// Full Constructor
	MultiYieldSurfaceHardeningSoftening(int tag, int classTag, double r0,
		double K0, double G0, double P0, double m0, int T0,
		DataDrivenNestedSurfaces* ys, int dtype, int itype, bool verbosity);


	// Null Constructor
	MultiYieldSurfaceHardeningSoftening(void);

	// Destructor
	~MultiYieldSurfaceHardeningSoftening(void);

	// iteration control
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// set material strain
	int setTrialStrain(const Vector& strain);
	int setTrialStrain(const Vector& strain, const Vector& rate);
	int setTrialStrainIncr(const Vector& strain);
	int setTrialStrainIncr(const Vector& strain, const Vector& rate);
		// used by the surface generator
	void setTheSurfaces(YieldSurfacePackage* theSurfaces);

	// return object info
	virtual NDMaterial* getCopy(void) = 0;
	virtual NDMaterial* getCopy(const char* type) = 0;
	const char* getType(void) const;
	int getDataDriver(void);		// return yield surface data driver preference
	int getOrder(void) const;
	Vector getState(void);
		// used by the yield surface objects
	double getGref(void);			// return reference shear modulus
	double getPref(void);			// return reference pressure
	double getTNYS(void);			// return total number of yield surfaces

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
	virtual int sendSelf(int commitTag, Channel& theChannel) = 0;
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) = 0;

	// miscellaneous
	virtual void Print(OPS_Stream& s, int flag = 0);
	virtual int setParameter(const char** argv, int argc, Parameter& param);
	virtual int updateParameter(int responseID, Information& eleInformation);

protected:
	// material parameters & constants
	int nOrd = 6;									// Material order (dimension) [3 or 6] (3D by default)
	int TNYS = 0;									// Total number of yield surfaces
	double rho = 0.0;								// Mass density
	double Kref = 0.0;								// Reference bulk modulus
	double Gref = 0.0;								// Reference shear modulus
	double Pref = 0.0;								// Reference pressure
	double Modn = 0.0;								// Modulus update power n
	MaterialStateVariables sv;						// Material state variables
	YieldSurfacePackage ys;							// Nested yield surface package (free memory before exit)
	DataDrivenNestedSurfaces* theData = nullptr;	// Pointer to the material response database (the glorious lookup table) (free memory before exit)

	// operational paramaters
	bool beVerbose = false;							// be verbose about internal processes (use for debugging) (no by default)

	// solution options
	bool use_data_driven_surface = false;			// yield surface flag [true: use data driven yield surfaces, false: use hyperbolic backbone]
	bool use_online_approach = false;				// yield surface update approach flag [true: online, false: offline]
	bool use_implex = false;						// integration type flag: impl-ex or implicit (latter by default)
	bool use_numerical_tangent = false;				// implicit tangent flag: numeric or elastoplastic (latter by default)
	int materialStage = 0;							// use updateMaterialStage [0 = linear elastic, 1 = elastoplastic, 2 = nonlinear elastic]
	int solution_strategy = 0;						// [0 = cutting plane algorithm, 1 = closest point projection]

protected:
	// the get methods
	double getK(void);																	// bulk modulus
	double getG(void);																	// shear modulus
	double getSizeYS(const int num_yield_surface);										// committed yield surface radius
	double getMeanStress(const Vector& stress);											// mean volumetric stress
	const Vector getStressVector(void);													// committed stress vector
	const Vector getStrainVector(void);													// committed strain vector
	const Vector getPlasticStrainVector(void);											// committed plastic strain vector
	const Vector getStressDeviator(const Vector& stress, int num_yield_surface);		// deviatoric stress tensor (Ziegler's rule)

	// yield surface operations (must be overloaded by the sub-class)
	virtual double yieldFunction(const Vector& stress, const int num_yield_surface, 	// Return yield function (F) value
		bool yield_stress) = 0;
	virtual Vector get_dF_dS(const Vector& stress, const int num_yield_surface) = 0;	// Return normal to the yield surface w.r.t stress
	virtual Vector get_dF_dA(const Vector& stress, const int num_yield_surface) = 0;	// Return normal to the yield surface w.r.t alpha (backstress)
	virtual Vector get_dH_dA(const Vector& stress, const int num_yield_surface) = 0;	// Return normal to the hardening potential w.r.t alpha (backstress)
	virtual Vector get_dP_dS(const Vector& stress, const int num_yield_surface) = 0;	// Return normal to the plastic potential w.r.t stress

	// material internal operations
		// update methods
	void updateStress(Vector& stress, const double lambda, const int num_yield_surface);			// Update stress
	void updateModulus(const Vector& stress, const int num_yield_surface);														// Update elastic modulus
	void updateInternal(bool do_implex, bool do_tangent);											// Update materal internal state
	void updatePlasticStrain(Vector& pstrain, const Vector& stress,									// Update plastic strain
		const double lambda, const int num_yield_surface);
	void updateFailureSurface(const Vector& stress);												// Update the final yield surface (alpha)
	void updateInnerYieldSurfaces(int num_yield_surface, const Vector& stress);						// Update the inner yield surfaces (alpha)
	void updateCurrentYieldSurface(int num_yield_surface, double lambda, const Vector& stress);		// Update current active yield surface (alpha)

		// compute methods
	void computeElastoplasticTangent(int num_yield_surface, const Vector& stress);				// Compute the algorithmic tangent operator
	double computePlasticLoadingFunction(const Vector& stress,									// Compute the consistency parameter (by linearizing F)
		const double yf_value, const int num_yield_surface);

		// return-mapping
	int cuttingPlaneAlgorithm(const Vector& sigma_trial);
	int closestPointProjection(const Vector& sigma_trial);

		// correct methods
	void correctStress(Vector& stress, const double lambda1, const double lambda2,				// After overshooting, correct the overshooting stress
		const int num_yield_surface);						
	void correctPlasticStrain(Vector& pstrain, const Vector& stress,							// After overshooting, correct the plastic strain
		const double lambda1, const double lambda2, const int num_yield_surface);	
	void correctPlasticLoadingFunction(Vector& stress, double& lambda1,							// After overshooting, correct the consistency parameter
		double& lambda2, const int num_yield_surface);

		// root search algortihm
	double zbrentStress(const Vector& start_stress, const Vector& end_stress, 
		const int num_yield_surface, const double x1, const double x2, const double tol);

};
#endif