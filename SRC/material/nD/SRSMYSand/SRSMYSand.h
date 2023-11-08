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
// Based on: PDMY02 and PM4Sand materials
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
	SRSMYSand(int tag,
		int nd,
		double rho,
		double refShearModul,
		double refBulkModul,
		double frictionAng,
		double peakShearStra,
		double refPress,
		double pressDependCoe,
		double phaseTransformAngle,
		double contractionParam1,
		double contractionParam3,
		double dilationParam1,
		double dilationParam3,
		int   numberOfYieldSurf = 20,
		double* gredu = 0,
		double contractionParam2 = 5.,
		double dilationParam2 = 3.,
		double liquefactionParam1 = 1.,
		double liquefactionParam2 = 0.,
		double e = 0.6,
		double volLimit1 = 0.9,
		double volLimit2 = 0.02,
		double volLimit3 = 0.7,
		double atm = 101.,
		double cohesi = 0.1,
		double hv = 0.,
		double pv = 1.
		//bool implex = false, 
		//bool verbosity = false	
	);

	// Null Constructor
	SRSMYSand(void);

	// Copy constructor
	SRSMYSand(const SRSMYSand&);

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
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* type);
	const char* getType(void) const;
	int getOrder(void) const;

	double getRho(void) { return rhox[matN]; };
	void getBackbone(Matrix&);

	// return stress & strain
	const Vector& getStress(void);
	const Vector& getStrain(void);
	const Vector& getCommittedStress(void);
	const Vector& getStressToRecord(int numOutput); // Added by Alborz Ghofrani - UW
	const Vector& getCommittedStrain(void);

	// return tangent 
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// return material response
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

	// parallel message passing
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// miscellaneous
	void Print(OPS_Stream& s, int flag = 0);
	void setSubTag(const int tag);
	int getSubTag(void) const;
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

private:
	// operational variables
	int nD = 3;
	int nDOF = 6;
	int e2p;
	int matN;
	int subTag = 0;
	bool useImplex = false;
	bool beVerbose = true;							// be verbose about internal processes (use for debugging) (no by default)

	std::vector<NestedSurface> theSurfaces;			// NOTE: surfaces[0] is not used
	std::vector<NestedSurface> committedSurfaces;
	MaterialStateVariables sv;						// Material state variables


	double modulusFactor;
	double initPress;
	double damage;
	double check;



	CTensor trialStress;
	CTensor trialStrain;
	CTensor updatedTrialStress;
	CTensor strainRate;
	CTensor subStrainRate;



	double pressureD;
	double strainPTOcta;
	double cumuDilateStrainOcta;
	double maxCumuDilateStrainOcta;
	double cumuTranslateStrainOcta;
	double pressureDCommitted;
	double cumuDilateStrainOctaCommitted;
	double maxCumuDilateStrainOctaCommitted;
	double cumuTranslateStrainOctaCommitted;
	double maxPress;
	double maxPressCommitted;

	// PPZ
	int onPPZ; //=-1 never reach PPZ before; =0 below PPZ; =1 on PPZ; =2 above PPZ
	double PPZSize;
	double prePPZStrainOcta;
	double oppoPrePPZStrainOcta;
	int onPPZCommitted;
	double PPZSizeCommitted;
	double prePPZStrainOctaCommitted;
	double oppoPrePPZStrainOctaCommitted;
	CTensor PPZPivot;
	CTensor PPZCenter;
	CTensor PivotStrainRate;
	CTensor PPZPivotCommitted;
	CTensor PPZCenterCommitted;
	CTensor PivotStrainRateCommitted;

	// Called by constructor
	void setUpSurfaces(double*);

	void updateInternal(const bool do_implex, const bool do_tangent);

	void elast2Plast(void);
	double yieldFunc(CTensor& stress, const std::vector<NestedSurface> surfaces, int surface_num);
	void deviatorScaling(CTensor& stress, const std::vector<NestedSurface> surfaces, int surfaceNum);
	void initSurfaceUpdate(void);
	void initStrainUpdate(void);



	// Return num_strain_subincre
	int setSubStrainRate(void);
	int isLoadReversal(CTensor&);
	int isCriticalState(CTensor& stress);
	void getContactStress(CTensor& contactStress);
	void getSurfaceNormal(CTensor& stress, CTensor& normal);
	double getModulusFactor(CTensor& stress);
	void updatePPZ(CTensor& stress);
	void PPZTranslation(const CTensor& contactStress);
	double getPPZLimits(int which, CTensor& contactStress);
	double getPlasticPotential(CTensor& stress, const CTensor& surfaceNormal);
	void setTrialStress(CTensor& stress);
	double getLoadingFunc(CTensor& contact, CTensor& surfaceNormal, double* plasticPotential, int crossedSurface);
	void updateActiveSurface(void);
	void updateInnerSurface(void);
	double sdevRatio(CTensor& stress, double residualPress);



	//return 1 if stress locked; o/w return 0.
	int stressCorrection(int crossedSurface);

	// Return 1 if crossing the active surface; return 0 o/w
	int  isCrossingNextSurface(void);
	double secondOrderEqn(double A, double B, double C, int i);






private:
	// user supplied
	static int matCount;
	static int* ndmx;  //num of dimensions (2 or 3)
	static int* loadStagex;  //=0 if elastic; =1 or 2 if plastic
	static double* rhox;  //mass density
	static double* refShearModulusx;
	static double* refBulkModulusx;
	static double* frictionAnglex;
	static double* peakShearStrainx;
	static double* refPressurex;
	static double* cohesionx;
	static double* pressDependCoeffx;
	static int* numOfSurfacesx;
	static double* phaseTransfAnglex;
	static double* contractParam1x;
	static double* contractParam2x;
	static double* contractParam3x;
	static double* dilateParam1x;
	static double* dilateParam2x;
	static double* liquefyParam1x;
	static double* liquefyParam2x;
	static double* dilateParam3x;
	static double* einitx;    //initial void ratio
	static double* volLimit1x;
	static double* volLimit2x;
	static double* volLimit3x;
	static double pAtm;
	static double* Hvx;
	static double* Pvx;

	// internal
	static double* residualPressx;
	static double* stressRatioPTx;
	double* mGredu;



};
#endif