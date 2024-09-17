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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PDMY02Implex.h,v $
// $Revision: 1.0 $
// $Date: 2024-01-05 12:00:00 $

// Written by:	    Onur Deniz Akan		(onur.akan@iusspavia.it)
// Based on:        ZHY's PressureDependMultiYield02.h
// Created:         December 2023
// Last Modified:
//
// Description: This file contains the implementation for the PDMY02Implex function.

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                           PDMY02Implex nD material                             |
 |                                                                                |
 +--------------------------------------------------------------------------------+
 |                                                                                |
 |         Authors: Onur Deniz Akan (IUSS),                                       |
 +                  Guido Camata, Enrico Spacone (UNICH),                         +
 |                  Carlo G. Lai (UNIPV) and Claudio Tamagnini (UNIPG)            |
 |                                                                                |
 +      Istituto Universitario di Studi Superiori di Pavia          (IUSS)        +
 |		Universita degli Studi 'G. d'Annunzio' Chieti - Pescara	    (UNICH)       |
 |      Universita degli Studi di Pavia                             (UNIPV)       |
 +		Universita degli Studi di Perugia                           (UNIPG)       +
 |                                                                                |
 |           Email: onur.akan@iusspavia.it                                        |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- December 2023                                                |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

#ifndef PDMY02Implex_h
#define PDMY02Implex_h

#include <NDMaterial.h>
#include <Matrix.h>
#include "soil/T2Vector.h"

class MultiYieldSurface;

class PDMY02Implex : public NDMaterial
{
public:
    // Initialization constructor
    PDMY02Implex(int tag,
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
        int   numberOfYieldSurf,
        double* gredu,
        double contractionParam2,
        double dilationParam2,
        double liquefactionParam1,
        double liquefactionParam2,
        double e,
        double volLimit1,
        double volLimit2,
        double volLimit3,
        double atm,
        double cohesi,
        double hv,
        double pv,
        int implexFlag,
        int verbosityFlag);

    // Default constructor
    PDMY02Implex();

    // Copy constructor
    PDMY02Implex(const PDMY02Implex&);

    // Destructor: clean up memory storage space.
    virtual ~PDMY02Implex();

    double getRho(void) { return rhox[matN]; };

    // Sets the values of the trial strain tensor.
    int setTrialStrain(const Vector& strain);

    // Sets the values of the trial strain and strain rate tensors.
    int setTrialStrain(const Vector& v, const Vector& r);

    int setTrialStrainIncr(const Vector& v);
    int setTrialStrainIncr(const Vector& v, const Vector& r);

    // Calculates current tangent stiffness.
    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);

    void getBackbone(Matrix&);

    // Calculates the corresponding stress increment (rate), for a given strain increment.
    const Vector& getStress(void);
    const Vector& getStrain(void);
    const Vector& getCommittedStress(bool isImplex = false);
    const Vector& getStressToRecord(int numOutput, bool isImplex = false); // Added by Alborz Ghofrani - UW
    const Vector& getCommittedStrain(void);
    const Vector& getImplexSateVariables(void);

    // Accepts the current trial strain values as being on the solution path, and updates
    // all model parameters related to stress/strain states. Return 0 on success.
    int commitState(void);

    // Revert the stress/strain states to the last committed states. Return 0 on success.
    int revertToLastCommit(void);

    int revertToStart(void) { return 0; }

    // Return an exact copy of itself.
    NDMaterial* getCopy(void);

    // Return a copy of itself if "code"="PDMY02Implex", otherwise return null.
    NDMaterial* getCopy(const char* code);

    // Return the string "PDMY02Implex".
    const char* getType(void) const;

    // Return ndm.
    int getOrder(void) const;

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);
    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& matInformation);
    void Print(OPS_Stream& s, int flag = 0);
    //void setCurrentStress(const Vector stress) { currentStress=T2Vector(stress); }
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int responseID, Information& eleInformation);


    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

protected:

private:
    // user supplied
    static int matCount;
    static int* ndmx;                       //num of dimensions (2 or 3)
    static int* loadStagex;                 //=0 if elastic; =1 or 2 if plastic
    static double* rhox;                    //mass density
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
    static double* einitx;                  //initial void ratio
    static double* volLimit1x;
    static double* volLimit2x;
    static double* volLimit3x;
    static double pAtm;
    static double* Hvx;
    static double* Pvx;
    static bool* doImplex;                  // solve implicit by default
    static bool* beVerbose;                 // give internal solution info

    // internal
    static double* residualPressx;
    static double* stressRatioPTx;
    static Matrix theTangent;
    double* mGredu;

    // implex state variables
    double lambda = 0.0;                    // step n+1 plastic multiplier
    double lambda_commit = 0.0;             // step n   plastic multiplier 
    double lambda_commit_old = 0.0;         // step n-1 plastic multiplier
    double chi = 0.0;                       // step n+1 plastic distortion deformability
    double chi_commit = 0.0;                // step n   plastic distortion deformability
    double chi_commit_old = 0.0;            // step n-1 plastic distortion deformability
    double iota = 0.0;
    double iota_commit = 0.0;
    double iota_commit_old = 0.0;
    double ksi = 0.0;
    double ksi_commit = 0.0;
    double ksi_commit_old = 0.0;
    double kappa = 0.0;                     // step n+1 normalized plastic dilatancy
    double kappa_commit = 0.0;              // step n   normalized plastic dilatancy
    double kappa_commit_old = 0.0;          // step n-1 normalized plastic dilatancy
    double dtime_n = 0.0;                   // time factor
    double dtime_n_commit = 0.0;            // committed time factor
    bool dtime_is_user_defined = false;
    bool dtime_first_set = false;

    double maxCumuDilateStrainOcta_bar = 0.0;
    double maxCumuDilateStrainOcta_bar_commit = 0.0;
    double cumuDilateStrainOcta_bar = 0.0;
    double cumuDilateStrainOcta_bar_commit = 0.0;
    bool isCS = false;
    bool isCS_commit = false;
    double shearLoading_bar = 0.0;
    double shearLoading_bar_commit = 0.0;
    double currentRatio_bar = 0.0;
    double currentRatio_bar_commit = 0.0;
    double trialRatio_bar = 0.0;
    double trialRatio_bar_commit = 0.0;
    double contactStressVolume_bar = 0.0;
    double contactStressVolume_bar_commit = 0.0;
    double contactRatio_bar = 0.0;
    double contactRatio_bar_commit = 0.0;
    int onPPZ_bar = 0;
    int onPPZ_bar_commit = 0;

    // new results
    T2Vector currentStressImplex;               // extrapolated stress as a material result

    int matN;
    int e2p;
    MultiYieldSurface* theSurfaces; // NOTE: surfaces[0] is not used
    MultiYieldSurface* committedSurfaces;
    int    activeSurfaceNum;
    int    committedActiveSurf;
    double modulusFactor;
    double initPress;
    double damage;
    double check;
    T2Vector currentStress;
    T2Vector trialStress;
    T2Vector updatedTrialStress;
    T2Vector currentStrain;
    T2Vector strainRate;
    static T2Vector subStrainRate;

    double pressureD;
    int onPPZ; //=-1 never reach PPZ before; =0 below PPZ; =1 on PPZ; =2 above PPZ
    double strainPTOcta;
    double PPZSize;
    double cumuDilateStrainOcta;
    double maxCumuDilateStrainOcta;
    double cumuTranslateStrainOcta;
    double prePPZStrainOcta;
    double oppoPrePPZStrainOcta;
    static T2Vector trialStrain;
    T2Vector PPZPivot;
    T2Vector PPZCenter;
    Vector PivotStrainRate;

    double pressureDCommitted;
    int onPPZCommitted;
    double PPZSizeCommitted;
    double cumuDilateStrainOctaCommitted;
    double maxCumuDilateStrainOctaCommitted;
    double cumuTranslateStrainOctaCommitted;
    double prePPZStrainOctaCommitted;
    double oppoPrePPZStrainOctaCommitted;
    T2Vector PPZPivotCommitted;
    T2Vector PPZCenterCommitted;
    Vector PivotStrainRateCommitted;
    static Matrix workM66;
    static Vector workV6;
    static T2Vector workT2V;
    double maxPress;

    // move stress computations here
    int implicitSress(void);
    int updateImplex(void);     // explicit stage

    void elast2Plast(void);
    // Called by constructor
    void setUpSurfaces(double*);
    double yieldFunc(const T2Vector& stress, const MultiYieldSurface* surfaces,
        int surface_num);
    void deviatorScaling(T2Vector& stress, const MultiYieldSurface* surfaces,
        int surfaceNum);
    void initSurfaceUpdate(void);
    void initStrainUpdate(void);

    // Return num_strain_subincre
    int setSubStrainRate(void);
    int isLoadReversal(const T2Vector&);
    int isCriticalState(const T2Vector& stress);
    void getContactStress(T2Vector& contactStress);
    void getSurfaceNormal(const T2Vector& stress, T2Vector& normal);
    double getModulusFactor(T2Vector& stress);
    void updatePPZ(const T2Vector& stress);
    void PPZTranslation(const T2Vector& contactStress);
    double getPPZLimits(int which, const T2Vector& contactStress);
    double getPlasticPotential(const T2Vector& stress, const T2Vector& surfaceNormal);
    void setTrialStress(T2Vector& stress);
    double getLoadingFunc(const T2Vector& contact, const T2Vector& surfaceNormal,
        double* plasticPotential, int crossedSurface);
    //return 1 if stress locked; o/w return 0.
    int stressCorrection(int crossedSurface);
    void updateActiveSurface(void);
    void updateInnerSurface(void);

    // Return 1 if crossing the active surface; return 0 o/w
    int  isCrossingNextSurface(void);
};

#endif