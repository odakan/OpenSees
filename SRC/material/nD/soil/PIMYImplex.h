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

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PIMYImplex.h,v $
// $Revision: 1.0 $
// $Date: 2024-01-05 12:00:0 $

// Written by:	    Onur Deniz Akan		(onur.akan@iusspavia.it)
// Based on:        ZHY's PressureIndependMultiYield.h
// Created:         December 2023
// Last Modified:
//
// Description: This file contains the implementation for the PIMYImplex function.

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                            PIMYImplex nD material                              |
 |                                                                                |
 +--------------------------------------------------------------------------------+
 |                                                                                |
 |         Authors: Onur Deniz Akan (IUSS),                                       |
 +                  Guido Camata, Enrico Spacone (UNICH)                          +
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
 |  Created       -- December 2023                                                |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

#ifndef PIMYImplex_h
#define PIMYImplex_h

#include <NDMaterial.h>
#include "soil/T2Vector.h"
#include <Matrix.h>

class MultiYieldSurface;

class PIMYImplex : public NDMaterial
{
public:
    // Initialization constructor
    PIMYImplex(int tag,
        double rho,
        double refShearModul,
        double refBulkModul,
        double cohesi,
        double peakShearStra,
        double frictionAng,
        double refPress,
        double pressDependCoe,
        int   numberOfYieldSurf,
        double* gredu,
        int implexFlag,
        int verboseFlag);

    // Default constructor
    PIMYImplex();

    // Copy constructor
    PIMYImplex(const PIMYImplex&);

    // Destructor: clean up memory storage space.
    virtual ~PIMYImplex();

    const char* getClassType(void) const { return "PIMYImplex"; };

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
    const Vector& getImplexSateVariables(void);
    const Vector& getDissipatedEnergy(void);
    const Vector& getCommittedStrain(void);

    // Accepts the current trial strain values as being on the solution path, and updates 
    // all model parameters related to stress/strain states. Return 0 on success.
    int commitState(void);

    // Revert the stress/strain states to the last committed states. Return 0 on success.
    int revertToLastCommit(void);

    int revertToStart(void) { return 0; }

    // Return an exact copy of itself.
    NDMaterial* getCopy(void);

    // Return a copy of itself if "code"="PIMYImplex", otherwise return null.
    NDMaterial* getCopy(const char* code);

    // Return the string "PIMYImplex".
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

    static int matCount;
    static int* loadStagex;  //=0 if elastic; =1 if plastic

    // user supplied
    static int* ndmx;  //num of dimensions (2 or 3)
    static double* rhox;
    static double* frictionAnglex;
    static double* peakShearStrainx;
    static double* refPressurex;
    static double* cohesionx;
    static double* pressDependCoeffx;
    static int* numOfSurfacesx;
    static bool* doImplex;                 // solve implicit by default
    static bool* beVerbose;                // give internal solution info

    // internal
    static double* residualPressx;
    static Matrix theTangent;  //classwise member
    int e2p;
    int matN;
    double refShearModulus;
    double refBulkModulus;
    MultiYieldSurface* theSurfaces; // NOTE: surfaces[0] is not used  
    MultiYieldSurface* committedSurfaces;
    int    activeSurfaceNum;
    int    committedActiveSurf;
    T2Vector currentStress;
    T2Vector trialStress;
    T2Vector currentStrain;
    T2Vector strainRate;
    static T2Vector subStrainRate;
    double* mGredu;

    // implex state variables
    double lambda = 0.0;                     // step n+1 plastic multiplier
    double lambda_commit = 0.0;              // step n   plastic multiplier 
    double lambda_commit_old = 0.0;          // step n-1 plastic multiplier
    double dtime_n = 0.0;                    // time factor
    double dtime_n_commit = 0.0;             // committed time factor
    bool dtime_is_user_defined = false;
    bool dtime_first_set = false;

    // new results
    T2Vector currentStressImplex;            // extrapolated stress as a material result
    double lambda_bar = 0.0;
    double plasticStressNorm = 0.0;
    double plasticMultiplier = 0.0;
    double spentEnergy = 0.0;

    // move stress computation here
    int implicitSress(void);

    void elast2Plast(void);
    // Called by constructor
    void setUpSurfaces(double*);

    double yieldFunc(const T2Vector& stress, const MultiYieldSurface* surfaces,
        int surface_num);

    void deviatorScaling(T2Vector& stress, const MultiYieldSurface* surfaces,
        int surfaceNum, int count = 0);

    void initSurfaceUpdate(void);

    void paramScaling(void);

    // Return num_strain_subincre
    int setSubStrainRate(void);

    int isLoadReversal(void);

    void getContactStress(T2Vector& contactStress);

    void getSurfaceNormal(const T2Vector& stress, Vector& surfaceNormal);

    void setTrialStress(T2Vector& stress);

    double getLoadingFunc(const T2Vector& contact,
        const Vector& surfaceNormal, int crossedSurface);

    void stressCorrection(int crossedSurface);

    void updateActiveSurface(void);

    void updateInnerSurface(void);

    // Return 1 if crossing the active surface; return 0 o/w
    int  isCrossingNextSurface(void);

};

#endif