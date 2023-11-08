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

#include "SRSMYSand.h"
#include <MaterialResponse.h>

using tc = CTensor::Constants;

namespace tools {
    // keep track of all material instances local to the current process
    class MaterialList {
    public:
        size_t count = 0;
        int next_tag = 0;
        std::vector<std::unique_ptr<SRSMYSand>> mptr;

    public:
        MaterialList(void) = default;
        ~MaterialList(void) = default;

        void append(std::unique_ptr<SRSMYSand> matptr) {
            matptr->setSubTag(next_tag);
            next_tag++;
            mptr.push_back(std::move(matptr));
            count = mptr.size();
        }

        void remove(std::unique_ptr<SRSMYSand> matptr) {
            auto it = std::find(mptr.begin(), mptr.end(), matptr);
            if (it != mptr.end()) {
                mptr.erase(it);
            }
            // clean up, if the list is empty or, if not, update size
            if (mptr.empty()) { cleanup(); }
            else { count = mptr.size(); }
        }

        void cleanup(void) {
            count = 0;
            next_tag = 0;
            mptr.clear();
        }
    };

    // allocate static global storage for all instances
    class GlobalStorage {
    public:
        int size = 0;
        Matrix Ct;	//stiffness
        Matrix Ce;	//elasic stiffness
        Vector s;	//stress
        Vector e;	//strain


    public:
        GlobalStorage() = default;
        GlobalStorage& resize(int N) {
            if (N != size) {
                Ct.resize(N, N);
                Ce.resize(N, N);
                s.resize(N);
                e.resize(N);
            }
            return *this;
        }
    };


    static GlobalStorage& getGlobalStorage(int N) {
        static std::map<int, GlobalStorage> gsmap;
        return gsmap[N].resize(N);
    }
}


int SRSMYSand::matCount = 0;
int* SRSMYSand::loadStagex = 0;  //=0 if elastic; =1 if plastic
int* SRSMYSand::ndmx = 0;  //num of dimensions (2 or 3)
double* SRSMYSand::rhox = 0;
double* SRSMYSand::refShearModulusx = 0;
double* SRSMYSand::refBulkModulusx = 0;
double* SRSMYSand::frictionAnglex = 0;
double* SRSMYSand::peakShearStrainx = 0;
double* SRSMYSand::refPressurex = 0;
double* SRSMYSand::cohesionx = 0;
double* SRSMYSand::pressDependCoeffx = 0;
int* SRSMYSand::numOfSurfacesx = 0;
double* SRSMYSand::phaseTransfAnglex = 0;
double* SRSMYSand::contractParam1x = 0;
double* SRSMYSand::contractParam2x = 0;
double* SRSMYSand::contractParam3x = 0;
double* SRSMYSand::dilateParam1x = 0;
double* SRSMYSand::dilateParam2x = 0;
double* SRSMYSand::liquefyParam1x = 0;
double* SRSMYSand::liquefyParam2x = 0;
double* SRSMYSand::dilateParam3x = 0;
double* SRSMYSand::einitx = 0;    //initial void ratio
double* SRSMYSand::volLimit1x = 0;
double* SRSMYSand::volLimit2x = 0;
double* SRSMYSand::volLimit3x = 0;
double* SRSMYSand::residualPressx = 0;
double* SRSMYSand::stressRatioPTx = 0;
double* SRSMYSand::Hvx = 0;
double* SRSMYSand::Pvx = 0;
double SRSMYSand::pAtm = 101.;
const	double pi = 3.14159265358979;

// Public methods
    // full constructor
SRSMYSand::SRSMYSand(int tag, int nd,
    double r, double refShearModul,
    double refBulkModul, double frictionAng,
    double peakShearStra, double refPress,
    double pressDependCoe,
    double phaseTransformAng,
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
    double ei,
    double volLim1, double volLim2, double volLim3,
    double atm, double cohesi,
    double hv, double pv)
    :NDMaterial(tag, ND_TAG_SRSMYSand),
    trialStress(6, 2), updatedTrialStress(6, 2), strainRate(6, 1),
    PPZPivot(6, 2), PPZCenter(6, 2), PPZPivotCommitted(6, 2), PPZCenterCommitted(6, 2),
    PivotStrainRate(6, 1), PivotStrainRateCommitted(6, 1), check(0)
{
    opserr << "SRSMYSand::SRSMYSand() - in!\n";
    // handle dimension
    nD = OPS_GetNDM();
    if (nD != 2 && nD != 3) {
        opserr << "FATAL! SRSMYSand::SRSMYSand() - erroneous dimension recieved: needs to either 2 or 3" << "\n";
        exit(-1);
    }
    nDOF = (nD == 2) * (3) + (nD == 3) * (6);
    // do some checks
    if (refShearModul <= 0) {
        opserr << "FATAL:SRSMYSand:: refShearModulus <= 0" << endln;
        exit(-1);
    }
    if (refBulkModul <= 0) {
        opserr << "FATAL:SRSMYSand:: refBulkModulus <= 0" << endln;
        exit(-1);
    }
    if (frictionAng <= 0.) {
        opserr << "FATAL:SRSMYSand:: frictionAngle <= 0" << endln;
        exit(-1);
    }
    if (frictionAng >= 90.) {
        opserr << "FATAL:SRSMYSand:: frictionAngle >= 90" << endln;
        exit(-1);
    }
    if (phaseTransformAng <= 0.) {
        opserr << "FATAL:SRSMYSand:: phaseTransformAng " << phaseTransformAng << "<= 0" << endln;
        exit(-1);
    }
    if (cohesi < 0) {
        opserr << "WARNING:SRSMYSand:: cohesion < 0" << endln;
        opserr << "Will reset cohesion to 0.3." << endln;
        cohesi = 0.3;
    }
    if (peakShearStra <= 0) {
        opserr << "FATAL:SRSMYSand:: peakShearStra <= 0" << endln;
        exit(-1);
    }
    if (refPress <= 0) {
        opserr << "FATAL:SRSMYSand:: refPress <= 0" << endln;
        exit(-1);
    }
    if (pressDependCoe < 0) {
        opserr << "WARNING:SRSMYSand:: pressDependCoe < 0" << endln;
        opserr << "Will reset pressDependCoe to 0.5." << endln;
        pressDependCoe = 0.5;
    }
    if (numberOfYieldSurf <= 0) {
        opserr << "WARNING:SRSMYSand:: numberOfSurfaces " << numberOfYieldSurf << "<= 0" << endln;
        opserr << "Will use 10 yield surfaces." << endln;
        numberOfYieldSurf = 10;
    }
    if (numberOfYieldSurf > 100) {
        opserr << "WARNING:SRSMYSand::SRSMYSand: numberOfSurfaces > 100" << endln;
        opserr << "Will use 100 yield surfaces." << endln;
        numberOfYieldSurf = 100;
    }
    if (volLim1 < 0) {
        opserr << "WARNING:SRSMYSand:: volLim1 < 0" << endln;
        opserr << "Will reset volLimit to 0.8" << endln;
        volLim1 = 0.8;
    }
    if (r < 0) {
        opserr << "FATAL:SRSMYSand:: rho <= 0" << endln;
        exit(-1);
    }
    if (ei < 0) {
        opserr << "FATAL:SRSMYSand:: e <= 0" << endln;
        exit(-1);
    }
 
    if (matCount % 20 == 0) {
        int* temp1 = loadStagex;
        int* temp2 = ndmx;
        double* temp3 = rhox;
        double* temp4 = refShearModulusx;
        double* temp5 = refBulkModulusx;
        double* temp6 = frictionAnglex;
        double* temp7 = peakShearStrainx;
        double* temp8 = refPressurex;
        double* temp9 = cohesionx;
        double* temp10 = pressDependCoeffx;
        int* temp11 = numOfSurfacesx;
        double* temp12 = residualPressx;
        double* temp13 = phaseTransfAnglex;
        double* temp14 = contractParam1x;
        double* temp14a = contractParam2x;
        double* temp14b = contractParam3x;
        double* temp15 = dilateParam1x;
        double* temp16 = dilateParam2x;
        double* temp17 = liquefyParam1x;
        double* temp18 = liquefyParam2x;
        double* temp19 = dilateParam3x;
        double* temp20 = einitx;    //initial void ratio
        double* temp21 = volLimit1x;
        double* temp22 = volLimit2x;
        double* temp23 = volLimit3x;
        double* temp24 = stressRatioPTx;
        double* temp25 = Hvx;
        double* temp26 = Pvx;

        loadStagex = new int[matCount + 20];
        ndmx = new int[matCount + 20];
        rhox = new double[matCount + 20];
        refShearModulusx = new double[matCount + 20];
        refBulkModulusx = new double[matCount + 20];
        frictionAnglex = new double[matCount + 20];
        peakShearStrainx = new double[matCount + 20];
        refPressurex = new double[matCount + 20];
        cohesionx = new double[matCount + 20];
        pressDependCoeffx = new double[matCount + 20];
        numOfSurfacesx = new int[matCount + 20];
        residualPressx = new double[matCount + 20];
        phaseTransfAnglex = new double[matCount + 20];
        contractParam1x = new double[matCount + 20];
        contractParam2x = new double[matCount + 20];
        contractParam3x = new double[matCount + 20];
        dilateParam1x = new double[matCount + 20];
        dilateParam2x = new double[matCount + 20];
        liquefyParam1x = new double[matCount + 20];
        liquefyParam2x = new double[matCount + 20];
        dilateParam3x = new double[matCount + 20];
        einitx = new double[matCount + 20];    //initial void ratio
        volLimit1x = new double[matCount + 20];
        volLimit2x = new double[matCount + 20];
        volLimit3x = new double[matCount + 20];
        stressRatioPTx = new double[matCount + 20];
        Hvx = new double[matCount + 20];
        Pvx = new double[matCount + 20];

        for (int i = 0; i < matCount; i++) {
            loadStagex[i] = temp1[i];
            ndmx[i] = temp2[i];
            rhox[i] = temp3[i];
            refShearModulusx[i] = temp4[i];
            refBulkModulusx[i] = temp5[i];
            frictionAnglex[i] = temp6[i];
            peakShearStrainx[i] = temp7[i];
            refPressurex[i] = temp8[i];
            cohesionx[i] = temp9[i];
            pressDependCoeffx[i] = temp10[i];
            numOfSurfacesx[i] = temp11[i];
            residualPressx[i] = temp12[i];
            phaseTransfAnglex[i] = temp13[i];
            contractParam1x[i] = temp14[i];
            contractParam2x[i] = temp14a[i];
            contractParam3x[i] = temp14b[i];
            dilateParam1x[i] = temp15[i];
            dilateParam2x[i] = temp16[i];
            liquefyParam1x[i] = temp17[i];
            liquefyParam2x[i] = temp18[i];
            dilateParam3x[i] = temp19[i];
            einitx[i] = temp20[i];    //initial void ratio
            volLimit1x[i] = temp21[i];
            volLimit2x[i] = temp22[i];
            volLimit3x[i] = temp23[i];
            stressRatioPTx[i] = temp24[i];
            Hvx[i] = temp25[i];
            Pvx[i] = temp26[i];
        }

        if (matCount > 0) {
            delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
            delete[] temp5; delete[] temp6; delete[] temp7; delete[] temp8;
            delete[] temp9; delete[] temp10; delete[] temp11; delete[] temp12;
            delete[] temp13; delete[] temp14; delete[] temp14a; delete[] temp14b;
            delete[] temp15; delete[] temp16;
            delete[] temp17; delete[] temp18; delete[] temp19; delete[] temp20;
            delete[] temp21; delete[] temp22; delete[] temp23; delete[] temp24;
            delete[] temp25; delete[] temp26;
        }
    }
    
    ndmx[matCount] = nD;
    loadStagex[matCount] = 0;   //default
    refShearModulusx[matCount] = refShearModul;
    refBulkModulusx[matCount] = refBulkModul;
    frictionAnglex[matCount] = frictionAng;
    peakShearStrainx[matCount] = peakShearStra;
    refPressurex[matCount] = -refPress;  //compression is negative
    cohesionx[matCount] = cohesi;
    pressDependCoeffx[matCount] = pressDependCoe;
    numOfSurfacesx[matCount] = numberOfYieldSurf;
    rhox[matCount] = r;
    phaseTransfAnglex[matCount] = phaseTransformAng;
    contractParam1x[matCount] = contractionParam1;
    contractParam2x[matCount] = contractionParam2;
    contractParam3x[matCount] = contractionParam3;
    dilateParam1x[matCount] = dilationParam1;
    dilateParam2x[matCount] = dilationParam2;
    volLimit1x[matCount] = volLim1;
    volLimit2x[matCount] = volLim2;
    volLimit3x[matCount] = volLim3;
    liquefyParam1x[matCount] = liquefactionParam1;
    liquefyParam2x[matCount] = liquefactionParam2;
    dilateParam3x[matCount] = dilationParam3;
    einitx[matCount] = ei;
    Hvx[matCount] = hv;
    Pvx[matCount] = pv;

    residualPressx[matCount] = 0.;
    stressRatioPTx[matCount] = 0.;

    matN = matCount;
    matCount++;
    pAtm = atm;

    size_t numOfSurfaces = numOfSurfacesx[matN];
    initPress = refPressurex[matN];

    e2p = sv.nYs_commit = sv.nYs = 0;
    onPPZCommitted = onPPZ = -1;
    PPZSizeCommitted = PPZSize = 0.;
    pressureDCommitted = pressureD = modulusFactor = 0.;
    cumuDilateStrainOctaCommitted = cumuDilateStrainOcta = 0.;
    maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta = 0.;
    cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta = 0.;
    prePPZStrainOctaCommitted = prePPZStrainOcta = 0.;
    oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta = 0.;
    maxPress = 0.;
    damage = 0.;
    theSurfaces.clear();
    committedSurfaces.clear();
    theSurfaces.resize(numOfSurfaces + 1);
    committedSurfaces.resize(numOfSurfaces + 1);
    mGredu = gredu;
    setUpSurfaces(gredu);  // residualPress and stressRatioPT are calculated inside.

    // reset internal parameters
    revertToStart();
    opserr << "SRSMYSand::SRSMYSand() - out!\n";
}


SRSMYSand::SRSMYSand()
    :NDMaterial(0, ND_TAG_SRSMYSand), trialStress(),strainRate(), PPZPivot(), 
    PPZCenter(), PivotStrainRate(6, 1), PivotStrainRateCommitted(6, 1),
    PPZPivotCommitted(), PPZCenterCommitted(), theSurfaces(0), committedSurfaces(0)
{
    //does nothing
}


SRSMYSand::SRSMYSand(const SRSMYSand& other)
    :NDMaterial(other.getTag(), ND_TAG_SRSMYSand), trialStress(other.trialStress),
    strainRate(other.strainRate), check(0),PPZPivot(other.PPZPivot), PPZCenter(other.PPZCenter), 
    updatedTrialStress(other.updatedTrialStress), PPZPivotCommitted(other.PPZPivotCommitted), 
    PPZCenterCommitted(other.PPZCenterCommitted), PivotStrainRate(other.PivotStrainRate), 
    PivotStrainRateCommitted(other.PivotStrainRateCommitted)
{
    sv.sig = other.sv.sig;
    sv.eps = other.sv.eps;



    matN = other.matN;

    int numOfSurfaces = numOfSurfacesx[matN];

    e2p = other.e2p;
    strainPTOcta = other.strainPTOcta;
    modulusFactor = other.modulusFactor;
    sv.nYs = other.sv.nYs;
    sv.nYs_commit = other.sv.nYs_commit;
    pressureDCommitted = other.pressureDCommitted;
    onPPZCommitted = other.onPPZCommitted;
    PPZSizeCommitted = other.PPZSizeCommitted;
    cumuDilateStrainOctaCommitted = other.cumuDilateStrainOctaCommitted;
    maxCumuDilateStrainOctaCommitted = other.maxCumuDilateStrainOctaCommitted;
    cumuTranslateStrainOctaCommitted = other.cumuTranslateStrainOctaCommitted;
    prePPZStrainOctaCommitted = other.prePPZStrainOctaCommitted;
    oppoPrePPZStrainOctaCommitted = other.oppoPrePPZStrainOctaCommitted;
    pressureD = other.pressureD;
    onPPZ = other.onPPZ;
    PPZSize = other.PPZSize;
    cumuDilateStrainOcta = other.cumuDilateStrainOcta;
    maxCumuDilateStrainOcta = other.maxCumuDilateStrainOcta;
    cumuTranslateStrainOcta = other.cumuTranslateStrainOcta;
    prePPZStrainOcta = other.prePPZStrainOcta;
    oppoPrePPZStrainOcta = other.oppoPrePPZStrainOcta;
    initPress = other.initPress;
    maxPress = other.maxPress;
    damage = other.damage;
    theSurfaces.clear();
    committedSurfaces.clear();
    theSurfaces.resize(numOfSurfaces + 1);
    committedSurfaces.resize(numOfSurfaces + 1);
    for (int i = 1; i <= numOfSurfaces; i++) {
        committedSurfaces[i] = other.committedSurfaces[i];
        theSurfaces[i] = other.theSurfaces[i];
    }
}


SRSMYSand::~SRSMYSand()
{
    opserr << "SRSMYSand::~SRSMYSand()\n";
}


void SRSMYSand::elast2Plast(void)
{
    opserr << "SRSMYSand::elast2Plast() - in!\n";
    int loadStage = loadStagex[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (loadStage != 1 || e2p == 1)
        return;

    e2p = 1;

    if (sv.sig.trace() > 0.0) {
        if (beVerbose) { opserr << "WARNING! SRSMYSand::elast2Plast(): material in tension\n"; }
        sv.sig = sv.sig.deviator();  // assume zero mean pressure
    }
    
    // Active surface is 0, return
    sv.sig = sv.sig.deviator();
    if (sv.sig.norm() == 0.) return;
    
    // Find active surface
    while (yieldFunc(sv.sig, committedSurfaces, ++sv.nYs_commit) > 0) {
        if (sv.nYs_commit == numOfSurfaces) {
            if (beVerbose) { opserr << "WARNING:SRSMYSand::elast2Plast(): stress out of failure surface\n"; }
            deviatorScaling(sv.sig, committedSurfaces, numOfSurfaces);
            initSurfaceUpdate();
            return;
        }
    }

    sv.nYs_commit--;
    initSurfaceUpdate();
    opserr << "SRSMYSand::elast2Plast() - out!\n";
}


int SRSMYSand::setTrialStrain(const Vector& strain)
{
    sv.eps = CTensor(strain, 1); // 1: covariant
    if (beVerbose)  opserr << "SRSMYSand::setTrialStrain() - " << sv.eps << "\n";
    // implex time step
    if (!sv.dtime_is_user_defined) {
        sv.dtime_n = ops_Dt;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
    }
    // trigger material update
    if (useImplex) {
        if (beVerbose) {
            opserr << "SRSMYSand::setTrialStrain() -nD Material "
                << getTag() << "::" << getSubTag() << " -> IMPL-EX: explicit stage...\n";
        }
        updateInternal(true, true);
        sv.sig_implex = sv.sig; // save stress for output
    }
    else {
        if (beVerbose) {
            opserr << "SRSMYSand::setTrialStrain() - nD Material "
                << getTag() << "::" << getSubTag() << " -> IMPLICIT...\n";
        }
        updateInternal(false, true);
    }
    opserr << "SRSMYSand::setTrialStrain() - out!\n";
    return 0;
}


int SRSMYSand::setTrialStrain(const Vector& strain, const Vector& rate)
{
    return setTrialStrain(strain);
}


int SRSMYSand::setTrialStrainIncr(const Vector& strain)
{
    sv.eps += CTensor(strain, 1); // 1: covariant
    if (beVerbose)  opserr << "SRSMYSand::setTrialStrainIncr() - " << sv.eps << "\n";
    // implex time step
    if (!sv.dtime_is_user_defined) {
        sv.dtime_n = ops_Dt;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
    }
    // trigger material update
    if (useImplex) {
        if (beVerbose) {
            opserr << "SRSMYSand::setTrialStrainIncr() -nD Material "
                << getTag() << "::" << getSubTag() << " -> IMPL-EX: explicit stage...\n";
        }
        updateInternal(true, true);
        sv.sig_implex = sv.sig; // save stress for output
    }
    else {
        if (beVerbose) {
            opserr << "SRSMYSand::setTrialStrainIncr() - nD Material "
                << getTag() << "::" << getSubTag() << " -> IMPLICIT...\n";
        }
        updateInternal(false, true);
    }
    opserr << "SRSMYSand::setTrialStrainIncr() - out!\n";
    return 0;
}


int SRSMYSand::setTrialStrainIncr(const Vector& strain, const Vector& rate)
{
    return setTrialStrainIncr(strain);
}


const Matrix& SRSMYSand::getTangent(void)
{
    opserr << "SRSMYSand::getTangent() - in!\n";
    auto& gs = tools::getGlobalStorage(nDOF);
    Matrix& stiff = gs.Ct;
    CTensor normal(nDOF, 2); // 2: contravariant

    int loadStage = loadStagex[matN];
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refPressure = refPressurex[matN];
    double residualPress = residualPressx[matN];
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 3;

    if (loadStage == 1 && e2p == 0) {
        initPress = sv.sig.trace() / nD;
        elast2Plast();
    }
    if (loadStage == 2 && initPress == refPressure)
        initPress = sv.sig.trace() / nD;

    if (loadStage == 0 || loadStage == 2) {  //linear elastic
        double factor;
        if (loadStage == 0)
            factor = 1.0;
        else {
            factor = (initPress - residualPress) / (refPressure - residualPress);
            if (factor <= 1.e-10) factor = 1.e-10;
            else factor = pow(factor, pressDependCoeff);
            factor = (1.e-10 > factor) ? 1.e-10 : factor;
        }
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                stiff(i, j) = 0.;
                if (i == j)
                    stiff(i, j) += refShearModulus * factor;
                if (i < 3 && j < 3 && i == j)
                    stiff(i, j) += refShearModulus * factor;
                if (i < 3 && j < 3)
                    stiff(i, j) += (refBulkModulus - 2. * refShearModulus / 3.) * factor;
            }
    }
    else {
        double coeff1, coeff2, coeff3, coeff4;
        double factor = getModulusFactor(updatedTrialStress);
        double shearModulus = factor * refShearModulus;
        double bulkModulus = factor * refBulkModulus;

        // volumetric plasticity
        if (Hvx[matN] != 0. && trialStress.trace() / nD <= maxPress
            && strainRate.trace() / nD < 0. && loadStage == 1) {
            double tp = fabs(trialStress.trace() / nD - residualPress);
            bulkModulus = (bulkModulus * Hvx[matN] * pow(tp, Pvx[matN])) / (bulkModulus + Hvx[matN] * pow(tp, Pvx[matN]));
        }

        if (loadStage != 0 && sv.nYs > 0) {
            //	 opserr << "PDMY02::getTang() - 5\n";
            factor = getModulusFactor(trialStress);
            shearModulus = factor * refShearModulus;
            bulkModulus = factor * refBulkModulus;
            getSurfaceNormal(trialStress, normal);
            double volume = normal.trace() / nD;
            normal = normal.deviator();
            double Ho = 9. * bulkModulus * volume * volume + 2. * shearModulus * (normal % normal);
            double plastModul = factor * theSurfaces[sv.nYs].modulus();
            coeff1 = 9. * bulkModulus * bulkModulus * volume * volume / (Ho + plastModul);
            coeff2 = 4. * shearModulus * shearModulus / (Ho + plastModul);
        }

        else {
            coeff1 = coeff2 = coeff3 = coeff4 = 0.;
            normal.Zero();
        }

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                stiff(i, j) = -coeff2 * normal(i) * normal(j);
                if (i == j) stiff(i, j) += shearModulus;
                if (i < 3 && j < 3 && i == j) stiff(i, j) += shearModulus;
                if (i < 3 && j < 3) stiff(i, j) += (bulkModulus - 2. * shearModulus / 3. - coeff1);
            }
    }

    opserr << "SRSMYSand::getTangent() - out!\n";
    return stiff;
}


const Matrix& SRSMYSand::getInitialTangent(void)
{
    opserr << "SRSMYSand::getInitialTangent() - in!\n";
    auto& gs = tools::getGlobalStorage(nDOF);
    Matrix& stiff = gs.Ce;

    int loadStage = loadStagex[matN];
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refPressure = refPressurex[matN];
    double residualPress = residualPressx[matN];
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 3;

    if (loadStage == 1 && e2p == 0) {
        initPress = sv.sig.trace() / nD;
        elast2Plast();
    }
    if (loadStage == 2 && initPress == refPressure)
        initPress = sv.sig.trace() / nD;
    double factor;
    if (loadStage == 0)
        factor = 1.;
    else if (loadStage == 2) {
        factor = (initPress - residualPress) / (refPressure - residualPress);
        if (factor <= 1.e-10) factor = 1.e-10;
        else factor = pow(factor, pressDependCoeff);
        factor = (1.e-10 > factor) ? 1.e-10 : factor;
    }
    else if (loadStage == 1)
        factor = getModulusFactor(sv.sig);

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            stiff(i, j) = 0.;
            if (i == j) stiff(i, j) += refShearModulus * factor;
            if (i < 3 && j < 3 && i == j) stiff(i, j) += refShearModulus * factor;
            if (i < 3 && j < 3) stiff(i, j) += (refBulkModulus - 2. * refShearModulus / 3.) * factor;
        }

    opserr << "SRSMYSand::getInitialTangent() - out!\n";
    return stiff;
    /* non-symmetric stiffness
        workM(0,2) = theTangent(0,3);
        workM(1,2) = theTangent(1,3);
        workM(2,0) = theTangent(3,0);
        workM(2,1) = theTangent(3,1);*/
}


const Vector& SRSMYSand::getStress(void)
{
    auto& gs = tools::getGlobalStorage(nDOF);
    Vector& stress = gs.s;
    if (beVerbose)  opserr << "SRSMYSand::getStress() - " << sv.sig << "\n";
    stress = sv.sig.makeVector();
    return stress;
}


const Vector& SRSMYSand::getStrain(void)
{
    auto& gs = tools::getGlobalStorage(nDOF);
    Vector& strain = gs.e;
    if (beVerbose)  opserr << "SRSMYSand::getStrain() - " << sv.eps << "\n";
    strain = sv.eps.makeVector();
    return strain;
}


int SRSMYSand::commitState(void)
{
    opserr << "SRSMYSand::commitState() - in!\n";

    // do the implicit correction if impl-ex
    if (useImplex) {
        // update material internal variables
        if (beVerbose) {
            opserr << "SRSMYSand::commitState() - nD Material "
                << getTag() << "::" << getSubTag() << " -> IMPL-EX: implicit correction...\n";
        }
        updateInternal(false, false);  // explicit_phase?, do_tangent?
    }

    // commit internal variables
    if (sv.commitState()) {
        opserr << "FATAL: SRSMYSand::commitState() - material internal variables failed commit!\n";
        exit(-1);
    }

    int loadStage = loadStagex[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (loadStage == 1) {
        sv.nYs_commit = sv.nYs;
        for (int i = 1; i <= numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
        pressureDCommitted = pressureD;
        onPPZCommitted = onPPZ;
        PPZSizeCommitted = PPZSize;
        cumuDilateStrainOctaCommitted = cumuDilateStrainOcta;
        maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta;
        cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta;
        prePPZStrainOctaCommitted = prePPZStrainOcta;
        oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta;
        PPZPivotCommitted = PPZPivot;
        PivotStrainRateCommitted = PivotStrainRate;
        PPZCenterCommitted = PPZCenter;
        if (sv.sig.trace() / nD < maxPress) {
            maxPressCommitted = maxPress;
            maxPress = sv.sig.trace() / nD;
        }
    }

    if (beVerbose) {
        opserr << "SRSMYSand::commitState() - nD Material " << getTag() << "::" << getSubTag() << " ->\n";
        sv.printStats(false);
    }

    // done
    opserr << "SRSMYSand::commitState() - out!\n";
    return 0;
}


int SRSMYSand::revertToLastCommit(void)
{
    // restore committed internal variables
    if (sv.revertToLastCommit()) {
        opserr << "FATAL: SRSMYSand::revertToLastCommit() - material internal variables failed reverting to last commit!\n";
        exit(-1);
    }

    int loadStage = loadStagex[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (loadStage == 1) {
        sv.nYs = sv.nYs_commit;
        for (int i = 1; i <= numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
        pressureD = pressureDCommitted;
        onPPZ = onPPZCommitted;
        PPZSize = PPZSizeCommitted;
        cumuDilateStrainOcta = cumuDilateStrainOctaCommitted;
        maxCumuDilateStrainOcta = maxCumuDilateStrainOctaCommitted;
        cumuTranslateStrainOcta = cumuTranslateStrainOctaCommitted;
        prePPZStrainOcta = prePPZStrainOctaCommitted;
        oppoPrePPZStrainOcta = oppoPrePPZStrainOctaCommitted;
        PPZPivot = PPZPivotCommitted;
        PivotStrainRate = PivotStrainRateCommitted;
        PPZCenter = PPZCenterCommitted;
        maxPress = maxPressCommitted;
    }

    if (beVerbose) { opserr << "SRSMYSand::revertToLastCommit() - nD Material " << getTag() << "::" << getSubTag() << "\n"; }

    // done
    return 0;
}

int SRSMYSand::revertToStart(void) {
    // reset state variables
    sv = MaterialStateVariables(nD);
    loadStagex[matN] = 0;
    // done
    return 0;
}

NDMaterial* SRSMYSand::getCopy(void)
{
    opserr << "SRSMYSand::getCopy() - in!\n";
    SRSMYSand* copy = new SRSMYSand(*this);
    opserr << "SRSMYSand::getCopy() - out!\n";
    return copy;
}


NDMaterial* SRSMYSand::getCopy(const char* code)
{
    opserr << "SRSMYSand::getCopy(char*) - in!\n";
    if (strcmp(code, "PlaneStrain") == 0 || strcmp(code, "ThreeDimensional") == 0) {
        SRSMYSand* copy = new SRSMYSand(*this);
        return copy;
    }
    else {
        opserr << "ERROR SRSMYSand::getCopy -- cannot make copy for type " << code << endln;
        return 0;
    }
    opserr << "SRSMYSand::getCopy(char*) - out!\n";
}


const char* SRSMYSand::getType(void) const
{
    opserr << "SRSMYSand::getType() - in!\n";
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;

    opserr << "SRSMYSand::getType() - out!\n";
    return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int SRSMYSand::getOrder(void) const
{
    opserr << "SRSMYSand::getOrder() - in!\n";
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;

    return (ndm == 2) ? 3 : 6;
    opserr << "SRSMYSand::getOrder() - out!\n";
}

int SRSMYSand::setParameter(const char** argv, int argc, Parameter& param)
{
    opserr << "SRSMYSand::setParameter() - in!\n";
    /*if (argc < 1)
      return -1;

    if (strcmp(argv[0],"updateMaterialStage") == 0) {
      if (argc < 2)
        return -1;
      int matTag = atoi(argv[1]);
      if (this->getTag() == matTag)
        return param.addObject(1, this);
      else
        return -1;
    }

    else if (strcmp(argv[0],"shearModulus") == 0)
      return param.addObject(10, this);

    else if (strcmp(argv[0],"bulkModulus") == 0)
      return param.addObject(11, this);

    return -1;*/

    if (argc < 2)
        return -1;

    int theMaterialTag;
    theMaterialTag = atoi(argv[1]);

    opserr << "SRSMYSand::setParameter() - out!\n";

    // check for material tag
    if (theMaterialTag == this->getTag()) {

        if (strcmp(argv[0], "updateMaterialStage") == 0) {
            return param.addObject(1, this);
        }
        else if (strcmp(argv[0], "shearModulus") == 0) {
            return param.addObject(10, this);
        }
        else if (strcmp(argv[0], "bulkModulus") == 0) {
            return param.addObject(11, this);
        }
        else if (strcmp(argv[0], "frictionAngle") == 0) {
            return param.addObject(12, this);
        }
        else if (strcmp(argv[0], "cohesion") == 0) {
            return param.addObject(13, this);
        }
    }
    return -1;
}

int SRSMYSand::updateParameter(int responseID, Information& info)
{
    opserr << "SRSMYSand::updateParameter() - in!\n";
    if (responseID == 1) {
        loadStagex[matN] = info.theInt;
    }
    else if (responseID == 10) {
        refShearModulusx[matN] = info.theDouble;
    }
    else if (responseID == 11) {
        refBulkModulusx[matN] = info.theDouble;
    }
    else if (responseID == 12) {
        frictionAnglex[matN] = info.theDouble;
        setUpSurfaces(mGredu);
        initSurfaceUpdate();
    }
    else if (responseID == 13) {
        cohesionx[matN] = info.theDouble;
        setUpSurfaces(mGredu);
        initSurfaceUpdate();
    }

    // used by BBarFourNodeQuadUP element
    else if (responseID == 20 && ndmx[matN] == 2)
        ndmx[matN] = 0;

    opserr << "SRSMYSand::updateParameter() - out!\n";
    return 0;
}


int SRSMYSand::sendSelf(int commitTag, Channel& theChannel)
{
    opserr << "SRSMYSand::sendSelf() - in!\n";
    // ndmx[matCount] = nd;
    // loadStagex[matCount] = 0;   //default
     //refShearModulusx[matCount] = refShearModul;
     //refBulkModulusx[matCount] = refBulkModul;
     //frictionAnglex[matCount] = frictionAng;
     //peakShearStrainx[matCount] = peakShearStra;
     //refPressurex[matCount] = -refPress;  //compression is negative
     //cohesionx[matCount] = cohesi;
     //pressDependCoeffx[matCount] = pressDependCoe;
     //numOfSurfacesx[matCount] = numberOfYieldSurf;
     // rhox[matCount] = r;
     //phaseTransfAnglex[matCount] = phaseTransformAng;
     //contractParam1x[matCount] = contractionParam1;
     //contractParam2x[matCount] = contractionParam2;

     //dilateParam1x[matCount] = dilationParam1;
     //dilateParam2x[matCount] = dilationParam2;
     //volLimit1x[matCount] = volLim1;
     //volLimit2x[matCount] = volLim2;
     //volLimit3x[matCount] = volLim3;
     //liquefyParam1x[matCount] = liquefactionParam1;
     //liquefyParam2x[matCount] = liquefactionParam2;
     //dilateParam3x[matCount] = dilationParam3;
     //einitx[matCount] = ei;

       /*
     contractParam3x[matCount] = contractionParam3;
     Hvx[matCount] = hv;
     Pvx[matCount] = pv;
     */

    int loadStage = loadStagex[matN];
    int ndm = ndmx[matN];
    double rho = rhox[matN];
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];
    double frictionAngle = frictionAnglex[matN];
    double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];
    double phaseTransfAngle = phaseTransfAnglex[matN];
    double stressRatioPT = stressRatioPTx[matN];
    double contractParam1 = contractParam1x[matN];
    double contractParam2 = contractParam2x[matN];
    double dilateParam1 = dilateParam1x[matN];
    double dilateParam2 = dilateParam2x[matN];
    double liquefyParam1 = liquefyParam1x[matN];
    double liquefyParam2 = liquefyParam2x[matN];
    double dilateParam3 = dilateParam3x[matN];
    double einit = einitx[matN];
    double volLimit1 = volLimit1x[matN];
    double volLimit2 = volLimit2x[matN];
    double volLimit3 = volLimit3x[matN];

    double contractionParam3 = contractParam3x[matN];
    double hv = Hvx[matN];
    double Pv = Pvx[matN];

    int i, res = 0;

    static ID idData(6);
    idData(0) = this->getTag();
    idData(1) = numOfSurfaces;
    idData(2) = loadStage;
    idData(3) = ndm;
    idData(4) = matN;
    idData(5) = matCount;

    res += theChannel.sendID(this->getDbTag(), commitTag, idData);
    if (res < 0) {
        opserr << "SRSMYSand::sendSelf -- could not send ID\n";
        return res;
    }

    Vector data(69 + numOfSurfaces * 8);
    data(0) = rho;
    data(1) = einit;
    data(2) = refShearModulus;
    data(3) = refBulkModulus;
    data(4) = frictionAngle;
    data(5) = peakShearStrain;
    data(6) = refPressure;
    data(7) = cohesion;
    data(8) = pressDependCoeff;
    data(9) = phaseTransfAngle;
    data(10) = contractParam1;
    data(11) = dilateParam1;
    data(12) = dilateParam2;
    data(13) = volLimit1;
    data(14) = volLimit2;
    data(15) = volLimit3;
    data(16) = pAtm;
    data(17) = liquefyParam1;
    data(18) = liquefyParam2;
    data(19) = dilateParam3;
    data(20) = residualPress;
    data(21) = stressRatioPT;
    data(22) = e2p;
    data(23) = sv.nYs_commit;
    data(24) = strainPTOcta;
    data(25) = pressureDCommitted;
    data(26) = onPPZCommitted;
    data(27) = PPZSizeCommitted;
    data(28) = cumuDilateStrainOctaCommitted;
    data(29) = maxCumuDilateStrainOctaCommitted;
    data(30) = cumuTranslateStrainOctaCommitted;
    data(31) = prePPZStrainOctaCommitted;
    data(32) = oppoPrePPZStrainOctaCommitted;

    data(33) = initPress;
    data(34) = contractParam2;
    data(35) = contractionParam3;// = contractParam3x[matN];
    data(36) = hv; //  = Hvx[matN];
    data(37) = Pv; //  = Pvx[matN];

    for (i = 0; i < 6; i++) data(i + 38) = sv.sig(i);
    for (i = 0; i < 6; i++) data(i + 44) = sv.eps(i);
    for (i = 0; i < 6; i++) data(i + 50) = PPZPivotCommitted(i);
    for (i = 0; i < 6; i++) data(i + 56) = PPZCenterCommitted(i);

    for (i = 0; i < numOfSurfaces; i++) {
        int k = 62 + i * 8;
        data(k) = committedSurfaces[i + 1].size();
        data(k + 1) = committedSurfaces[i + 1].modulus();
        data(k + 2) = committedSurfaces[i + 1].center()(0);
        data(k + 3) = committedSurfaces[i + 1].center()(1);
        data(k + 4) = committedSurfaces[i + 1].center()(2);
        data(k + 5) = committedSurfaces[i + 1].center()(3);
        data(k + 6) = committedSurfaces[i + 1].center()(4);
        data(k + 7) = committedSurfaces[i + 1].center()(5);
    }

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "SRSMYSand::sendSelf -- could not send Vector\n";
        return res;
    }
    opserr << "SRSMYSand::sendSelf() - out!\n";
    return res;
}


int SRSMYSand::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    opserr << "SRSMYSand::recvSelf() - in!\n";
    int i, res = 0;

    static ID idData(6);
    res += theChannel.recvID(this->getDbTag(), commitTag, idData);
    if (res < 0) {
        opserr << "SRSMYSand::recvelf -- could not recv ID\n";

        return res;
    }

    this->setTag(idData(0));
    int numOfSurfaces = idData(1);
    int loadStage = idData(2);
    int ndm = idData(3);
    matN = idData(4);

    int otherMatCount = idData(5);

    Vector data(69 + idData(1) * 8);
    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "SRSMYSand::recvSelf -- could not recv Vector\n";
        return res;
    }

    double rho = data(0);
    double einit = data(1);
    double refShearModulus = data(2);
    double refBulkModulus = data(3);
    double frictionAngle = data(4);
    double peakShearStrain = data(5);
    double refPressure = data(6);
    double cohesion = data(7);
    double pressDependCoeff = data(8);
    double phaseTransfAngle = data(9);
    double contractParam1 = data(10);
    double dilateParam1 = data(11);
    double dilateParam2 = data(12);
    double volLimit1 = data(13);
    double volLimit2 = data(14);
    double volLimit3 = data(15);
    pAtm = data(16);
    double liquefyParam1 = data(17);
    double liquefyParam2 = data(18);
    double dilateParam3 = data(19);
    double residualPress = data(20);
    double stressRatioPT = data(21);
    e2p = data(22);
    sv.nYs_commit = data(23);
    strainPTOcta = data(24);
    pressureDCommitted = data(25);
    onPPZCommitted = data(26);
    PPZSizeCommitted = data(27);
    cumuDilateStrainOctaCommitted = data(28);
    maxCumuDilateStrainOctaCommitted = data(29);
    cumuTranslateStrainOctaCommitted = data(30);
    prePPZStrainOctaCommitted = data(31);
    oppoPrePPZStrainOctaCommitted = data(32);

    initPress = data(33);
    double contractParam2 = data(34);
    double contractParam3 = data(35); //  = contractionParam3;// = contractParam3x[matN];
    double hv = data(36); // =  hv; //  = Hvx[matN];
    double Pv = data(37); // = Pv; //  = Pvx[matN];


    for (i = 0; i < 6; i++) sv.sig(i) = data(i + 38);
    for (i = 0; i < 6; i++) sv.eps(i) = data(i + 44);
    for (i = 0; i < 6; i++) PPZPivotCommitted(i) = data(i + 50);
    for (i = 0; i < 6; i++) PPZCenterCommitted(i) = data(i + 56);

    theSurfaces.clear();
    committedSurfaces.clear();
    theSurfaces.resize(numOfSurfaces + 1);
    committedSurfaces.resize(numOfSurfaces + 1);

    CTensor center(6, 2);   // 2: contravariant
    for (i = 0; i < numOfSurfaces; i++) {
        int k = 62 + i * 8;
        center(0) = data(k + 2);
        center(1) = data(k + 3);
        center(2) = data(k + 4);
        center(3) = data(k + 5);
        center(4) = data(k + 6);
        center(5) = data(k + 7);
        committedSurfaces[i + 1].setData(center, data(k), data(k + 1));
    }

    int* temp1, * temp2, * temp11;
    double* temp3, * temp4, * temp5, * temp6, * temp7, * temp8, * temp9, * temp10, * temp12;
    double* temp13, * temp14, * temp15, * temp16, * temp17, * temp18, * temp19, * temp20;
    double* temp14a, * temp14b;
    double* temp21, * temp22, * temp23, * temp24, * temp25, * temp26;

    if (matCount < otherMatCount) {  // allocate memory if not enough

        temp1 = loadStagex;
        temp2 = ndmx;
        temp3 = rhox;
        temp4 = refShearModulusx;
        temp5 = refBulkModulusx;
        temp6 = frictionAnglex;
        temp7 = peakShearStrainx;
        temp8 = refPressurex;
        temp9 = cohesionx;
        temp10 = pressDependCoeffx;
        temp11 = numOfSurfacesx;
        temp12 = residualPressx;
        temp13 = phaseTransfAnglex;
        temp14 = contractParam1x;
        temp14a = contractParam2x;
        temp14b = contractParam3x;
        temp15 = dilateParam1x;
        temp16 = dilateParam2x;
        temp17 = liquefyParam1x;
        temp18 = liquefyParam2x;
        temp19 = dilateParam3x;
        temp20 = einitx;    //initial void ratio
        temp21 = volLimit1x;
        temp22 = volLimit2x;
        temp23 = volLimit3x;
        temp24 = stressRatioPTx;
        temp25 = Hvx;
        temp26 = Pvx;

        loadStagex = new int[otherMatCount];
        ndmx = new int[otherMatCount];
        rhox = new double[otherMatCount];
        refShearModulusx = new double[otherMatCount];
        refBulkModulusx = new double[otherMatCount];
        frictionAnglex = new double[otherMatCount];
        peakShearStrainx = new double[otherMatCount];
        refPressurex = new double[otherMatCount];
        cohesionx = new double[otherMatCount];
        pressDependCoeffx = new double[otherMatCount];
        numOfSurfacesx = new int[otherMatCount];
        residualPressx = new double[otherMatCount];
        phaseTransfAnglex = new double[otherMatCount];
        contractParam1x = new double[otherMatCount];
        contractParam2x = new double[otherMatCount];
        contractParam3x = new double[otherMatCount];
        dilateParam1x = new double[otherMatCount];
        dilateParam2x = new double[otherMatCount];
        liquefyParam1x = new double[otherMatCount];
        liquefyParam2x = new double[otherMatCount];
        dilateParam3x = new double[otherMatCount];
        einitx = new double[otherMatCount];    //initial void ratio
        volLimit1x = new double[otherMatCount];
        volLimit2x = new double[otherMatCount];
        volLimit3x = new double[otherMatCount];
        stressRatioPTx = new double[otherMatCount];
        Hvx = new double[otherMatCount];
        Pvx = new double[otherMatCount];

        if (matCount > 0) {
            for (int i = 0; i < matCount; i++) {
                loadStagex[i] = temp1[i];
                ndmx[i] = temp2[i];
                rhox[i] = temp3[i];
                refShearModulusx[i] = temp4[i];
                refBulkModulusx[i] = temp5[i];
                frictionAnglex[i] = temp6[i];
                peakShearStrainx[i] = temp7[i];
                refPressurex[i] = temp8[i];
                cohesionx[i] = temp9[i];
                pressDependCoeffx[i] = temp10[i];
                numOfSurfacesx[i] = temp11[i];
                residualPressx[i] = temp12[i];
                phaseTransfAnglex[i] = temp13[i];
                contractParam1x[i] = temp14[i];
                contractParam2x[i] = temp14a[i];
                contractParam3x[i] = temp14b[i];
                dilateParam1x[i] = temp15[i];
                dilateParam2x[i] = temp16[i];
                liquefyParam1x[i] = temp17[i];
                liquefyParam2x[i] = temp18[i];
                dilateParam3x[i] = temp19[i];
                einitx[i] = temp20[i];    //initial void ratio
                volLimit1x[i] = temp21[i];
                volLimit2x[i] = temp22[i];
                volLimit3x[i] = temp23[i];
                stressRatioPTx[i] = temp24[i];
                Hvx[i] = temp25[i];
                Pvx[i] = temp26[i];
            }

            delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
            delete[] temp5; delete[] temp6; delete[] temp7; delete[] temp8;
            delete[] temp9; delete[] temp10; delete[] temp11; delete[] temp12;
            delete[] temp13; delete[] temp14; delete[] temp15; delete[] temp16;
            delete[] temp14a;
            delete[] temp17; delete[] temp18; delete[] temp19; delete[] temp20;
            delete[] temp21; delete[] temp22; delete[] temp23; delete[] temp24;
            delete[] temp25; delete[] temp26;
        }
        matCount = otherMatCount;
    }

    loadStagex[matN] = loadStage;
    ndmx[matN] = ndm;
    rhox[matN] = rho;
    residualPressx[matN] = residualPress;
    numOfSurfacesx[matN] = numOfSurfaces;
    refPressurex[matN] = refPressure;
    pressDependCoeffx[matN] = pressDependCoeff;
    refShearModulusx[matN] = refShearModulus;
    refBulkModulusx[matN] = refBulkModulus;
    frictionAnglex[matN] = frictionAngle;
    cohesionx[matN] = cohesion;
    peakShearStrainx[matN] = peakShearStrain;
    phaseTransfAnglex[matN] = phaseTransfAngle;
    stressRatioPTx[matN] = stressRatioPT;
    contractParam1x[matN] = contractParam1;
    contractParam2x[matN] = contractParam2;
    contractParam3x[matN] = contractParam3;
    dilateParam1x[matN] = dilateParam1;
    dilateParam2x[matN] = dilateParam2;
    liquefyParam1x[matN] = liquefyParam1;
    liquefyParam2x[matN] = liquefyParam2;
    dilateParam3x[matN] = dilateParam3;
    einitx[matN] = einit;
    volLimit1x[matN] = volLimit1;
    volLimit2x[matN] = volLimit2;
    volLimit3x[matN] = volLimit3;

    opserr << "SRSMYSand::recvSelf() - out!\n";
    return res;
}


Response*
SRSMYSand::setResponse(const char** argv, int argc, OPS_Stream& s)
{
    opserr << "SRSMYSand::setResponse()\n";
    // begin change by Alborz Ghofrani - UW --- get only 6 components of stress
    if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
        if ((argc > 1) && (atoi(argv[1]) > 2) && (atoi(argv[1]) < 8))
            return new MaterialResponse(this, 2 + atoi(argv[1]), this->getStressToRecord(atoi(argv[1])));
        else
            return new MaterialResponse(this, 1, this->getCommittedStress());
    // end change by Alborz Ghofrani - UW

    else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
        return new MaterialResponse(this, 2, this->getCommittedStrain());

    else if (strcmp(argv[0], "tangent") == 0)
        return new MaterialResponse(this, 3, this->getTangent());

    else if (strcmp(argv[0], "backbone") == 0) {
        int numOfSurfaces = numOfSurfacesx[matN];
        Matrix curv(numOfSurfaces + 1, (argc - 1) * 2);
        for (int i = 1; i < argc; i++)
            curv(0, (i - 1) * 2) = atoi(argv[i]);
        return new MaterialResponse(this, 4, curv);
    }
    else
        return 0;
}


void SRSMYSand::getBackbone(Matrix& bb)
{
    opserr << "SRSMYSand::getBackbone() - in!\n";
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    double vol, conHeig, scale, factor, shearModulus, stress1,
        stress2, strain1, strain2, plastModulus, elast_plast, gre;

    for (int k = 0; k < bb.noCols() / 2; k++) {
        vol = bb(0, k * 2);
        if (vol <= 0.) {
            opserr << k << "\nNDMaterial " << this->getTag()
                << ": invalid confinement for backbone recorder, " << vol << endln;
            continue;
        }
        conHeig = vol + residualPress;
        scale = -conHeig / (refPressure - residualPress);
        factor = pow(scale, pressDependCoeff);
        shearModulus = factor * refShearModulus;

        for (int i = 1; i <= numOfSurfaces; i++) {
            if (i == 1) {
                stress2 = committedSurfaces[i].size() * conHeig / sqrt(3.0);
                strain2 = stress2 / shearModulus;
                bb(1, k * 2) = strain2; bb(1, k * 2 + 1) = shearModulus;
            }
            else {
                stress1 = stress2; strain1 = strain2;
                plastModulus = factor * committedSurfaces[i - 1].modulus();
                elast_plast = 2 * shearModulus * plastModulus / (2 * shearModulus + plastModulus);
                stress2 = committedSurfaces[i].size() * conHeig / sqrt(3.0);
                strain2 = 2 * (stress2 - stress1) / elast_plast + strain1;
                gre = stress2 / strain2;
                bb(i, k * 2) = strain2; bb(i, k * 2 + 1) = gre;
            }
        }
    }
    opserr << "SRSMYSand::getBackbone() - out!\n";

}


int SRSMYSand::getResponse(int responseID, Information& matInfo)
{
    opserr << "SRSMYSand::getResponse() - in!\n";
    switch (responseID) {
    case -1:
        return -1;
    case 1:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getCommittedStress();
        return 0;
    case 2:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getCommittedStrain();
        return 0;
    case 3:
        if (matInfo.theMatrix != 0)
            *(matInfo.theMatrix) = getTangent();
        return 0;
    case 4:
        if (matInfo.theMatrix != 0)
            getBackbone(*(matInfo.theMatrix));
        return 0;
        // begin change by Alborz Ghofrani UW --- get 6 components of stress
    case 5:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getStressToRecord(3);
        return 0;
    case 6:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getStressToRecord(4);
        return 0;
    case 7:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getStressToRecord(5);
        return 0;
    case 8:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getStressToRecord(6);
        return 0;
    case 9:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = getStressToRecord(7);
        return 0;
        // end change by Alborz Ghofrani UW
    default:
        return -1;
    }
}

void SRSMYSand::Print(OPS_Stream& s, int flag) { s << "SRSMYSand" << endln; }
void SRSMYSand::setSubTag(const int tag) { subTag = tag; }
int SRSMYSand::getSubTag(void) const { return subTag; }

const Vector& SRSMYSand::getCommittedStress(void)
{
    opserr << "SRSMYSand::getCommittedStress()\n";
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;
    int numOfSurfaces = numOfSurfacesx[matN];
    double residualPress = residualPressx[matN];

    double deviatorRatio = sdevRatio(sv.sig, residualPress);
    double scale = deviatorRatio / committedSurfaces[numOfSurfaces].size();
    if (loadStagex[matN] != 1) scale = 0.;

    static Vector temp7(sv.sig.length());
    for (int i = 0; i < sv.sig.length(); i++)
    {
        temp7[i] = sv.sig(i);
    }
    temp7[sv.sig.length()] = scale;
    return temp7;
}

// begin change by Alborz Ghofrani - UW --- get 6 components of stress
const Vector&
SRSMYSand::getStressToRecord(int numOutput)
{
    opserr << "SRSMYSand::getStressToRecord()\n";
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;

    if (ndm == 3) {
        static Vector temp7(7);
        temp7 = this->getCommittedStress();
        if (numOutput == 6)
        {
            static Vector temp6(6);
            temp6[0] = temp7[0];
            temp6[1] = temp7[1];
            temp6[2] = temp7[2];
            temp6[3] = temp7[3];
            temp6[4] = temp7[4];
            temp6[5] = temp7[5];
            return temp6;
        }
        else if (numOutput == 7)
        {
            return temp7;
        }
        else {
            opserr << "Wrong number of stress components to record!" << endln;
            return temp7;
        }
    }

    else {
        static Vector temp5(5);
        temp5 = this->getCommittedStress();
        if (numOutput == 3)
        {
            static Vector temp3(3);
            temp3[0] = temp5[0];
            temp3[1] = temp5[1];
            temp3[2] = temp5[3];
            return temp3;
        }
        else if (numOutput == 4)
        {
            static Vector temp4(4);
            temp4[0] = temp5[0];
            temp4[1] = temp5[1];
            temp4[2] = temp5[2];
            temp4[3] = temp5[3];
            return temp4;
        }
        else if (numOutput == 5)
        {
            return temp5;
        }
        else {
            opserr << "Wrong number of stress components to record!" << endln;
            return temp5;
        }
    }
}
// end change by Alborz Ghofrani - UW 

const Vector& SRSMYSand::getCommittedStrain(void)
{
    auto& gs = tools::getGlobalStorage(nDOF);
    Vector& strain = gs.e;
    if (beVerbose)  opserr << "SRSMYSand::getCommittedStrain() - " << sv.eps_commit << "\n";
    strain = sv.eps_commit.makeVector();
    return strain;
}


// NOTE: surfaces[0] is not used
void SRSMYSand::setUpSurfaces(double* gredu)
{
    opserr << "SRSMYSand::setUpSurfaces() - in!\n";
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    double frictionAngle = frictionAnglex[matN];
    double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];
    double phaseTransfAngle = phaseTransfAnglex[matN];
    double stressRatioPT = stressRatioPTx[matN];

    double refStrain, peakShear, coneHeight;
    double stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
    double ratio1, ratio2;

    // intialize
    CTensor workV6(6, 2);
    CTensor workT2V(6, 2);

    if (gredu == 0) {
        double sinPhi = sin(frictionAngle * pi / 180.);
        double Mnys = 6. * sinPhi / (3. - sinPhi);
        double sinPhiPT = sin(phaseTransfAngle * pi / 180.);
        stressRatioPT = 6. * sinPhiPT / (3. - sinPhiPT);
        // tao = cohesion * sqrt(8)/3.
        residualPress = 2 * cohesion / Mnys;
        // a small nonzero residualPress for numerical purpose only
        if (residualPress < 0.0001 * pAtm) residualPress = 0.0001 * pAtm;
        coneHeight = -(refPressure - residualPress);
        peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
        refStrain = (peakShearStrain * peakShear)
            / (refShearModulus * peakShearStrain - peakShear);

        double stressInc = peakShear / numOfSurfaces;

        for (int ii = 1; ii <= numOfSurfaces; ii++) {
            stress1 = ii * stressInc;
            stress2 = stress1 + stressInc;
            ratio1 = 3. * stress1 / sqrt(2.) / coneHeight;
            ratio2 = 3. * stress2 / sqrt(2.) / coneHeight;
            strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
            strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);

            if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
                double ratio = (ratio2 - stressRatioPT) / (ratio2 - ratio1);
                strainPTOcta = strain2 - ratio * (strain2 - strain1);
            }

            size = ratio1;
            elasto_plast_modul = 2. * (stress2 - stress1) / (strain2 - strain1);
            if ((2. * refShearModulus - elasto_plast_modul) <= 0)
                plast_modul = UP_LIMIT;
            else
                plast_modul = (2. * refShearModulus * elasto_plast_modul) /
                (2. * refShearModulus - elasto_plast_modul);
            if (plast_modul < 0) plast_modul = 0;
            if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
            if (ii == numOfSurfaces) plast_modul = 0;
            workV6.Zero();
            //opserr<<ii<<" "<<size<<" "<<plast_modul<<endln;
            committedSurfaces[ii] = NestedSurface(workV6, size, plast_modul);
        }  // ii
    }
    else {  //user defined surfaces
        int ii = 2 * (numOfSurfaces - 1);
        double tmax = refShearModulus * gredu[ii] * gredu[ii + 1];
        double Mnys = -(sqrt(3.) * tmax - 2. * cohesion) / refPressure;
        residualPress = 2 * cohesion / Mnys;
        if (residualPress < 0.0001 * pAtm) residualPress = 0.0001 * pAtm;
        coneHeight = -(refPressure - residualPress);

        double sinPhi = 3 * Mnys / (6 + Mnys);
        if (sinPhi < 0. || sinPhi>1.) {
            opserr << "\nNDMaterial " << this->getTag() << ": Invalid friction angle, please modify ref. pressure or G/Gmax curve." << endln;
            exit(-1);
        }

        frictionAngle = asin(sinPhi) * 180 / pi;
        opserr << "\nNDMaterial " << this->getTag() << ": Friction angle is " << frictionAngle << "\n" << endln;
        if (phaseTransfAngle > frictionAngle) {
            opserr << "\nNDMaterial " << this->getTag() << ": phase Transformation Angle > friction Angle,"
                << "will set phase Transformation Angle = friction Angle.\n" << endln;
            phaseTransfAngle = frictionAngle;
        }
        double sinPhiPT = sin(phaseTransfAngle * pi / 180.);
        stressRatioPT = 6. * sinPhiPT / (3. - sinPhiPT);

        for (int i = 1; i < numOfSurfaces; i++) {
            int ii = 2 * (i - 1);
            strain1 = gredu[ii];
            stress1 = refShearModulus * gredu[ii + 1] * strain1;
            strain2 = gredu[ii + 2];
            stress2 = refShearModulus * gredu[ii + 3] * strain2;

            ratio1 = sqrt(3.) * stress1 / coneHeight;
            ratio2 = sqrt(3.) * stress2 / coneHeight;
            if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
                double ratio = (ratio2 - stressRatioPT) / (ratio2 - ratio1);
                // gamma_oct = sqrt(6)/3*gamma12
                strainPTOcta = sqrt(6.) / 3 * (strain2 - ratio * (strain2 - strain1));
            }

            size = ratio1;
            elasto_plast_modul = 2. * (stress2 - stress1) / (strain2 - strain1);

            if ((2. * refShearModulus - elasto_plast_modul) <= 0)
                plast_modul = UP_LIMIT;
            else
                plast_modul = (2. * refShearModulus * elasto_plast_modul) /
                (2. * refShearModulus - elasto_plast_modul);
            if (plast_modul <= 0) {
                opserr << "\nNDMaterial " << this->getTag() << ": Surface " << i
                    << " has plastic modulus < 0.\n Please modify G/Gmax curve.\n" << endln;
                exit(-1);
            }
            if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;

            workV6.Zero();
            //opserr<<size<<" "<<i<<" "<<plast_modul<<" "<<gredu[ii]<<" "<<gredu[ii+1]<<endln;
            committedSurfaces[i] = NestedSurface(workV6, size, plast_modul);

            if (i == (numOfSurfaces - 1)) {
                plast_modul = 0;
                size = ratio2;
                //opserr<<size<<" "<<i+1<<" "<<plast_modul<<" "<<gredu[ii+2]<<" "<<gredu[ii+3]<<endln;
                committedSurfaces[i + 1] = NestedSurface(workV6, size, plast_modul);
            }
        }
    }

    residualPressx[matN] = residualPress;
    frictionAnglex[matN] = frictionAngle;
    cohesionx[matN] = cohesion;
    phaseTransfAnglex[matN] = phaseTransfAngle;
    stressRatioPTx[matN] = stressRatioPT;
    opserr << "SRSMYSand::setUpSurfaces() - out!\n";
}


double SRSMYSand::yieldFunc(CTensor& stress, const std::vector<NestedSurface> surfaces, int surfaceNum)
{
    opserr << "SRSMYSand::yieldFunc() - in!\n";
    double residualPress = residualPressx[matN];
    double currentPress = stress.trace() / nD;
    double coneHeight = currentPress - residualPress;
    CTensor zeta = stress.deviator() -  surfaces[surfaceNum].center() * coneHeight;;
    double sz = surfaces[surfaceNum].size() * coneHeight;
    double yieldf = 3. / 2. * (zeta % zeta) - sz * sz;
    opserr << "SRSMYSand::yieldFunc() - out!\n";
    return yieldf;
}


void SRSMYSand::deviatorScaling(CTensor& stress, std::vector<NestedSurface> surfaces, int surfaceNum)
{
    opserr << "SRSMYSand::deviatorScaling() - in!\n";
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    
    double diff = yieldFunc(stress, surfaces, surfaceNum);
    double coneHeight = stress.trace() / nD - residualPress;
    
    if (surfaceNum < numOfSurfaces && diff < 0.) {
        double sz = -surfaces[surfaceNum].size() * coneHeight;
        double deviaSz = sqrt(sz * sz + diff);
        CTensor devia = stress.deviator();
        CTensor zeta = devia - surfaces[surfaceNum].center() * coneHeight;
        double coeff = (sz - deviaSz) / deviaSz;
        if (coeff < 1.e-13) coeff = 1.e-13;
        devia += zeta * coeff;
        stress.setData(devia, stress.trace() / nD);
        deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
    }

    if (surfaceNum == numOfSurfaces && fabs(diff) > LOW_LIMIT) {
        double sz = -surfaces[surfaceNum].size() * coneHeight;
        CTensor devia = stress.deviator() * sz / sqrt(diff + sz * sz);
        stress.setData(devia, stress.trace() / nD);
    }
    opserr << "SRSMYSand::deviatorScaling() - out!\n";
}


void SRSMYSand::initSurfaceUpdate(void)
{
    opserr << "SRSMYSand::initSurfaceUpdate() - in!\n";
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (sv.nYs_commit == 0) return;

    CTensor devia = sv.sig.deviator();
    CTensor center(nDOF, 2); // 2: contravariant

    double coneHeight = -(sv.sig.trace() / nD - residualPress);

    double Ms = sqrt(3. / 2. * (devia % devia));

    if (sv.nYs_commit < numOfSurfaces) { // failure surface can't move
        center = devia * (1. - committedSurfaces[sv.nYs_commit].size() * coneHeight / Ms);

        //workV6 = workV6 / -coneHeight;
        center /= -coneHeight;
        committedSurfaces[sv.nYs_commit].setCenter(center);
    }

    for (int i = 1; i < sv.nYs_commit; i++) {
        center = devia * (1. - committedSurfaces[i].size() * coneHeight / Ms);
        center /= -coneHeight;
        committedSurfaces[i].setCenter(center);
        theSurfaces[i] = committedSurfaces[i];
    }
    sv.nYs = sv.nYs_commit;
    opserr << "SRSMYSand::initSurfaceUpdate() - out!\n";
}


void SRSMYSand::initStrainUpdate(void)
{
    opserr << "SRSMYSand::initStrainUpdate() - in!\n";
    
    CTensor deviator;

    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];
    double stressRatioPT = stressRatioPTx[matN];

    // elastic strain state
    double stressRatio = sdevRatio(sv.sig, residualPress);
    double ratio = (-sv.sig.trace() / nD + residualPress) / (-refPressure + residualPress);
    ratio = pow(ratio, 1. - pressDependCoeff);
    modulusFactor = getModulusFactor(sv.sig);
    double shearCoeff = 1. / (2. * refShearModulus * modulusFactor);
    double bulkCoeff = 1. / (3. * refBulkModulus * modulusFactor);

    //sv.eps = sv.sig.deviator()*shearCoeff
    //              + sv.sig.trace() / nD*bulkCoeff;

    // modified fmk as discussed with z.yang
    deviator = sv.sig.deviator() * shearCoeff;
    sv.eps.setData(deviator, sv.sig.trace() / nD * bulkCoeff);

    double octalStrain = sv.eps.octahedral();
    if (octalStrain <= LOW_LIMIT) octalStrain = LOW_LIMIT;

    // plastic strain state (scaled from elastic strain)
    double scale, PPZLimit;
    if (stressRatio >= stressRatioPT) {  //above PT
        onPPZCommitted = 2;
        prePPZStrainOctaCommitted = strainPTOcta * ratio;
        PPZLimit = getPPZLimits(1, sv.sig);
        scale = sqrt(prePPZStrainOctaCommitted + PPZLimit) / octalStrain;
    }
    else {  // below PT
        onPPZCommitted = -1;
        prePPZStrainOctaCommitted = octalStrain;
        if (prePPZStrainOctaCommitted > strainPTOcta * ratio)
            prePPZStrainOctaCommitted = strainPTOcta * ratio;
        scale = sqrt(prePPZStrainOctaCommitted) / octalStrain;
    }
    //sv.eps.setData(sv.eps.deviator()*scale, sv.eps.trace() / nD);
    deviator = sv.eps.deviator() * scale;
    sv.eps.setData(deviator, sv.eps.trace() / nD);
    PPZPivotCommitted = sv.eps;
    opserr << "SRSMYSand::initStrainUpdate() - out!\n";
}


double SRSMYSand::getModulusFactor(CTensor& stress)
{
    opserr << "SRSMYSand::getModulusFactor() - in!\n";
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];

    double conHeig = stress.trace() / nD - residualPress;
    double scale = conHeig / (refPressure - residualPress);
    scale = pow(scale, pressDependCoeff);

    opserr << "SRSMYSand::getModulusFactor() - out!\n";
    return (1.e-10 > scale) ? 1.e-10 : scale;
}


void SRSMYSand::setTrialStress(CTensor& stress)
{
    opserr << "SRSMYSand::setTrialStress() - in!\n";
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];

    CTensor deviator(nDOF, 2);  // 2: contravariant

    modulusFactor = getModulusFactor(stress);
    deviator = stress.deviator();
    deviator += subStrainRate.deviator() * 2 * refShearModulus * modulusFactor;

    double B = refBulkModulus * modulusFactor;

    if (Hvx[matN] != 0. && trialStress.trace() / nD <= maxPress
        && subStrainRate.trace() / nD < 0. && loadStagex[matN] == 1) {
        double tp = fabs(trialStress.trace() / nD - residualPressx[matN]);
        B = (B * Hvx[matN] * pow(tp, Pvx[matN])) / (B + Hvx[matN] * pow(tp, Pvx[matN]));
    }

    double volume = stress.trace() / nD + subStrainRate.trace() / nD * 3. * B;

    if (volume > 0.) volume = 0.;
    trialStress.setData(deviator, volume);
    opserr << "SRSMYSand::setTrialStress() - out!\n";
}


int SRSMYSand::setSubStrainRate(void)
{
    opserr << "SRSMYSand::setSubStrainRate() - in!\n";

    CTensor deviator;
    double residualPress = residualPressx[matN];
    double refShearModulus = refShearModulusx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (strainRate.norm() == 0) return 0;

    double elast_plast_modulus;
    double conHeig = -(sv.sig.trace() / nD - residualPress);
    double factor = getModulusFactor(sv.sig);
    if (sv.nYs == 0)
        elast_plast_modulus = 2 * refShearModulus * factor;
    else {
        double plast_modulus = theSurfaces[sv.nYs].modulus() * factor;
        elast_plast_modulus = 2 * refShearModulus * factor * plast_modulus
            / (2 * refShearModulus * factor + plast_modulus);
    }
    deviator = strainRate.deviator() * elast_plast_modulus;

    double singleCross = theSurfaces[numOfSurfaces].size() * conHeig / numOfSurfaces;
    double totalCross = 3. * deviator.octahedral() / sqrt(2.);
    int numOfSub = totalCross / singleCross + 1;
    if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;

    int numOfSub1 = strainRate.octahedral() / 1.0e-5;
    int numOfSub2 = strainRate.trace() / nD / 1.e-5;
    if (numOfSub1 > numOfSub) numOfSub = numOfSub1;
    if (numOfSub2 > numOfSub) numOfSub = numOfSub2;

    deviator = strainRate * (1.0 / numOfSub);

    subStrainRate = deviator;

    opserr << "SRSMYSand::setSubStrainRate() - out!\n";
    return numOfSub;
}


void SRSMYSand::getContactStress(CTensor& contactStress)
{
    opserr << "SRSMYSand::getContactStress() - in!\n";
    double residualPress = residualPressx[matN];
    double conHeig = trialStress.trace() / nD - residualPress;
    CTensor center(6, 2); // 2: contravariant
    CTensor deviator(6, 2); // 2: contravariant
    center = theSurfaces[sv.nYs].center();
    deviator = trialStress.deviator() - center * conHeig;
    double Ms = sqrt(3. / 2. * (deviator % deviator));
    deviator = deviator * theSurfaces[sv.nYs].size() * (-conHeig) / Ms + center * conHeig;
    contactStress.setData(deviator, trialStress.trace() / nD);
    opserr << "SRSMYSand::getContactStress() - out!\n";
}


int SRSMYSand::isLoadReversal(CTensor& stress)
{
    opserr << "SRSMYSand::isLoadReversal() - in!\n";
    if (sv.nYs == 0) return 0;
    CTensor normal(6, 2); // 2: contravariant
    
    getSurfaceNormal(stress, normal);

    //if (((trialStress - sv.sig)
    //	&& workT2V) < 0) return 1;
    CTensor temp = trialStress;
    temp -= sv.sig;

    if ((temp % normal) < 0) return 1;

    opserr << "SRSMYSand::isLoadReversal() - out!\n";
    return 0;
}


void
SRSMYSand::getSurfaceNormal(CTensor& stress, CTensor& normal)
{
    opserr << "SRSMYSand::getSurfaceNormal() - in!\n";
    double residualPress = residualPressx[matN];
    double conHeig = stress.trace() / nD - residualPress;
    CTensor deviator = stress.deviator();
    CTensor center = theSurfaces[sv.nYs].center();
    center = theSurfaces[sv.nYs].center();
    double sz = theSurfaces[sv.nYs].size();
    double volume = conHeig * ((center % center) - 2. / 3. * sz * sz) - (deviator % center);
    deviator += center * -conHeig;
    deviator *= 3.0;
    normal.setData(deviator, volume);
    normal.Normalize();
    opserr << "SRSMYSand::getSurfaceNormal() - out!\n";
}


double SRSMYSand::getPlasticPotential(CTensor& contactStress,
    const CTensor& surfaceNormal)
{
    opserr << "SRSMYSand::getPlasticPotential() - in!\n";
    double residualPress = residualPressx[matN];
    double stressRatioPT = stressRatioPTx[matN];
    double contractParam1 = contractParam1x[matN];
    double contractParam2 = contractParam2x[matN];
    double contractParam3 = contractParam3x[matN];
    double dilateParam1 = dilateParam1x[matN];
    double dilateParam2 = dilateParam2x[matN];
    CTensor deviator_ratio(nDOF, 2);  // 2: contravariant

    double plasticPotential, contractRule, shearLoading, angle;

    double contactRatio = sdevRatio(contactStress, residualPress);
    double factorPT = contactRatio / stressRatioPT;

    double currentRatio = sdevRatio(updatedTrialStress, residualPress);
    double trialRatio = sdevRatio(trialStress, residualPress);

    shearLoading = updatedTrialStress.deviator() % trialStress.deviator();
    //shearLoading = sv.sig.deviator() % trialStress.deviator();

    if (factorPT >= 1. && trialRatio >= currentRatio && shearLoading >= 0.) {  //dilation
        updatePPZ(contactStress);
        if (onPPZ == 1)
            plasticPotential = 0.;
        else if (onPPZ == 2) {
            factorPT -= 1.0;
            double dilateParam3 = dilateParam3x[matN];
            double ppp = pow((fabs(contactStress.trace() / nD) + fabs(residualPress)) / pAtm, -dilateParam3);
            plasticPotential = ppp * factorPT * (factorPT) * (dilateParam1 + pow(cumuDilateStrainOcta, dilateParam2));
            if (plasticPotential < 0.) plasticPotential = -plasticPotential;
            if (plasticPotential > 5.0e4) plasticPotential = 5.0e4;
        }
        else {
            opserr << "FATAL: Wrong onPPZ value: " << onPPZ << endln;
            exit(-1);
        }
    }
    else {  //contraction
        if (currentRatio == 0.) angle = 1.0;
        else {
            deviator_ratio = trialStress.deviator();
            deviator_ratio /= (fabs(trialStress.trace() / nD) + fabs(residualPress));
            deviator_ratio -= updatedTrialStress.deviator() / (fabs(updatedTrialStress.trace() / nD) + fabs(residualPress));
            CTensor temp = deviator_ratio;
            temp = temp.deviator();
            double tempnorm = temp.norm();
            CTensor updatedTrialStressDev = updatedTrialStress.deviator();
            if (tempnorm == 0.) angle = 1.0;
            else angle = (updatedTrialStressDev % deviator_ratio) / tempnorm / updatedTrialStressDev.norm();
        }
        factorPT = factorPT * angle - 1.0;

        contractRule = pow((fabs(contactStress.trace() / nD) + fabs(residualPress)) / pAtm, contractParam3);
        if (contractRule < 0.1) contractRule = 0.1;

        plasticPotential = -(factorPT) * (factorPT) * (contractParam1 + maxCumuDilateStrainOcta * contractParam2) * contractRule;

        if (plasticPotential > 0.) plasticPotential = -plasticPotential;
        if (onPPZ > 0) onPPZ = 0;
        if (onPPZ != -1) PPZTranslation(contactStress);
    }

    if (isCriticalState(contactStress)) plasticPotential = 0;
    opserr << "SRSMYSand::getPlasticPotential() - out!\n";
    return plasticPotential;
}


int SRSMYSand::isCriticalState(CTensor& stress)
{
    opserr << "SRSMYSand::isCriticalState() - in!\n";
    double einit = einitx[matN];
    double volLimit1 = volLimit1x[matN];
    double volLimit2 = volLimit2x[matN];
    double volLimit3 = volLimit3x[matN];

    double vol = trialStrain.trace() / nD * 3.0;
    double etria = einit + vol + vol * einit;
    vol = sv.eps.trace() / nD * 3.0;
    double ecurr = einit + vol + vol * einit;

    double ecr1, ecr2;
    if (volLimit3 != 0.) {
        ecr1 = volLimit1 - volLimit2 * pow(fabs(-stress.trace() / nD / pAtm), volLimit3);
        ecr2 = volLimit1 - volLimit2 * pow(fabs(-updatedTrialStress.trace() / nD / pAtm), volLimit3);
    }
    else {
        ecr1 = volLimit1 - volLimit2 * log(fabs(-stress.trace() / nD / pAtm));
        ecr2 = volLimit1 - volLimit2 * log(fabs(-updatedTrialStress.trace() / nD / pAtm));
    }

    if (ecurr < ecr2 && etria < ecr1) return 0;
    if (ecurr > ecr2 && etria > ecr1) return 0;
    opserr << "SRSMYSand::isCriticalState() - out!\n";
    return 1;
}


void SRSMYSand::updatePPZ(CTensor& contactStress)
{
    opserr << "SRSMYSand::updatePPZ() - in!\n";
    double liquefyParam1 = liquefyParam1x[matN];
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff = pressDependCoeffx[matN];
    double liquefyParam2 = liquefyParam2x[matN];

    // onPPZ=-1 may not be needed. can start with onPPZ=0  ****

    double temp = strainRate.deviator() % PivotStrainRateCommitted;
    check = strainRate.deviator()(3);

    if (onPPZ < 1) {
        damage = 0.0;
        if ((maxPress - sv.sig.trace() / nD) / (maxPress - residualPress) > 0.)
            damage = pow((maxPress - sv.sig.trace() / nD) / (maxPress - residualPress), 0.25);
    }

    // PPZ inactive if liquefyParam1==0.
    if (liquefyParam1 == 0. || (onPPZ < 1 && damage < 0.)) {
        if (onPPZ == 2) {
            PPZPivot = trialStrain;
            //PivotStrainRate = strainRate.deviator();
            cumuDilateStrainOcta += subStrainRate.octahedral();
        }
        else if (onPPZ != 2) {
            onPPZ = 2;
            PPZPivot = trialStrain;
            PivotStrainRate = strainRate.deviator();
            if (temp < 0.) cumuDilateStrainOcta = 0.;
            //if (temp < 0.) maxCumuDilateStrainOcta = 0.;
        }
        return;
    }

    // dilation: calc. cumulated dilative strain
    if (onPPZ == 2) {
        PPZPivot = trialStrain;
        //PivotStrainRate = strainRate.deviator();
        cumuDilateStrainOcta += subStrainRate.octahedral();

        double zzz = 0.;
        if (damage > zzz) zzz = damage;
        maxCumuDilateStrainOcta += zzz * liquefyParam1 * subStrainRate.octahedral();
        return;
    }

    if (onPPZ == -1 || onPPZ == 0) {

        // moved to the opposite direction, update prePPZStrainOcta and
        // oppoPrePPZStrainOcta (they are small values, may consider drop these
        // parameters altogether).

        if (temp < 0.) {
            double volume = -contactStress.trace() / nD;
            oppoPrePPZStrainOcta = prePPZStrainOcta;
            double ratio = (volume + residualPress) / (-refPressure + residualPress);
            ratio = pow(ratio, 1. - pressDependCoeff);
            prePPZStrainOcta = ratio * strainPTOcta;
            if (oppoPrePPZStrainOcta == 0.) oppoPrePPZStrainOcta = prePPZStrainOcta;
        }
    }

    //double dd = maxCumuDilateStrainOcta > cumuDilateStrainOcta ? maxCumuDilateStrainOcta : cumuDilateStrainOcta;

    //PPZSize = (cumuTranslateStrainOcta+dd)/2.;
    PPZSize = (cumuTranslateStrainOcta + maxCumuDilateStrainOcta) / 2.;
    //(prePPZStrainOcta+oppoPrePPZStrainOcta+cumuTranslateStrainOcta+maxCumuDilateStrainOcta)/2.;

 // calc. new PPZ center.
    CTensor PPZ;
    if (onPPZ == 0 || (onPPZ == 1 && temp < 0.0)) {
        // new center lies on the vector of PPZPivot-PPZCenter,
        // its distance from PPZPivot is (PPZSize-cumuTranslateStrainOcta).

        //workT2V.setData(PPZPivot - PPZCenter);
        PPZ = PPZPivot - PPZCenter;

        double coeff;
        if (PPZ.octahedral() == 0.) coeff = 0.;
        else coeff = (PPZSize - cumuTranslateStrainOcta) / PPZ.octahedral();
        PPZCenter = PPZPivot - PPZ * coeff;
    }

    //workT2V.setData(trialStrain - PPZCenter);
    PPZ = trialStrain - PPZCenter;

    //outside PPZ
    //if (workT2V.octahedral(1) > PPZSize && temp > 0. || PPZLimit==0.) {
    if (PPZ.octahedral() > PPZSize) {

        // rezero cumuDilateStrainOcta only when load reverses
        // (partial unloading does not erase cumuDilateStrainOcta.

        //if (temp < 0.) {
        cumuDilateStrainOcta = 0.;
        //if (PPZLimit == 0.) maxCumuDilateStrainOcta = 0.;
   //}
        onPPZ = 2;
        PPZPivot = trialStrain;
        PivotStrainRate = strainRate.deviator();
        cumuTranslateStrainOcta = 0.;
    }
    else {  //inside PPZ
        if (onPPZ == 0 || onPPZ == 1) PPZTranslation(contactStress);
        //if (onPPZ == 0) onPPZ = 1;
        if (onPPZ == -1 || onPPZ == 0) onPPZ = 1;
    }
    opserr << "SRSMYSand::updatePPZ() - out!\n";
}


void SRSMYSand::PPZTranslation(const CTensor& contactStress)
{
    opserr << "SRSMYSand::PPZTranslation() - in!\n";
    double liquefyParam1 = liquefyParam1x[matN];
    double liquefyParam2 = liquefyParam2x[matN];
    double residualPress = residualPressx[matN];

    //cumuDilateStrainOcta -= subStrainRate.octahedral(1);
    //if (cumuDilateStrainOcta < 0.) cumuDilateStrainOcta = 0.;

    if (liquefyParam1 == 0.) return;

    //Amount of translation is proportional to the amount of previous unloading,
    //and is also limited by the amount of previous dilation (no translation
    //if there was no dilation), as damage is really caused by dilation.
    damage = 0.0;
    if ((maxPress - sv.sig.trace() / nD) / (maxPress - residualPress) > 0.)
        damage = pow((maxPress - sv.sig.trace() / nD) / (maxPress - residualPress), 0.25);

    double zzz = 0.;
    if (damage > zzz) zzz = damage;

    double temp = strainRate.deviator() % PivotStrainRateCommitted;

    if (temp < 0.0) {  //update only when load reverses

        //workT2V.setData(trialStrain.deviator()-PPZPivot.deviator());
        CTensor deviator = trialStrain.deviator();
        deviator -= PPZPivot.deviator();

        temp = deviator.octahedral();
        if (cumuTranslateStrainOcta < zzz * liquefyParam2 * temp)
            cumuTranslateStrainOcta = zzz * liquefyParam2 * temp;
        //if (maxCumuDilateStrainOcta == 0.) temp = 0.; //PPZTransLimit;
          //\\// else temp = PPZTransLimit*cumuDilateStrainOcta/maxCumuDilateStrainOcta;
        //else temp = dilateParam3*cumuDilateStrainOcta/maxCumuDilateStrainOcta;
        //if (cumuTranslateStrainOcta > temp) cumuTranslateStrainOcta = temp;

        //\\// cumuTranslateStrainOcta = dilateParam3*cumuDilateStrainOcta;
    }
    opserr << "SRSMYSand::PPZTranslation() - out!\n";
}


double SRSMYSand::getPPZLimits(int which, CTensor& contactStress)
{
    opserr << "SRSMYSand::getPPZLimits() - in!\n";
    double liquefyParam1 = liquefyParam1x[matN];
    double liquefyParam2 = liquefyParam2x[matN];
    double dilateParam3 = dilateParam3x[matN];

    double PPZLimit, temp;
    double volume = -contactStress.trace() / nD;

    if (volume >= liquefyParam1) PPZLimit = 0.;
    else {
        temp = volume * pi / liquefyParam1 / 2.;
        // liquefyParam3 = 3.0 default
        PPZLimit = liquefyParam2 * pow(cos(temp), 3.);
        //\\//
        PPZLimit = 0.0;
    }

    opserr << "SRSMYSand::getPPZLimits() - out!\n";
    if (which == 1)
        return PPZLimit;
    else if (which == 2)
        return dilateParam3 * PPZLimit;
    else {
        opserr << "FATAL:SRSMYSand::getPPZLimits: unknown argument value" << endln;
        exit(-1);
        return 0.0;
    }
}


double SRSMYSand::getLoadingFunc(CTensor& contactStress, CTensor& surfaceNormal, double* plasticPotential, int crossedSurface)
{
    opserr << "SRSMYSand::getLoadingFunc() - in!\n";
    int numOfSurfaces = numOfSurfacesx[matN];
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];

    double loadingFunc, limit;
    double modul = theSurfaces[sv.nYs].modulus();
    double temp1 = 2. * refShearModulus * modulusFactor * (surfaceNormal.deviator() % surfaceNormal.deviator());
    double temp2 = 9. * refBulkModulus * modulusFactor * surfaceNormal.trace() / nD * (*plasticPotential);

    //for the first crossing
    double temp = temp1 + temp2 + modul * modulusFactor;
    if (sv.nYs == numOfSurfaces)
        limit = theSurfaces[sv.nYs - 1].modulus() * modulusFactor / 2.;
    else limit = modul * modulusFactor / 2.;
    if (temp < limit) {
        (*plasticPotential) = (temp2 + limit - temp) / (9. * refBulkModulus * modulusFactor
            * surfaceNormal.trace() / nD);
        temp = limit;
    }
    //loadingFunc = (surfaceNormal
    //	             && (trialStress.deviator()-contactStress.deviator()))/temp;
    CTensor deviator = trialStress.deviator();
    deviator -= contactStress.deviator();
    loadingFunc = (surfaceNormal % deviator) / temp;

    if (loadingFunc < 0.) loadingFunc = 0;

    //for more than one crossing
    if (crossedSurface) {
        temp = (theSurfaces[sv.nYs - 1].modulus() - modul)
            / theSurfaces[sv.nYs - 1].modulus();
        loadingFunc *= temp;
    }

    opserr << "SRSMYSand::getLoadingFunc() - out!\n";
    return loadingFunc;
}


int SRSMYSand::stressCorrection(int crossedSurface)
{
    opserr << "SRSMYSand::stressCorrection() - in!\n";
    double refShearModulus = refShearModulusx[matN];
    double refBulkModulus = refBulkModulusx[matN];

    static CTensor contactStress;
    getContactStress(contactStress);
    static CTensor surfNormal;
    getSurfaceNormal(contactStress, surfNormal);
    double plasticPotential = getPlasticPotential(contactStress, surfNormal);
    double tVolume = trialStress.trace() / nD;
    double loadingFunc = getLoadingFunc(contactStress, surfNormal,
        &plasticPotential, crossedSurface);
    double volume = tVolume - plasticPotential * 3 * refBulkModulus * modulusFactor * loadingFunc;

    //workV6 = trialStress.deviator()
    //	         - surfNormal.deviator()*2*refShearModulus*modulusFactor*loadingFunc;
    CTensor deviator = trialStress.deviator();

    if (volume > 0. && volume != tVolume) {
        double coeff = tVolume / (tVolume - volume);
        coeff *= -2 * refShearModulus * modulusFactor * loadingFunc;
        deviator += surfNormal.deviator() * coeff;
        volume = 0.;
    }
    else if (volume > 0.) {
        volume = 0.;
    }
    else {
        double coeff = -2 * refShearModulus * modulusFactor * loadingFunc;
        deviator += surfNormal.deviator() * coeff;
    }
    /*
      if (volume>0.)volume = 0.;
        double coeff = -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
  */
    trialStress.setData(deviator, volume);
    deviatorScaling(trialStress, theSurfaces, sv.nYs);

    if (isCrossingNextSurface()) {
        sv.nYs++;
        return stressCorrection(1);  //recursive call
    }

    opserr << "SRSMYSand::stressCorrection() - out!\n";
    return 0;
}


void SRSMYSand::updateActiveSurface(void)
{
    opserr << "SRSMYSand::updateActiveSurface() - in!\n";
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

    if (sv.nYs == numOfSurfaces) return;

    double A, B, C, X;
    static CTensor t1(6, 2);        // 2: contravariant
    static CTensor t2(6, 2);        // 2: contravariant
    static CTensor center(6, 2);    // 2: contravariant
    static CTensor outcenter(6, 2); // 2: contravariant
    double conHeig = trialStress.trace() / nD - residualPress;
    center = theSurfaces[sv.nYs].center();
    double size = theSurfaces[sv.nYs].size();
    outcenter = theSurfaces[sv.nYs + 1].center();
    double outsize = theSurfaces[sv.nYs + 1].size();

    //t1 = trialStress.deviator() - center*conHeig;
    //t2 = (center - outcenter)*conHeig;
    t1 = trialStress.deviator();
    t1 += center * -conHeig;
    t2 = center;
    t2 -= outcenter;
    t2 *= conHeig;

    A = t1 % t1;
    B = 2. * (t1 % t2);
    C = (t2 % t2) - 2. / 3. * outsize * outsize * conHeig * conHeig;
    X = secondOrderEqn(A, B, C, 0);
    if (fabs(X - 1.) < LOW_LIMIT) X = 1.;
    if (X < 1.) return;

    if (X < 1.) {
        //t2 = trialStress.deviator() - outcenter*conHeig;
        t2 = trialStress.deviator();
        t2 += outcenter * -conHeig;

        double xx1 = (t2 % t2) - 2. / 3. * outsize * outsize * conHeig * conHeig;
        double xx2 = (t1 % t1) - 2. / 3. * size * size * conHeig * conHeig;
        opserr << "FATAL:SRSMYSand::updateActiveSurface(): error in Direction of surface motion." << endln;
        opserr << "X-1= " << X - 1 << " A= " << A << " B= " << B << " C= " << C << " M= " << sv.nYs << " low_limit=" << LOW_LIMIT << endln;
        opserr << "diff1= " << xx1 << " diff2= " << xx2 << " p= " << conHeig << " size= " << size << " outs= " << outsize << endln;
        exit(-1);
    }

    //workV6 = (t1 * X + center*conHeig) * (1. - size / outsize)
    //	     - (center - outcenter * size / outsize) * conHeig;

    CTensor workV6 = t1 * X;
    workV6 += center * conHeig;
    workV6 *= (1.0 - size / outsize);
    t2 = center;
    t2 += outcenter * (-size / outsize);
    t2 *= conHeig;
    workV6 -= t2;

    CTensor workT2V = workV6.deviator();
    if (workT2V.norm() < LOW_LIMIT) return;

    workV6 = workT2V;
    A = conHeig * conHeig * (workV6 % workV6);
    B = 2 * conHeig * (t1 % workV6);
    if (fabs(B) < LOW_LIMIT) B = 0.;
    C = (t1 % t1) - 2. / 3. * size * size * conHeig * conHeig;
    if (fabs(C) < LOW_LIMIT || fabs(C) / (t1 % t1) < LOW_LIMIT) return;
    if (B > 0. || C < 0.) {
        opserr << "FATAL:SRSMYSand::updateActiveSurface(): error in surface motion.\n"
            << "A= " << A << " B= " << B << " C= " << C << " (t1&&t1)= " << (t1 % t1) << endln;
        exit(-1);
    }
    X = secondOrderEqn(A, B, C, 1);

    center += workV6 * -X;
    theSurfaces[sv.nYs].setCenter(center);
    opserr << "SRSMYSand::updateActiveSurface() - out!\n";
}


void SRSMYSand::updateInnerSurface(void)
{
    opserr << "SRSMYSand::updateInnerSurface() - in!\n";
    double residualPress = residualPressx[matN];

    if (sv.nYs <= 1) return;
    CTensor devia(6, 2);     // 2: contravariant
    CTensor center(6, 2);    // 2: contravariant
    CTensor temp;

    double conHeig = sv.sig.trace() / nD - residualPress;
    devia = sv.sig.deviator();
    center = theSurfaces[sv.nYs].center();
    double size = theSurfaces[sv.nYs].size();

    for (int i = 1; i < sv.nYs; i++) {
        temp = center * conHeig;
        temp -= devia;
        temp *= theSurfaces[i].size() / size;
        temp += devia;

        //workV6 = devia - (devia - center*conHeig) * theSurfaces[i].size() / size;

        temp /= conHeig;
        theSurfaces[i].setCenter(temp);
    }
    opserr << "SRSMYSand::updateInnerSurface() - out!\n";
}

double SRSMYSand::sdevRatio(CTensor& stress, double residualPress) {
    opserr << "SRSMYSand::sdevRatio()\n";
    CTensor dev = stress.deviator();
    return sqrt(3. / 2. * (dev.norm())) / (abs(stress.trace() / nD) + abs(residualPress));
}

int SRSMYSand::isCrossingNextSurface(void)
{
    opserr << "SRSMYSand::isCrossingNextSurface()\n";
    int numOfSurfaces = numOfSurfacesx[matN];

    if (sv.nYs == numOfSurfaces) return 0;

    if (yieldFunc(trialStress, theSurfaces, sv.nYs + 1) > 0) return 1;

    return 0;
}

double SRSMYSand::secondOrderEqn(double A, double B, double C, int i)
{
    opserr << "SRSMYSand::secondOrderEqn() - in!\n";
    if (A == 0) {
        opserr << "FATAL:second_order_eqn: A=0." << endln;
        if (i == 0) opserr << " when finding reference point on outer surface." << endln;
        else opserr << " when moving active surface." << endln;
        exit(-1);
    }
    if (C == 0) return 0;
    if (B == 0) {
        if (C / A > 0) {
            opserr << "FATAL:second_order_eqn: Complex roots.\n";
            exit(-1);
        }
        return sqrt(-C / A);
    }

    double determ, val1, val2, val;
    determ = B * B - 4. * A * C;
    if (determ < 0) {
        opserr << "FATAL:second_order_eqn: Complex roots.\n";
        if (i == 0) opserr << " when finding reference point on outer surface." << endln;
        else opserr << " when moving active surface." << endln;
        opserr << "B2=" << B * B << " 4AC=" << 4. * A * C << endln;
        exit(-1);
    }

    if (B > 0) val1 = (-B - sqrt(determ)) / (2. * A);
    else val1 = (-B + sqrt(determ)) / (2. * A);
    val2 = C / (A * val1);

    if (val1 < 0 && val2 < 0) {
        if (fabs(val1) < LOW_LIMIT) val1 = 0.;
        else if (fabs(val2) < LOW_LIMIT) val2 = 0.;
    }

    if (val1 < 0 && val2 < 0) {
        opserr << "FATAL:second_order_eqn: Negative roots.\n";
        if (i == 0) opserr << " when finding reference point on outer surface." << endln;
        else opserr << " when moving active surface." << endln;
        opserr << "A=" << A << " B=" << B << " C=" << C << " det=" << determ <<
            " x1=" << val1 << " x2=" << val2 << endln;
        exit(-1);
    }

    if (val1 < 0) return  val2;
    else if (val2 < 0) return  val1;
    else {
        val = val1;
        if (val > val2)  val = val2;
        return val;
    }
    opserr << "SRSMYSand::secondOrderEqn() - out!\n";
}

void SRSMYSand::updateInternal(const bool do_implex, const bool do_tangent) {

    opserr << "SRSMYSand::updateInternal() - in!\n";

    int loadStage = loadStagex[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 3;

    int i, is;
    if (loadStage == 1 && e2p == 0) {
        initPress = sv.sig.trace() / nD;
        elast2Plast();
    }

    if (loadStage != 1) {  //linear elastic
        getTangent();
        auto& gs = tools::getGlobalStorage(nDOF);
        CTensor theTangent (gs.Ct, 2);  // 2: contravariant
        sv.sig += theTangent ^ (sv.eps - sv.eps_commit);

    }
    else {
        for (i = 1; i <= numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
        sv.nYs = sv.nYs_commit;
        pressureD = pressureDCommitted;
        onPPZ = onPPZCommitted;
        PPZSize = PPZSizeCommitted;
        cumuDilateStrainOcta = cumuDilateStrainOctaCommitted;
        maxCumuDilateStrainOcta = maxCumuDilateStrainOctaCommitted;
        cumuTranslateStrainOcta = cumuTranslateStrainOctaCommitted;
        prePPZStrainOcta = prePPZStrainOctaCommitted;
        oppoPrePPZStrainOcta = oppoPrePPZStrainOctaCommitted;
        PPZPivot = PPZPivotCommitted;
        PivotStrainRate = PivotStrainRateCommitted;
        PPZCenter = PPZCenterCommitted;

        subStrainRate = strainRate;
        setTrialStress(sv.sig);
        if (sv.nYs > 0 && isLoadReversal(sv.sig)) {
            updateInnerSurface();
            sv.nYs = 0;
        }

        if (sv.nYs == 0 && !isCrossingNextSurface()) {
            trialStrain = sv.eps;
        }
        else {
            int numSubIncre = setSubStrainRate();

            for (i = 0; i < numSubIncre; i++) {
                trialStrain = sv.eps_commit + subStrainRate * (i + 1);

                if (i == 0) {
                    updatedTrialStress = sv.sig;
                    setTrialStress(sv.sig);
                    is = isLoadReversal(sv.sig);
                }
                else {
                    updatedTrialStress = trialStress;
                    setTrialStress(trialStress);
                    is = isLoadReversal(trialStress);
                }

                if (sv.nYs > 0 && is) {
                    updateInnerSurface();
                    sv.nYs = 0;
                }
                if (sv.nYs == 0 && !isCrossingNextSurface()) continue;
                if (sv.nYs == 0) sv.nYs++;
                stressCorrection(0);
                updateActiveSurface();

                double refBulkModulus = refBulkModulusx[matN];
                //modulusFactor was calculated in setTrialStress
                double B = refBulkModulus * modulusFactor;
                //double deltaD = 3.*subStrainRate.trace() / nD
                //	         - (trialStress.trace() / nD-updatedTrialStress.trace() / nD)/B;
                //if (deltaD<0) deltaD /=2 ;
                //pressureD += deltaD;
                pressureD += 3. * subStrainRate.trace() / nD
                    - (trialStress.trace() / nD - updatedTrialStress.trace() / nD) / B;
                if (pressureD < 0.) pressureD = 0.;
                //opserr<<i<<" "<<sv.nYs<<" "<<is<<" "<<subStrainRate[3]<<endln;
            }
        }
    }

    opserr << "SRSMYSand::updateInternal() - out!\n";


}