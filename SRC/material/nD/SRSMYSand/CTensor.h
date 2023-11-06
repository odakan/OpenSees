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

// $Revision: 1.0 $
// $Date: 2023/11/01 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/SRSMYSand/CTensor.h,v $

#ifndef CTensor_h
#define CTensor_h

// Written: Onur Deniz Akan
// Created: 11/23
// Based on fmk's Matrix class and UWMaterials' symmetric tensor operations suite
//
// Description: This file contains the class definition for CTensor (Compressed Tensor).
// CTensor is a concrete class implementing the compressed tensor abstraction.
// CTensor class is used to provide the abstraction for the 2nd and 4th order
// symmetric tensors and tensor operations.
// CTensor object can store and manipulate 2nd and 4th order tensors in their
// corresponding (inital or resulting) compressed matrix representation.
// Rule of thumb:
//      Stress and Stiffness-like tensors -> contravariant
//      Strain-like tensors -> covariant
// For more information on compressed tensor operations, please see the following
// reference:
// 
// Helnwein, P. (2001). Some remarks on the compressed matrix representation of symmetric 
//		second-order and fourth-order tensors. Computer Methods in Applied Mechanics and 
//		Engineering, 190(22), 2753–2770. https://doi.org/10.1016/S0045-7825(00)00263-2
//

#include <Matrix.h>
#include <Vector.h>

class CTensor 
{
public:
    // constructors
    CTensor(void);
    CTensor(const CTensor& C);

    // second-order tensor
    CTensor(int nRows, int rep);
    CTensor(double* data, int nRows, int rep);
    CTensor(const Vector& V, int rep);

    // fourth-order tensor
    CTensor(int nRows, int nCols, int rep);
    CTensor(double* data, int nRows, int nCols, int rep);
    CTensor(const Matrix& M, int rep);

    // destructor
    ~CTensor();

    // utility methods
        // general
    void Zero(void);
    int setOrder(int ord);
    int makeRep(int rep);  // mind that makeRep may trigger switch-representation
    int getRep(void);
    int getOrder(void);
    inline int length(void) const;
    inline int noRows(void) const;
    inline int noCols(void) const;
    double trace(void);
    CTensor deviator(void);
    Vector makeVector(void);
    Matrix makeMatrix(void);

        // second-order tensor
    int setData(double* newData, int size, int rep);
    int resize(int nRows);
    int resize(const Vector& V);

        // fourth-order tensor
    int setData(double* newData, int nRows, int nCols, int rep);
    int resize(int nRows, int nCols);
    int resize(const Matrix& M);

        // from CTensor
    int setData(const CTensor& deviatoric, const double volumetric);

    // symmetric tensor operations
    int Normalize(void);
    double det(void) const;
    double norm(void) const;
    double operator%(const CTensor& C) const;       // % : [2-to-2] order double dot
    CTensor operator^(const CTensor& C) const;      // ^ : [2-to-4, 4-to-2 or 4-to-4] order double dot
    CTensor operator*(const CTensor& C) const;      // * : [2-to-2] order dyadic product
    CTensor& operator^=(const CTensor& C);
    CTensor& operator*=(const CTensor& C);
    CTensor& dot(const CTensor& C);                 // single dot product
    CTensor& invert(void);

    // special norms
    double J2(void);
    double octahedral(void);    // tangential component

    // fast and efficient operations
    //////////////////////////////////////////////////////////////////////////
    int addVector(double factThis, const Vector& other, double factOther);
    int addMatrixVector(double factThis, const Matrix& m, const Vector& v, double factOther);
    int addMatrixTransposeVector(double factThis, const Matrix& m, const Vector& v, double factOther);
    ////////////////////////////////////////////////////////////////////////////
    int addTensor(double factThis, const CTensor& other, double factOther);
    int addTensorTranspose(double factThis, const CTensor& other, double factOther);
    int addTensorProduct(double factThis, const CTensor& A, const CTensor& B, double factOther); // AB
    int addTensorTransposeProduct(double factThis, const CTensor& A, const CTensor& B, double factOther); // A'B
    int addTensorTripleProduct(double factThis, const CTensor& A, const CTensor& B, double factOther); // A'BA
    int addTensorTripleProduct(double factThis, const CTensor& A, const CTensor& B, const CTensor& C, double otherFact); //A'BC

    // overloaded operators
    inline double& operator()(int row);
    inline double operator()(int row) const;
    inline double& operator()(int row, int col);
    inline double operator()(int row, int col) const;
    CTensor operator()(const ID& rows, const ID& cols) const;
    CTensor& operator=(const CTensor& M);

    // ctensor-scalar operations
    CTensor operator+(double fact) const;
    CTensor operator-(double fact) const;
    CTensor operator*(double fact) const;
    CTensor operator/(double fact) const;
    CTensor& operator+=(double fact);
    CTensor& operator-=(double fact);
    CTensor& operator*=(double fact);
    CTensor& operator/=(double fact);

    // other ctensor-ctensor operations
    CTensor operator+(const CTensor& C) const;
    CTensor operator-(const CTensor& C) const;
    CTensor& operator+=(const CTensor& C);
    CTensor& operator-=(const CTensor& C);

    // output methods
    friend OPS_Stream& operator<<(OPS_Stream& s, const CTensor& C);

private:
    int dim = 0;
    int order = 0;
    int repr = -1;   // -1: Unset, Full: 0, Cov: 1, Contr: 2, CovContr: 3 and ContrCov: 4
    int numRows = 0;
    int numCols = 0;
    Matrix ct = Matrix(1, 1);

private:
    // switch-representation
    int toCov(void);
    int toFull(void);
    int toContr(void);
    int toCovContr(void);
    int toContrCov(void);

    // private utilities
    inline void matrix_dim(int nRows);

public:
    // tensor valued constants
    class Constants {
    public:
        static const CTensor I(const int dim, int rep);         // 2nd order identity tensor
        static const CTensor IIsymm(const int dim, int rep);    // 4th order symmetric tensor
        //static const CTensor IIskew(const int dim, int rep);    // 4th order skew-symmetric tensor
        static const CTensor IIvol(const int dim, int rep);     // 4th order volumetric operator
        static const CTensor IIdev(const int dim, int rep);     // 4th order deviatoric operator
    };
};
#endif