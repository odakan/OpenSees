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
// Based on: fmk's Vector and Matrix classes
//
// Description: This file contains the class definition for CTensor (Compressed Tensor).
// CTensor is a concrete class implementing the compressed tensor abstraction.
// CTensor class is used to provide the abstraction for the 2nd and 4th order
// material tensors and tensor operations. Depending on the symmetries,
// CTensor object can store and manipulate a tensor in its corresponding 
// compressed matrix represantation. For more information on compressed tensor 
// operations, please see the following reference:
// 
// Helnwein, P. (2001). Some remarks on the compressed matrix representation of symmetric 
//		second-order and fourth-order tensors. Computer Methods in Applied Mechanics and 
//		Engineering, 190(22), 2753–2770. https://doi.org/10.1016/S0045-7825(00)00263-2
//

#include <Vector.h>
#include <Matrix.h>

class CTensor {
public:
    // constructors
    CTensor(void);
    CTensor(const CTensor& T);

    // second-order tensor
    CTensor(int nrows, char* rep);
    CTensor(double* data, int nrows, char* rep);
    CTensor(const Vector& V, char* rep);

    // fourth-order tensor
    CTensor(int nrows, int ncols, char* rep);
    CTensor(double* data, int nrows, int ncols, char* rep);
    CTensor(const Matrix& M, char* rep);

    // destructor
    ~CTensor();

    // utility methods
        // general
    inline int length(void) const;
    inline int noRows(void) const;
    inline int noCols(void) const;
    void Zero(void);
    double Norm(void) const;
    double pNorm(int p) const;
    int Normalize(void);

    // second-order tensor
    int setData(double* newData, int size);
    int resize(int numRow);
    int resize(const Vector& V);
    Vector diagonal() const;

    // fourth-order tensor
    int setData(double* newData, int nRows, int nCols);
    int resize(int numRow, int numCol);
    int resize(const Matrix& M);

//#ifdef USE_CXX11
//    CTensor(CTensor&& T);
//    CTensor(Vector&& V, char* rep);
//    CTensor(Matrix&& M, char* rep);
//    int resize(Vector&& V);
//    int resize(Matrix&& M);
//#endif

    // Symmetric Tensor Operations
    double GetTrace(const Vector& v);
    Vector GetDevPart(const Vector& aV);
    double DoubleDot2_2_Contr(const Vector& v1, const Vector& v2);
    double DoubleDot2_2_Cov(const Vector& v1, const Vector& v2);
    double DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2);
    double GetNorm_Contr(const Vector& v);
    double GetNorm_Cov(const Vector& v);
    Matrix Dyadic2_2(const Vector& v1, const Vector& v2);
    Vector DoubleDot4_2(const Matrix& m1, const Vector& v1);
    Vector DoubleDot2_4(const Vector& v1, const Matrix& m1);
    Matrix DoubleDot4_4(const Matrix& m1, const Matrix& m2);
    Vector ToContraviant(const Vector& v1);
    Vector ToCovariant(const Vector& v1);

    // overloaded operators
    inline double& operator()(int row);
    inline double operator()(int row) const;
    inline double& operator()(int row, int col);
    inline double operator()(int row, int col) const;

    Matrix operator()(const ID& rows, const ID& cols) const;
    Matrix& operator=(const Matrix& M);

    // matrix operations which will preserve the derived type and
    // which can be implemented efficiently without many constructor calls.

    // matrix-scalar operations
    Matrix& operator+=(double fact);
    Matrix& operator-=(double fact);
    Matrix& operator*=(double fact);
    Matrix& operator/=(double fact);

    // matrix operations which generate a new Matrix. They are not the
    // most efficient to use, as constructors must be called twice. They
    // however are useful for matlab like expressions involving Matrices.

    // matrix-scalar operations
    Matrix operator+(double fact) const;
    Matrix operator-(double fact) const;
    Matrix operator*(double fact) const;
    Matrix operator/(double fact) const;

    // matrix-vector operations
    Vector operator*(const Vector& V) const;
    Vector operator^(const Vector& V) const;

    // matrix-matrix operations
    Matrix operator+(const Matrix& M) const;
    Matrix operator-(const Matrix& M) const;
    Matrix operator*(const Matrix& M) const;
    //     Matrix operator/(const Matrix &M) const;    
    Matrix operator^(const Matrix& M) const;
    Matrix& operator+=(const Matrix& M);
    Matrix& operator-=(const Matrix& M);

    // methods to read/write to/from the matrix
    void Output(OPS_Stream& s) const;
    //    void Input(istream &s);


private:

    static double CTENSOR_NOT_VALID_ENTRY;
    static double* ctensorWork;
    static int* intWork;
    static int sizeDoubleWork;
    static int sizeIntWork;

    int order = 0;
    int numRows = 0;
    int numCols = 0;
    int dataSize = 0;
    double* theData = nullptr;
    char* representation = nullptr;
    int fromFree = 0;

public:
    // tensor valued constants
    class Constants {
        static const CTensor I(const int N);
        static const CTensor I4(const int N);
        static const CTensor IIvol(const int N);
        static const CTensor IIdev(const int N);
    };
};
#endif