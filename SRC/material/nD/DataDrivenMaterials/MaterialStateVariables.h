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


#ifndef MaterialStateVariables_h
#define MaterialStateVariables_h


#include <map>
#include <stdio.h>
#include <stdlib.h> 
#include <Vector.h>
#include <Matrix.h>


namespace tools {
	class GlobalStorage {
	public:
		int size = 0;
		Matrix K;	//stiffness
		Matrix K0;	//initial stiffness
		Vector p;	//stress
		Vector u;	//strain


	public:
		GlobalStorage() = default;
		GlobalStorage& resize(int N) {
			if (N != size) {
				K.resize(N, N);
				K0.resize(N, N);
				p.resize(N);
				u.resize(N);
			}
			return *this;
		}
	};


	static GlobalStorage& getGlobalStorage(int N) {
		static std::map<int, GlobalStorage> gsmap;
		return gsmap[N].resize(N);
	}

	inline double macaulay(double A) {
		double val = (A + abs(A)) * 0.5;
		return val;
	}
}


class MaterialStateVariables {
public:
	// state variables
	Vector eps = Vector(6);					// material strain 
	Vector eps_commit = Vector(6);			// committed material strain
	Vector eps_commit_old = Vector(6);		// previously committed material strain
	Vector sig = Vector(6);					// effective stress
	Vector sig_commit = Vector(6);			// committed effective stress
	Vector sig_commit_old = Vector(6);		// previously committed effective stress
	Vector xs = Vector(6);					// equivalent plastic strain
	Vector xs_commit = Vector(6);			// committed equivalent plastic strain
	Vector xs_commit_old = Vector(6);		// previously committed equivalent plastic strain
	double dtime_n = 0.0;					// time factor
	double dtime_n_commit = 0.0;			// committed time factor
	bool dtime_is_user_defined = false;
	bool dtime_first_set = false;
	// modulus
	Matrix Ce = Matrix(6, 6);				// elastic modulus matrix
	Matrix Cep = Matrix(6, 6);				// tangent modulus matrix
	Vector sig_implex = Vector(3);			// only for output
	double Kmod = 0;						// Current bulk modulus
	double Gmod = 0;						// Current shear modulus
	double Hmod = 0;						// Current plastic modulus

public:
	// null constructor
	MaterialStateVariables(void) = default;

	// full constructors
	MaterialStateVariables(const MaterialStateVariables&) = default;

	// destructor
	~MaterialStateVariables(void);

	// operator overloading
	MaterialStateVariables& operator= (const MaterialStateVariables&) = default;

	// pack and unpack state variables (vectorize) for message passing
	Vector& pack(void);
	void unpack(Vector& data);

	void printStats(bool detail);

private:
	// un/pack methods
	int getSize(void);
	void vectorvector(int& N, Vector& A, Vector& B, bool forward = true);
	void vectormatrix(int& N, Matrix& A, Vector& B, bool forward = true);
};
#endif