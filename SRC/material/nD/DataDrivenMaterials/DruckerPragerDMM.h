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

#ifndef DruckerPragerDMM_h
#define DruckerPragerDMM_h

#include "MultiYieldSurfaceHardeningSoftening.h"

class DruckerPragerDMM : public MultiYieldSurfaceHardeningSoftening
{
public:
	//full constructor
	DruckerPragerDMM(int tag, double rho,
		double Kref, double Gref, double Pref, double modn,
		DataDrivenNestedSurfaces* theData,
		int dataDriverType, int integrationType, bool verbosity);

	//null constructor
	DruckerPragerDMM();

	//destructor
	~DruckerPragerDMM();

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
	double yieldFunction(const Vector& stress, const int num_yield_surface, bool yield_stress);		// Return yield function (F) value
	Vector get_dF_dS(const Vector& stress, const int num_yield_surface);							// Return normal to the yield surface w.r.t stress
	Vector get_dF_dA(const Vector& stress, const int num_yield_surface);							// Return normal to the yield surface w.r.t alpha (backstress)
	Vector get_dH_dA(const Vector& stress, const int num_yield_surface);							// Return normal to the hardening potential w.r.t alpha (backstress)
	Vector get_dP_dS(const Vector& stress, const int num_yield_surface);							// Return normal to the plastic potential w.r.t stress

	// closest point projection methods
	double linearizedFlow(const double dlambda);
	double yieldf(const Vector& stress, const int num_ys, bool yield_stress);
	Vector Qi(const Vector& stress, const int num_ys);
	Vector Di(const Vector& stress, const int num_ys);
};
#endif