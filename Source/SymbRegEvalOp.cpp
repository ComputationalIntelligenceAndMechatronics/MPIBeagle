/*
 *  Open BEAGLE
 *  Copyright (C) 2001-2004 by Christian Gagne and Marc Parizeau
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Contact:
 *  Laboratoire de Vision et Systemes Numeriques
 *  Departement de genie electrique et de genie informatique
 *  Universite Laval, Quebec, Canada, G1K 7P4
 *  http://vision.gel.ulaval.ca
 *
 */

/*!
 *  \file   SymbRegEvalOp.cpp
 *  \brief  Implementation of the class SymbRegEvalOp.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.3 $
 *  $Date: 2009/11/26 14:41:50 $
 */

#include "beagle/GP.hpp"
#include "SymbRegEvalOp.hpp"

#include <cmath>

using namespace Beagle;

/*!
 *  \brief Construct a new symbolic regression evaluation operator.
 *  \param inName Name of the evaluation operator.
 */
SymbRegEvalOp::SymbRegEvalOp(std::string inName) :
#ifdef WITHOUT_MPI
	GP::EvaluationOp(inName),
#else
	MPI::GP::EvaluationOp(inName),
#endif
	mX(0), mY(0)
{ }


/*!
 *  \brief Evaluate the individual fitness for the symbolic regression problem.
 *  \param inIndividual Individual to evaluate.
 *  \param ioContext Evolutionary context.
 *  \return Handle to the fitness measure,
 */
Fitness::Handle SymbRegEvalOp::evaluate(GP::Individual& inIndividual, GP::Context& ioContext)
{ 
	double lSquareError = 0.0;
	for(unsigned int i=0; i<mX.size(); i++) {
		setValue("X", mX[i], ioContext);
		Double lResult;
		inIndividual.run(lResult, ioContext);
		double lError = mY[i]-lResult;
		lSquareError += (lError*lError);
	}
	double lMSE  = lSquareError / mX.size();
	double lRMSE = sqrt(lMSE);
	double lFitness = (1.0 / (lRMSE + 1.0));
	return new FitnessSimple(lFitness);
}


/*!
 * \brief Post-initialize the operator by sampling the function to regress.
 * \param ioSystem System to use to sample.
 */
void SymbRegEvalOp::postInit(System& ioSystem)
{
#ifdef WITHOUT_MPI
	GP::EvaluationOp::postInit(ioSystem);
#else
	MPI::GP::EvaluationOp::postInit(ioSystem);
#endif
	
	//  for(unsigned int i=0; i<20; i++) {
	//    mX.push_back(ioSystem.getRandomizer().rollUniform(-1.0,1.0));
	//    mY.push_back(mX[i]*(mX[i]*(mX[i]*(mX[i]+1.0)+1.0)+1.0));
	//  }
	for(unsigned int i=0; i<20; i++) {
		mX.push_back(i);
		mY.push_back(mX[i]*(mX[i]*(mX[i]*(mX[i]+1.0)+1.0)+1.0));
	}
}
