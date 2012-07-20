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
 *  \file   MaxFctEvalOp.cpp
 *  \brief  Implementation of the class MaxFctEvalOp.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.4 $
 *  $Date: 2009/11/26 14:41:50 $
 */

#include "beagle/GA.hpp"
#include "MaxFctEvalOp.hpp"
#include <cmath>

using namespace Beagle;

/*!
 *  \brief Construct the individual evaluation operator for maximising the function.
 *  \param inEncoding Encoding of the variables (in bits).
 *  \param inName Name of the operator.
 */
MaxFctEvalOp::MaxFctEvalOp(UIntArray inEncoding, std::string inName) : EvaluationOp(inName)
{
  if(inEncoding.size() != 5)
  throw ValidationException("Size of the encoding vector is different from 5!");
  for(unsigned int i=0; i<5; i++) {
    GA::BitString::DecodingKey lKey(-200.0, 200.0, inEncoding[i]);
    mDecodingKeys.push_back(lKey);
  }
}


/*!
 *  \brief Evaluate the fitness of the given individual.
 *  \param inIndividual Current individual to evaluate.
 *  \param ioContext Evolutionary context.
 *  \return Handle to the fitness value of the individual.
 */
Fitness::Handle MaxFctEvalOp::evaluate(Individual& inIndividual, Context& ioContext)
{
  Beagle_AssertM(inIndividual.size() == 1);
  GA::BitString::Handle lBitString = castHandleT<GA::BitString>(inIndividual[0]);
//  std::vector<double> lX;
	Beagle::DoubleArray lX;
//  lBitString->decodeGray(mDecodingKeys, lX);
  lBitString->decode(mDecodingKeys, lX);
  double lU   = 10.0;
  double lSum = 0.0;
  for(unsigned int i=0; i<5; i++) {
    lSum += ((lX[i])*(lX[i])) + (lU*lU);
    lU += lX[i];
  }
  lSum += (lU*lU);
  double lF = 161.8 / lSum;
  return new FitnessSimple(lF);
}
