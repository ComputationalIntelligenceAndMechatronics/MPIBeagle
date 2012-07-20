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
 *  \file   beagle/GA/src/EvolverBitString.cpp
 *  \brief  Source code of class GA::EvolverBitString.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.3 $
 *  $Date: 2009/11/26 14:41:50 $
 */

#include "beagle/GA.hpp"
#include "MPI_GA_EvolverBitString.hpp"

#include <string>

using namespace Beagle;


/*!
 *  \brief Construct a GA Generational evolver.
 *  \param inInitSize Number of bits in the GA bit strings
 */
Beagle::MPI::GA::EvolverBitString::EvolverBitString(unsigned int inInitSize)
{
  addOperator(new Beagle::GA::InitBitStrOp(inInitSize));
	addOperator(new Beagle::GA::CrossoverOnePointBitStrOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsBitStrOp);
	addOperator(new Beagle::GA::CrossoverUniformBitStrOp);
	addOperator(new Beagle::GA::MutationFlipBitStrOp);
}


/*!
 *  \brief Construct a GA Generational evolver.
 *  \param inEvalOp Evaluation operator.
 *  \param inInitSize Number of bits in the GA bit strings.
 */
Beagle::MPI::GA::EvolverBitString::EvolverBitString(MPI::EvaluationOp::Handle inEvalOp, unsigned int inInitSize)
: Beagle::MPI::Evolver(inEvalOp)
{
  addOperator(inEvalOp);
	addOperator(new Beagle::GA::InitBitStrOp(inInitSize));
	addOperator(new Beagle::GA::CrossoverOnePointBitStrOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsBitStrOp);
	addOperator(new Beagle::GA::CrossoverUniformBitStrOp);
	addOperator(new Beagle::GA::MutationFlipBitStrOp);

  addBootStrapOp("IfThenElseOp");
  IfThenElseOp::Handle lITE = castHandleT<IfThenElseOp>(getBootStrapSet().back());
  lITE->setConditionTag("ms.restart.file");
  lITE->setConditionValue("");
  lITE->insertPositiveOp("GA-InitBitStrOp", getOperatorMap());
  lITE->insertPositiveOp(inEvalOp->getName(), getOperatorMap());
  lITE->insertPositiveOp("StatsCalcFitnessSimpleOp", getOperatorMap());
  lITE->insertNegativeOp("MilestoneReadOp", getOperatorMap());
  addBootStrapOp("TermMaxGenOp");
  addBootStrapOp("MilestoneWriteOp");

  addMainLoopOp("SelectTournamentOp");
  addMainLoopOp("GA-CrossoverOnePointBitStrOp");
  addMainLoopOp("GA-MutationFlipBitStrOp");
  addMainLoopOp(inEvalOp->getName());
  addMainLoopOp("MigrationRandomRingOp");
  addMainLoopOp("StatsCalcFitnessSimpleOp");
  addMainLoopOp("TermMaxGenOp");
  addMainLoopOp("MilestoneWriteOp");
}


/*!
 *  \brief Construct a GA Generational evolver.
 *  \param inInitSize Size of the GA bit strings.
 *  \deprecated Use EvolverBitString(EvaluationOp::Handle,unsigned int) constructor instead.
 *  \throw Beagle::RunTimeException If init size vector has more than one value.
 */

Beagle::MPI::GA::EvolverBitString::EvolverBitString(UIntArray inInitSize)

{
	if(inInitSize.size()==0) addOperator(new Beagle::GA::InitBitStrOp(0));
	else if(inInitSize.size()==1) addOperator(new Beagle::GA::InitBitStrOp(inInitSize[0]));
  else {
    std::ostringstream lOSS;
    lOSS << "Initialization of bit string individuals with more than one bit string ";
    lOSS << "is no more valid. You should use individuals made of one bit string, or ";
    lOSS << "define your own bit string initialization operator.";
    throw Beagle_RunTimeExceptionM(lOSS.str());
  }
	addOperator(new Beagle::GA::CrossoverOnePointBitStrOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsBitStrOp);
	addOperator(new Beagle::GA::CrossoverUniformBitStrOp);
	addOperator(new Beagle::GA::MutationFlipBitStrOp);
}


/*!
 *  \brief Construct a GA Generational evolver.
 *  \param inEvalOp Evaluation operator.
 *  \param inInitSize Size of the GA bit strings.
 *  \deprecated Use EvolverBitString(EvaluationOp::Handle,unsigned int) constructor instead.
 *  \throw Beagle::RunTimeException If init size vector has more than one value.
 */
Beagle::MPI::GA::EvolverBitString::EvolverBitString(EvaluationOp::Handle inEvalOp, UIntArray inInitSize) : Beagle::MPI::Evolver(inEvalOp)
{
  addOperator(inEvalOp);
	if(inInitSize.size()==0) addOperator(new Beagle::GA::InitBitStrOp(0));
	else if(inInitSize.size()==1) addOperator(new Beagle::GA::InitBitStrOp(inInitSize[0]));
  else {
    std::ostringstream lOSS;
    lOSS << "Initialization of bit string individuals with more than one bit string ";
    lOSS << "is no more valid. You should use individuals made of one bit string, or ";
    lOSS << "define your own bit string initialization operator.";
    throw Beagle_RunTimeExceptionM(lOSS.str());
  }
	addOperator(new Beagle::GA::CrossoverOnePointBitStrOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsBitStrOp);
	addOperator(new Beagle::GA::CrossoverUniformBitStrOp);
	addOperator(new Beagle::GA::MutationFlipBitStrOp);

  addBootStrapOp("IfThenElseOp");
  IfThenElseOp::Handle lITE = castHandleT<IfThenElseOp>(getBootStrapSet().back());
  lITE->setConditionTag("ms.restart.file");
  lITE->setConditionValue("");
  lITE->insertPositiveOp("GA-InitBitStrOp", getOperatorMap());
  lITE->insertPositiveOp(inEvalOp->getName(), getOperatorMap());
  lITE->insertPositiveOp("StatsCalcFitnessSimpleOp", getOperatorMap());
  lITE->insertNegativeOp("MilestoneReadOp", getOperatorMap());
  addBootStrapOp("TermMaxGenOp");
  addBootStrapOp("MilestoneWriteOp");

  addMainLoopOp("SelectTournamentOp");
  addMainLoopOp("GA-CrossoverOnePointBitStrOp");
  addMainLoopOp("GA-MutationFlipBitStrOp");
  addMainLoopOp(inEvalOp->getName());
  addMainLoopOp("MigrationRandomRingOp");
  addMainLoopOp("StatsCalcFitnessSimpleOp");
  addMainLoopOp("TermMaxGenOp");
  addMainLoopOp("MilestoneWriteOp");
}

