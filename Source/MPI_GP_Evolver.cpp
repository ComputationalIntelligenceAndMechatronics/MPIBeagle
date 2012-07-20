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
 *  \file   beagle/GP/src/Evolver.cpp
 *  \brief  Source code of class GP::Evolver.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.2 $
 *  $Date: 2009/11/26 09:50:20 $
 */

#include "beagle/GP.hpp"
#include "MPI_GP_Evolver.hpp"

#include <string>

using namespace Beagle;


/*!
 *  \brief Construct a GP Generational evolver.
 */
Beagle::MPI::GP::Evolver::Evolver()
{
	addOperator(new Beagle::GP::InitGrowOp);
	addOperator(new Beagle::GP::InitFullOp);
	addOperator(new Beagle::GP::InitHalfOp);
	addOperator(new Beagle::GP::CrossoverOp);
	addOperator(new Beagle::GP::MutationStandardOp);
	addOperator(new Beagle::GP::MutationShrinkOp);
	addOperator(new Beagle::GP::MutationSwapOp);
	addOperator(new Beagle::GP::MutationSwapSubtreeOp);
	addOperator(new Beagle::GP::InitGrowConstrainedOp);
	addOperator(new Beagle::GP::InitFullConstrainedOp);
	addOperator(new Beagle::GP::InitHalfConstrainedOp);
	addOperator(new Beagle::GP::CrossoverConstrainedOp);
	addOperator(new Beagle::GP::MutationStandardConstrainedOp);
	addOperator(new Beagle::GP::MutationShrinkConstrainedOp);
	addOperator(new Beagle::GP::MutationSwapConstrainedOp);
	addOperator(new Beagle::GP::MutationSwapSubtreeConstrainedOp);
	addOperator(new Beagle::GP::StatsCalcFitnessKozaOp);
	addOperator(new Beagle::GP::StatsCalcFitnessSimpleOp);
	addOperator(new Beagle::GP::StatsCalcFitnessSimpleOp("GP-StatsCalcFitnessSimpleMinOp"));
	addOperator(new Beagle::GP::PrimitiveUsageStatsOp);
	addOperator(new Beagle::GP::TermMaxHitsOp);
}


/*!
 *  \brief Construct a GP Generational evolver.
 *  \param inEvalOp GP evaluation operator.
 */
Beagle::MPI::GP::Evolver::Evolver(Beagle::EvaluationOp::Handle inEvalOp) 
: Beagle::MPI::Evolver(inEvalOp)
{
	addOperator(inEvalOp);
	addOperator(new Beagle::GP::InitGrowOp);
	addOperator(new Beagle::GP::InitFullOp);
	addOperator(new Beagle::GP::InitHalfOp);
	addOperator(new Beagle::GP::CrossoverOp);
	addOperator(new Beagle::GP::MutationStandardOp);
	addOperator(new Beagle::GP::MutationShrinkOp);
	addOperator(new Beagle::GP::MutationSwapOp);
	addOperator(new Beagle::GP::MutationSwapSubtreeOp);
	addOperator(new Beagle::GP::InitGrowConstrainedOp);
	addOperator(new Beagle::GP::InitFullConstrainedOp);
	addOperator(new Beagle::GP::InitHalfConstrainedOp);
	addOperator(new Beagle::GP::CrossoverConstrainedOp);
	addOperator(new Beagle::GP::MutationStandardConstrainedOp);
	addOperator(new Beagle::GP::MutationShrinkConstrainedOp);
	addOperator(new Beagle::GP::MutationSwapConstrainedOp);
	addOperator(new Beagle::GP::MutationSwapSubtreeConstrainedOp);
	addOperator(new Beagle::GP::StatsCalcFitnessKozaOp);
	addOperator(new Beagle::GP::StatsCalcFitnessSimpleOp);
	addOperator(new Beagle::GP::PrimitiveUsageStatsOp);
	addOperator(new Beagle::GP::TermMaxHitsOp);
	
//	addBootStrapOp("IfThenElseOp");
//	IfThenElseOp::Handle lITE = castHandleT<IfThenElseOp>(getBootStrapSet().back());
//	lITE->setConditionTag("ms.restart.file");
//	lITE->setConditionValue("");
//	lITE->insertPositiveOp("GP-InitHalfOp", getOperatorMap());
//	lITE->insertPositiveOp(inEvalOp->getName(), getOperatorMap());
//	lITE->insertPositiveOp("GP-StatsCalcFitnessSimpleOp", getOperatorMap());
//	lITE->insertNegativeOp("MilestoneReadOp", getOperatorMap());
//	addBootStrapOp("TermMaxGenOp");
//	addBootStrapOp("MilestoneWriteOp");
//	
//	addMainLoopOp("SelectTournamentOp");
//	addMainLoopOp("GP-CrossoverOp");
//	addMainLoopOp("GP-MutationStandardOp");
//	addMainLoopOp("GP-MutationShrinkOp");
//	addMainLoopOp("GP-MutationSwapOp");
//	addMainLoopOp("GP-MutationSwapSubtreeOp");
//	addMainLoopOp(inEvalOp->getName());
//	addMainLoopOp("MigrationRandomRingOp");
//	addMainLoopOp("GP-StatsCalcFitnessSimpleOp");
//	addMainLoopOp("TermMaxGenOp");
//	addMainLoopOp("MilestoneWriteOp");
}

