/*
 *  MPI_GA_EvolverFloatVector.cpp
 *  Copyright 2010 Jean-Francois Dupuis.
 *
 *  This file is part of MPIBeagle.
 *  
 *  MPIBeagle is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MPIBeagle is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MPIBeagle.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  This file was created by Jean-Francois Dupuis on 14/01/10.
 */

#include "MPI_GA_EvolverFloatVector.hpp"
#include "beagle/GA.hpp"

using namespace Beagle;

/*!
 *  \brief Construct a real-valued GA evolver.
 *  \param inInitSize Size of the GA float vectors.
 */
Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(unsigned int inInitSize)
{
	Beagle_StackTraceBeginM();
	addOperator(new Beagle::GA::InitFltVecOp(inInitSize));
	addOperator(new Beagle::GA::InitCMAFltVecOp(inInitSize));
	addOperator(new Beagle::GA::CrossoverBlendFltVecOp);
	addOperator(new Beagle::GA::CrossoverSBXFltVecOp);
	addOperator(new Beagle::GA::CrossoverOnePointFltVecOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsFltVecOp);
	addOperator(new Beagle::GA::CrossoverUniformFltVecOp);
	addOperator(new Beagle::GA::MutationGaussianFltVecOp);
	addOperator(new Beagle::GA::MutationCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp("ga.cmaes.mulambdaratio",
												  "GA-MuWCommaLambdaCMAFltVecOp-2"));
	addOperator(new Beagle::GA::TermCMAOp);
	Beagle_StackTraceEndM("Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(unsigned int inInitSize)");
}


/*!
 *  \brief Construct a real-valued GA evolver.
 *  \param inEvalOp Evaluation operator.
 *  \param inInitSize Size of the GA float vectors.
 */
Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(EvaluationOp::Handle inEvalOp, unsigned int inInitSize)
{
	Beagle_StackTraceBeginM();
	addOperator(inEvalOp);
	addOperator(new Beagle::GA::InitFltVecOp(inInitSize));
	addOperator(new Beagle::GA::InitCMAFltVecOp(inInitSize));
	addOperator(new Beagle::GA::CrossoverBlendFltVecOp);
	addOperator(new Beagle::GA::CrossoverSBXFltVecOp);
	addOperator(new Beagle::GA::CrossoverOnePointFltVecOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsFltVecOp);
	addOperator(new Beagle::GA::CrossoverUniformFltVecOp);
	addOperator(new Beagle::GA::MutationGaussianFltVecOp);
	addOperator(new Beagle::GA::MutationCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp("ga.cmaes.mulambdaratio",
												  "GA-MuWCommaLambdaCMAFltVecOp-2"));
	addOperator(new Beagle::GA::TermCMAOp);
	
	addBootStrapOp("IfThenElseOp");
	IfThenElseOp::Handle lITE = castHandleT<IfThenElseOp>(getBootStrapSet().back());
	lITE->setConditionTag("ms.restart.file");
	lITE->setConditionValue("");
	lITE->insertPositiveOp("GA-InitFltVecOp", getOperatorMap());
	lITE->insertPositiveOp(inEvalOp->getName(), getOperatorMap());
	lITE->insertPositiveOp("StatsCalcFitnessSimpleOp", getOperatorMap());
	lITE->insertNegativeOp("MilestoneReadOp", getOperatorMap());
	addBootStrapOp("TermMaxGenOp");
	addBootStrapOp("MilestoneWriteOp");
	
	addMainLoopOp("SelectTournamentOp");
	addMainLoopOp("GA-CrossoverBlendFltVecOp");
	addMainLoopOp("GA-MutationGaussianFltVecOp");
	addMainLoopOp(inEvalOp->getName());
	addMainLoopOp("MigrationRandomRingOp");
	addMainLoopOp("StatsCalcFitnessSimpleOp");
	addMainLoopOp("TermMaxGenOp");
	addMainLoopOp("MilestoneWriteOp");
	Beagle_StackTraceEndM("Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(EvaluationOp::Handle inEvalOp, unsigned int inInitSize)");
}


/*!
 *  \brief Construct a real-valued GA evolver.
 *  \param inInitSize Size of the GA float vectors.
 *  \deprecated Use EvolverFloatVector(unsigned int) constructor instead.
 *  \throw Beagle::RunTimeException If init size vector has more than one value.
 */
Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(UIntArray inInitSize)
{
	Beagle_StackTraceBeginM();
	if(inInitSize.size()==0) {
		addOperator(new Beagle::GA::InitFltVecOp(0));
		addOperator(new Beagle::GA::InitCMAFltVecOp(0));
	}
	else if(inInitSize.size()==1) {
		addOperator(new Beagle::GA::InitFltVecOp(inInitSize[0]));
		addOperator(new Beagle::GA::InitCMAFltVecOp(inInitSize[0]));
	}
	else {
		std::ostringstream lOSS;
		lOSS << "Initialization of float vector individuals with more than one float vector ";
		lOSS << "is no more valid. You should use individuals made of one float vector, or ";
		lOSS << "define your own float vector initialization operator.";
		throw Beagle_RunTimeExceptionM(lOSS.str().c_str());
	}
	addOperator(new Beagle::GA::CrossoverBlendFltVecOp);
	addOperator(new Beagle::GA::CrossoverSBXFltVecOp);
	addOperator(new Beagle::GA::CrossoverOnePointFltVecOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsFltVecOp);
	addOperator(new Beagle::GA::CrossoverUniformFltVecOp);
	addOperator(new Beagle::GA::MutationGaussianFltVecOp);
	addOperator(new Beagle::GA::MutationCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp("ga.cmaes.mulambdaratio",
												  "GA-MuWCommaLambdaCMAFltVecOp-2"));
	addOperator(new Beagle::GA::TermCMAOp);
	
	Beagle_StackTraceEndM("Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(UIntArray inInitSize)");
}


/*!
 *  \brief Construct a real-valued GA evolver.
 *  \param inEvalOp Evaluation operator.
 *  \param inInitSize Size of the GA float vectors.
 *  \deprecated Use EvolverFloatVector(EvaluationOp::Handle,unsigned int) constructor instead.
 *  \throw Beagle::RunTimeException If init size vector has more than one value.
 */
Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(EvaluationOp::Handle inEvalOp, UIntArray inInitSize)
{
	Beagle_StackTraceBeginM();
	addOperator(inEvalOp);
	if(inInitSize.size()==0) {
		addOperator(new Beagle::GA::InitFltVecOp(0));
		addOperator(new Beagle::GA::InitCMAFltVecOp(0));
	}
	else if(inInitSize.size()==1) {
		addOperator(new Beagle::GA::InitFltVecOp(inInitSize[0]));
		addOperator(new Beagle::GA::InitCMAFltVecOp(inInitSize[0]));
	}
	else {
		std::ostringstream lOSS;
		lOSS << "Initialization of float vector individuals with more than one float vector ";
		lOSS << "is no more valid. You should use individuals made of one float vector, or ";
		lOSS << "define your own float vector initialization operator.";
		throw Beagle_RunTimeExceptionM(lOSS.str().c_str());
	}
	addOperator(new Beagle::GA::CrossoverBlendFltVecOp);
	addOperator(new Beagle::GA::CrossoverSBXFltVecOp);
	addOperator(new Beagle::GA::CrossoverOnePointFltVecOp);
	addOperator(new Beagle::GA::CrossoverTwoPointsFltVecOp);
	addOperator(new Beagle::GA::CrossoverUniformFltVecOp);
	addOperator(new Beagle::GA::MutationGaussianFltVecOp);
	addOperator(new Beagle::GA::MutationCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp);
	addOperator(new Beagle::GA::MuWCommaLambdaCMAFltVecOp("ga.cmaes.mulambdaratio",
												  "GA-MuWCommaLambdaCMAFltVecOp-2"));
	addOperator(new Beagle::GA::TermCMAOp);
	
	addBootStrapOp("IfThenElseOp");
	IfThenElseOp::Handle lITE = castHandleT<IfThenElseOp>(getBootStrapSet().back());
	lITE->setConditionTag("ms.restart.file");
	lITE->setConditionValue("");
	lITE->insertPositiveOp("GA-InitFltVecOp", getOperatorMap());
	lITE->insertPositiveOp(inEvalOp->getName(), getOperatorMap());
	lITE->insertPositiveOp("StatsCalcFitnessSimpleOp", getOperatorMap());
	lITE->insertNegativeOp("MilestoneReadOp", getOperatorMap());
	addBootStrapOp("TermMaxGenOp");
	addBootStrapOp("MilestoneWriteOp");
	
	addMainLoopOp("SelectTournamentOp");
	addMainLoopOp("GA-CrossoverBlendFltVecOp");
	addMainLoopOp("GA-MutationGaussianFltVecOp");
	addMainLoopOp(inEvalOp->getName());
	addMainLoopOp("MigrationRandomRingOp");
	addMainLoopOp("StatsCalcFitnessSimpleOp");
	addMainLoopOp("TermMaxGenOp");
	addMainLoopOp("MilestoneWriteOp");
	Beagle_StackTraceEndM("Beagle::MPI::GA::EvolverFloatVector::EvolverFloatVector(EvaluationOp::Handle inEvalOp, UIntArray inInitSize)");
}


/*!
 *  \brief Initialize the evolver, its operators and the system.
 *  \param ioSystem Handle to the system of the evolution.
 *  \param ioArgc Number of elements on the command-line.
 *  \param ioArgv Element on the command-line.
 */
void Beagle::MPI::GA::EvolverFloatVector::initialize(Beagle::System::Handle ioSystem, int ioArgc, char** ioArgv)
{
	Beagle_StackTraceBeginM();
	ioSystem->addComponent(new Beagle::GA::CMAHolder);
	Beagle::MPI::Evolver::initialize(ioSystem, ioArgc, ioArgv);
	Beagle_StackTraceEndM("void Beagle::MPI::GA::EvolverFloatVector::initialize(System::Handle ioSystem, int, char**)");
}


/*!
 *  \brief Initialize the evolver, its operators and the system.
 *  \param ioSystem Handle to the system of the evolution.
 *  \param inConfigFilename Configuration file from which system and evolver are read.
 */
void Beagle::MPI::GA::EvolverFloatVector::initialize(Beagle::System::Handle ioSystem, std::string inConfigFilename)
{
	
	Beagle_StackTraceBeginM();
	ioSystem->addComponent(new Beagle::GA::CMAHolder);
	Beagle::MPI::Evolver::initialize(ioSystem,inConfigFilename);
	Beagle_StackTraceEndM("void Beagle::MPI::GA::EvolverFloatVector::initialize(System::Handle ioSystem, std::string)");
}
