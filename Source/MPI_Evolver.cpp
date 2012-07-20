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
 *  \file   beagle/MPI/src/Beagle::MPI::Evolver.cpp
 *  \brief  Source code of class MPI::Beagle::MPI::Evolver.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.8 $
 *  $Date: 2010/01/19 13:28:36 $
 */

#include "beagle/Beagle.hpp"
#include <Util/Tokenizer.hpp>

#include "MPI_Evolver.hpp"
#include "mpi.h"
#include "CommunicationMPI.h"

#include <set>

#ifdef BEAGLE_HAVE_LIBZ
#include "gzstream.h"
#endif // BEAGLE_HAVE_LIBZ

extern int gArgc; 
extern char** gArgv;

using namespace Beagle;


/*!
 *  \brief Construct an evolver by adding common operators in it.
 */
Beagle::MPI::Evolver::Evolver()
{
	addOperator(new IfThenElseOp);
	addOperator(new MigrationRandomRingOp);
	addOperator(new MilestoneReadOp);
	addOperator(new MilestoneWriteOp);
	addOperator(new RegisterReadOp);
	addOperator(new SelectRandomOp);
	addOperator(new SelectRouletteOp);
	addOperator(new SelectTournamentOp);
	addOperator(new SelectParsimonyTournOp);
	addOperator(new StatsCalcFitnessSimpleOp);
	addOperator(new StatsCalcFitnessSimpleOp("StatsCalcFitnessSimpleMinOp"));
	addOperator(new StatsCalcFitnessMultiObjOp);
	addOperator(new StatsCalcFitnessMultiObjOp("StatsCalcFitnessMultiObjMinOp"));
	addOperator(new TermMaxGenOp);
	addOperator(new TermMaxFitnessOp);
	addOperator(new TermMinFitnessOp);
	addOperator(new TermMaxEvalsOp);
	addOperator(new GenerationalOp);
	addOperator(new SteadyStateOp);
	addOperator(new MuCommaLambdaOp);
	addOperator(new MuCommaLambdaOp("ec.mulambda.ratio","MuCommaLambdaOp-2"));
	addOperator(new MuPlusLambdaOp);
	addOperator(new NSGA2Op);
	addOperator(new NPGA2Op);
	addOperator(new ParetoFrontCalculateOp);
	addOperator(new DecimateOp);
	addOperator(new OversizeOp);
	addOperator(new RandomShuffleDemeOp);
	addOperator(new HierarchicalFairCompetitionOp);
}

Beagle::MPI::Evolver::Evolver(Beagle::EvaluationOp::Handle inEvalOp) {
	mEvaluator = inEvalOp;
	addOperator(new IfThenElseOp);
	addOperator(new MigrationRandomRingOp);
	addOperator(new MilestoneReadOp);
	addOperator(new MilestoneWriteOp);
	addOperator(new RegisterReadOp);
	addOperator(new SelectRandomOp);
	addOperator(new SelectRouletteOp);
	addOperator(new SelectTournamentOp);
	addOperator(new SelectParsimonyTournOp);
	addOperator(new StatsCalcFitnessSimpleOp);
	addOperator(new StatsCalcFitnessSimpleOp("StatsCalcFitnessSimpleMinOp"));
	addOperator(new StatsCalcFitnessMultiObjOp);
	addOperator(new StatsCalcFitnessMultiObjOp("StatsCalcFitnessMultiObjMinOp"));
	addOperator(new TermMaxGenOp);
	addOperator(new TermMaxFitnessOp);
	addOperator(new TermMinFitnessOp);
	addOperator(new TermMaxEvalsOp);
	addOperator(new GenerationalOp);
	addOperator(new SteadyStateOp);
	addOperator(new MuCommaLambdaOp);
	addOperator(new MuCommaLambdaOp("ec.mulambda.ratio","MuCommaLambdaOp-2"));
	addOperator(new MuPlusLambdaOp);
	addOperator(new NSGA2Op);
	addOperator(new NPGA2Op);
	addOperator(new ParetoFrontCalculateOp);
	addOperator(new DecimateOp);
	addOperator(new OversizeOp);
	addOperator(new RandomShuffleDemeOp);
	addOperator(new HierarchicalFairCompetitionOp);
}

Beagle::MPI::Evolver::~Evolver() { 
	//MPI_Finalize(); 
}

/*!
 *  \brief Evolve a given vivarium.
 *  \param ioVivarium Handle to the vivarium to evolve.
 */
void Beagle::MPI::Evolver::evolve(Vivarium::Handle ioVivarium)
{
//	if(mProcessSize->getWrappedValue() == 1) {
//		Beagle::Evolver::evolve(ioVivarium);
//	}
	
	if(mRank == 0) {
		// We are the evolver
		evolver(ioVivarium);
	}
	else {
		//We are a fitness evaluation client
		evaluater(ioVivarium);
	}
}

void Beagle::MPI::Evolver::evolver(Vivarium::Handle ioVivarium) {
	// Log welcome messages.
	Beagle_LogBasicM(
					 mSystemHandle->getLogger(),
					 "evolver", "Beagle::MPI::Evolver",
					 "Starting an evolution"
					 );
	
	mSystemHandle->getLogger().logCurrentTime(Logger::eBasic);
	
	Beagle_LogObjectM(
					  mSystemHandle->getLogger(),
					  Logger::eDetailed,
					  "evolver", "Beagle::MPI::Evolver",
					  mSystemHandle->getRegister()
					  );
	
	Beagle_LogObjectM(
					  mSystemHandle->getLogger(),
					  Logger::eDetailed,
					  "evolver", "Beagle::MPI::Evolver",
					  (*this)
					  );
	
	// Initialize the evolution context.
	Context::Handle lEvolContext =
	castObjectT<Context*>(mSystemHandle->getContextAllocator().allocate());
	lEvolContext->setSystemHandle(mSystemHandle);
	lEvolContext->setEvolverHandle(this);
	lEvolContext->setVivariumHandle(ioVivarium);
	lEvolContext->setDemeIndex(0);
	lEvolContext->setGeneration(0);
	lEvolContext->setContinueFlag(true);
	
	Beagle_LogTraceM(
					 mSystemHandle->getLogger(),
					 "evolver", "Beagle::MPI::Evolver",
					 std::string("Vivarium resized to ")+uint2str(mPopSize->size())+" demes"
					 );
	ioVivarium->resize(mPopSize->size());
	
	while( lEvolContext->getContinueFlag() ) {
		unsigned int lGeneration = lEvolContext->getGeneration();
		
		if(lGeneration == 0) {
			for(unsigned int i=lEvolContext->getDemeIndex(); i<ioVivarium->size(); i++) {
				lEvolContext->setDemeIndex(i);
				lEvolContext->setDemeHandle((*ioVivarium)[i]);
				
				Beagle_LogInfoM(
								mSystemHandle->getLogger(),
								"evolver", "Beagle::MPI::Evolver",
								std::string("Applying bootstrap operators to the ")+uint2ordinal(i+1)+
								std::string(" deme")
								);
				for(unsigned int j=0; j<mBootStrapSet.size(); j++) {
					Beagle_LogDetailedM(
										mSystemHandle->getLogger(),
										"evolver", "Beagle::MPI::Evolver",
										std::string("Applying \"")+mBootStrapSet[j]->getName()+std::string("\"")
										);
					mBootStrapSet[j]->operate(*(*ioVivarium)[i], *lEvolContext);
				}
				if(lEvolContext->getContinueFlag() == false) break;
				if(i != lEvolContext->getDemeIndex()) break;
				if(lGeneration != lEvolContext->getGeneration()) break;
				if(i == (ioVivarium->size()-1)) {
					lEvolContext->setGeneration(lGeneration+1);
					lEvolContext->setDemeIndex(0);
				}
			}
		}
		else {
			for(unsigned int i=lEvolContext->getDemeIndex(); i<ioVivarium->size(); i++) {
				lEvolContext->setDemeIndex(i);
				lEvolContext->setDemeHandle((*ioVivarium)[i]);
				Beagle_LogInfoM(
								mSystemHandle->getLogger(),
								"evolver", "Beagle::MPI::Evolver",
								std::string("Applying main-loop operators to the ")+uint2ordinal(i+1)+
								std::string(" deme")
								);
				for(unsigned int j=0; j<mMainLoopSet.size(); j++) {
					Beagle_LogDetailedM(
										mSystemHandle->getLogger(),
										"evolver", "Beagle::MPI::Evolver",
										std::string("Applying \"")+mMainLoopSet[j]->getName()+std::string("\"")
										);
					mMainLoopSet[j]->operate(*(*ioVivarium)[i], *lEvolContext);
				}
				if(lEvolContext->getContinueFlag() == false) break;
				if(i != lEvolContext->getDemeIndex()) break;
				if(lGeneration != lEvolContext->getGeneration()) break;
				if(i == (ioVivarium->size()-1)) {
					lEvolContext->setGeneration(lGeneration+1);
					lEvolContext->setDemeIndex(0);
				}
			}
		}
	}
	
	stopEvaluater();
	
	mSystemHandle->getLogger().logCurrentTime(Logger::eBasic);
	Beagle_LogBasicM(
					 mSystemHandle->getLogger(),
					 "evolver", "Beagle::MPI::Evolver",
					 "End of evolution"
					 );
}

void Beagle::MPI::Evolver::evaluater(Vivarium::Handle ioVivarium) {
	Context::Handle lEvolContext = castObjectT<Context*>(mSystemHandle->getContextAllocator().allocate());
	lEvolContext->setSystemHandle(mSystemHandle);
	lEvolContext->setEvolverHandle(this);
	lEvolContext->setVivariumHandle(ioVivarium);
	lEvolContext->setDemeIndex(0);
	lEvolContext->setGeneration(0);
	lEvolContext->setContinueFlag(true);
	ioVivarium->resize(1);
	
	int i = 0;
	lEvolContext->setDemeIndex(i);
	lEvolContext->setDemeHandle((*ioVivarium)[i]);
	
	
	//Start evaluation operator with dummy argument. They will be discarded
	//when looking at the rank
	mEvaluator->operate(*(*ioVivarium)[i], *lEvolContext);
}

void Beagle::MPI::Evolver::stopEvaluater() {
	Beagle_LogBasicM(
					 mSystemHandle->getLogger(),
					 "evolver", "Beagle::MPI::Evolver",
					 "Stopping the evaluators"
					 );
	for(unsigned int i = 1; i < mProcessSize->getWrappedValue(); ++i) {
		MPI_Send(NULL, 0, MPI_CHAR, i, eEvolutionEnd, MPI_COMM_WORLD);
	}
}

/*!
 *  \brief Initialize the evolver, its operators and the system.
 *  \param ioSystem Handle to the system of the evolution.
 *  \param ioArgc Number of elements on the command-line.
 *  \param ioArgv Element on the command-line.
 */
void Beagle::MPI::Evolver::initialize(System::Handle ioSystem, int& ioArgc, char** ioArgv)
{
	// Initialize MPI
	MPI_Init(&ioArgc, &ioArgv);

	// Get rank
	MPI_Comm_rank(MPI_COMM_WORLD, &mRank);
	
	// Get number of process
	int lSize;
	MPI_Comm_size(MPI_COMM_WORLD, &lSize);
	mProcessSize = new Int(lSize);
	if(ioSystem->getRegister().isRegistered("ec.mpi.size")) {
		ioSystem->getRegister().modifyEntry("ec.mpi.size", mProcessSize);
	} else {
		Register::Description lDescription(
										   "MPI Process Size",
										   "Int",
										   "\"\"",
										   "Specify the number of concurent process used to evaluate individuals"
										   );
		ioSystem->getRegister().addEntry("ec.mpi.size", mProcessSize, lDescription);
	}
	
	// Get system handle.
	mSystemHandle = ioSystem;
	
	// Parse command-line.
	parseCommandLine(*ioSystem, ioArgc, ioArgv);
	
	// Logging message.
	Beagle_LogDetailedM(
						ioSystem->getLogger(),
						"evolver", "Beagle::MPI::Evolver",
						"Initializing evolver"
						);
	
	// Add configuration dumper parameter.
	if(ioSystem->getRegister().isRegistered("ec.conf.dump")) {
		mConfigDumper =
		castHandleT<ConfigurationDumper>(ioSystem->getRegister().getEntry("ec.conf.dump"));
	}
	else {
		mConfigDumper = new ConfigurationDumper(*ioSystem, *this, "");
		std::string lLongDescripDump = "Filename used to dump the configuration. ";
		lLongDescripDump += "A configuration dump means that a configuration file is ";
		lLongDescripDump += "written with the evolver (including the composing operators) ";
		lLongDescripDump += "and the register (including the registered parameters and their ";
		lLongDescripDump += "default values). No evolution is conducted on a configuration dump. ";
		lLongDescripDump += "An empty string means no dump.";
		Register::Description lFileDumperDescription(
													 "Configuration dump filename",
													 "String",
													 "\"\"",
													 lLongDescripDump
													 );
		ioSystem->getRegister().addEntry("ec.conf.dump", mConfigDumper, lFileDumperDescription);
	}
	
	// Add configuration file name.
	if(ioSystem->getRegister().isRegistered("ec.conf.file")) {
		mFileName = castHandleT<String>(ioSystem->getRegister().getEntry("ec.conf.file"));
	}
	else {
		mFileName = new String("");
		std::string lLongDescripFN = "The name of a configuration file containing ";
		lLongDescripFN += "evolver and parameter values. A typical configuration file can ";
		lLongDescripFN += "be created with parameter \"ec.conf.dump\".";
		Register::Description lFileNameDescription(
												   "Configuration filename",
												   "String",
												   "\"\"",
												   lLongDescripFN
												   );
		ioSystem->getRegister().addEntry("ec.conf.file", mFileName, lFileNameDescription);
	}
	
	// Add population size parameter
	if(ioSystem->getRegister().isRegistered("ec.pop.size")) {
#ifdef BEAGLE_4
		mPopSize = castHandleT<IntegerVector>(ioSystem->getRegister().getEntry("ec.pop.size"));
	} else {
		mPopSize = new IntegerVector(1,100);
#else
		mPopSize = castHandleT<UIntArray>(ioSystem->getRegister().getEntry("ec.pop.size"));
	} else {
		mPopSize = new UIntArray(1,100);
#endif
		std::string lLongDescrip("Number of demes and size of each deme of the population. ");
		lLongDescrip += "The format of an IntegerVector is S1/S2/.../Sn, where Si is the ith value. ";
		lLongDescrip += "The size of the IntegerVector is the number of demes present in the ";
		lLongDescrip += "vivarium, while each value of the vector is the size of the corresponding ";
		lLongDescrip += "deme.";
		Register::Description lDescription(
										   "Vivarium and demes sizes",
										   "IntegerVector",
										   "100",
										   lLongDescrip
										   );
		ioSystem->getRegister().addEntry("ec.pop.size", mPopSize, lDescription);
	}

	// Initialize operators.
	initializeOperators(*ioSystem);

	// Initialize the system.
	ioSystem->initialize(ioArgc, ioArgv);
	
	//Make individual log filename for each process
	if(ioSystem->getRegister().isRegistered("lg.file.name")) {
		String::Handle lName = castHandleT<String>(ioSystem->getRegister().getEntry("lg.file.name"));
		
		std::istringstream lNameStream(lName->getWrappedValue());
#ifdef BEAGLE_4
		Tokenizer lTokenizer(lNameStream);
		lTokenizer.setWhiteSpace(".");
#else
		PACC::Tokenizer lTokenizer(lNameStream);
		lTokenizer.setDelimiters(".","");
#endif
		
		
		std::string lPrefix = lTokenizer.getNextToken();
		std::string lSuffix = lTokenizer.getNextToken();
		
		lName->setWrappedValue(lPrefix+"-"+int2str(mRank)+"."+lSuffix);
	} 
	
	// Calling system post-initialization hook.
	ioSystem->postInit();
	
	// Calling post-initialization hook of operators.
	postInitOperators(*ioSystem);
}

/*!
 *  \brief Initialize the evolver, its operators and the system.
 *  \param ioSystem Handle to the system of the evolution.
 *  \param inConfigFilename Filename containing configuration values.
 */
void Beagle::MPI::Evolver::initialize(System::Handle ioSystem, Beagle::string inConfigFilename)
{
	Beagle_StackTraceBeginM();

	// Get rank
	MPI_Comm_rank(MPI_COMM_WORLD, &mRank);

	// Get number of process
	int lSize;
	MPI_Comm_size(MPI_COMM_WORLD, &lSize);
	mProcessSize = new Int(lSize);
	if(ioSystem->getRegister().isRegistered("ec.mpi.size")) {
		ioSystem->getRegister().modifyEntry("ec.mpi.size", mProcessSize);
	} else {
		Register::Description lDescription(
										   "MPI Process Size",
										   "Int",
										   "\"\"",
										   "Specify the number of concurent process used to evaluate individuals"
										   );
		ioSystem->getRegister().addEntry("ec.mpi.size", mProcessSize, lDescription);
	}
	
	// Get system handle.
	mSystemHandle = ioSystem;
	
	// Reading evolver configuration, if any.
	if(inConfigFilename.empty() == false) readEvolverFile(inConfigFilename);
	
	// Logging message.
	Beagle_LogDetailedM(
						ioSystem->getLogger(),
						"evolver", "Beagle::Evolver",
						"Initializing evolver"
						);
	
	// Add configuration dumper parameter.
	if(ioSystem->getRegister().isRegistered("ec.conf.dump")) {
		mConfigDumper =
		castHandleT<ConfigurationDumper>(ioSystem->getRegister().getEntry("ec.conf.dump"));
	}
	else {
		mConfigDumper = new ConfigurationDumper(*ioSystem, *this, "");
		string lLongDescripDump = "Filename used to dump the configuration. ";
		lLongDescripDump += "A configuration dump means that a configuration file is ";
		lLongDescripDump += "written with the evolver (including the composing operators) ";
		lLongDescripDump += "and the register (including the registered parameters and their ";
		lLongDescripDump += "default values). No evolution is conducted on a configuration dump. ";
		lLongDescripDump += "An empty string means no dump.";
		Register::Description lFileDumperDescription(
													 "Configuration dump filename",
													 "String",
													 "\"\"",
													 lLongDescripDump
													 );
		ioSystem->getRegister().addEntry("ec.conf.dump", mConfigDumper, lFileDumperDescription);
	}
	
	// Add configuration file name.
	if(ioSystem->getRegister().isRegistered("ec.conf.file")) {
		mFileName = castHandleT<String>(ioSystem->getRegister().getEntry("ec.conf.file"));
	}
	else {
		mFileName = new String(inConfigFilename);
		string lDefaultFileName = string("\"") + inConfigFilename + string("\"");
		string lLongDescripFN = "The name of a configuration file containing ";
		lLongDescripFN += "evolver and parameter values. A typical configuration file can ";
		lLongDescripFN += "be created with parameter \"ec.conf.dump\".";
		Register::Description lFileNameDescription(
												   "Configuration filename",
												   "String",
												   lDefaultFileName,
												   lLongDescripFN
												   );
		ioSystem->getRegister().addEntry("ec.conf.file", mFileName, lFileNameDescription);
	}
	
	// Add population size parameter
	if(ioSystem->getRegister().isRegistered("ec.pop.size")) {
		mPopSize = castHandleT<UIntArray>(ioSystem->getRegister().getEntry("ec.pop.size"));
	} else {
		mPopSize = new UIntArray(1,100);
		string lLongDescrip("Number of demes and size of each deme of the population. ");
		lLongDescrip += "The format of an UIntArray is S1,S2,...,Sn, where Si is the ith value. ";
		lLongDescrip += "The size of the UIntArray is the number of demes present in the ";
		lLongDescrip += "vivarium, while each value of the vector is the size of the corresponding ";
		lLongDescrip += "deme.";
		Register::Description lDescription(
										   "Vivarium and demes sizes",
										   "UIntArray",
										   "100",
										   lLongDescrip
										   );
		ioSystem->getRegister().addEntry("ec.pop.size", mPopSize, lDescription);
	}
	
	// Initialize operators.
	initializeOperators(*ioSystem);
	
	// Initialize the system.
	ioSystem->initialize(inConfigFilename);
	
	//Make individual log filename for each process
	if(ioSystem->getRegister().isRegistered("lg.file.name")) {
		String::Handle lName = castHandleT<String>(ioSystem->getRegister().getEntry("lg.file.name"));
		
		std::istringstream lNameStream(lName->getWrappedValue());
#ifdef BEAGLE_4
		Tokenizer lTokenizer(lNameStream);
		lTokenizer.setWhiteSpace(".");
#else
		PACC::Tokenizer lTokenizer(lNameStream);
		lTokenizer.setDelimiters(".","");
#endif
		
		
		std::string lPrefix = lTokenizer.getNextToken();
		std::string lSuffix = lTokenizer.getNextToken();
		
		lName->setWrappedValue(lPrefix+"-"+int2str(mRank)+"."+lSuffix);
	} 
	
	// Calling system post-initialization hook.
	ioSystem->postInit();
	
	// Calling post-initialization hook of operators.
	postInitOperators(*ioSystem);
	
	Beagle_StackTraceEndM("void Evolver::initialize(System::Handle, string)");
}
