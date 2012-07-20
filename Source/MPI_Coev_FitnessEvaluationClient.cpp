/*
 *  MPI_Coev_FitnessEvaluationClient.cpp
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
 *  This file was created by Jean-Francois Dupuis on 15/01/10.
 */

#include "MPI_Coev_FitnessEvaluationClient.hpp"

#include <beagle/System.hpp>
#include <beagle/Context.hpp>
#include <mpi.h>
#include "CommunicationMPI.h"

#include "beagle/FitnessSimple.hpp"

using namespace Beagle;
using namespace std;

void Beagle::MPI::Coev::FitnessEvaluationClient::operate(Allocator::Handle inContextAllocator, 
														 vector<Genotype::Alloc::Handle>& inGenotypeAlloc,  
														 Fitness::Alloc::Handle inFitnessAlloc) {
	try {
		Vivarium::Handle lVivarium = new Vivarium(inGenotypeAlloc[0], inFitnessAlloc);
		lVivarium->resize(1);

		Context::Handle lEvolContext = castObjectT<Context*>(inContextAllocator->allocate());
		lEvolContext->setSystemHandle(mSystem);
		lEvolContext->setVivariumHandle(lVivarium);
		lEvolContext->setDemeHandle((*lVivarium)[0]);
		lEvolContext->setDemeIndex(0);
		lEvolContext->setGeneration(0);
		lEvolContext->setContinueFlag(true);
		



		int lMessageSize;
		int lNbIndividuals = 0;
		MPI_Status lStatus;
		int lSource;
		
		bool lDone = false;
		while(!lDone) {
			//Receive an individual to evaluate
			MPI_Recv(&lNbIndividuals, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
			lSource = lStatus.MPI_SOURCE;
			if(lStatus.MPI_TAG == eEvolutionEnd) {
				Beagle_LogDetailedM(
									lEvolContext->getSystem().getLogger(),
									"evaluation", "Beagle::MPIEvaluationOp",
									std::string("End of evolution received from process ") + int2str(lSource)
									);
				lDone = true;
			} else {
				Individual::Bag lIndividuals;
				for(unsigned int i = 0; i < lNbIndividuals; ++i) {
					MPI_Recv(&lMessageSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
					char *lMessage = new char[lMessageSize];
					MPI_Recv(lMessage, lMessageSize, MPI_CHAR, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
					
					//Parse the received individual
					std::istringstream lStreamIn(lMessage);
					PACC::XML::Document lXMLParser;
					lXMLParser.parse(lStreamIn);
					PACC::XML::ConstIterator lIndividualRootNode = lXMLParser.getFirstRoot(); 
					
					//Read the received individual
					lEvolContext->getDeme().resize(0);
					Individual::Handle lIndividual = new Individual(inGenotypeAlloc[i]);
					lIndividual->readWithContext(lIndividualRootNode,*lEvolContext);
					
					//Free message std::string
					delete [] lMessage;
					
					lIndividuals.push_back(lIndividual);
				}
				unsigned int lGeneration;
				MPI_Recv(&lGeneration, 1, MPI_INT, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
				lEvolContext->setGeneration(lGeneration);
				
				Beagle_LogTraceM(
								 lEvolContext->getSystem().getLogger(),
								 "evaluation", "Beagle::MPIEvaluationOp",
								 std::string("Evaluating individual send from process ") + int2str(lSource)
								 );
				
				
				//Evaluated the fitness of the received individual
				Fitness::Handle lFitness =  evaluate(lIndividuals,*lEvolContext);
				
				//Send back the fitness
				std::ostringstream lStreamOut;
				PACC::XML::Streamer lXMLStream(lStreamOut);
				
				lFitness->write(lXMLStream);
				//std::cout << "Sending fitness of size " << lStreamOut.str().size()+1 << ":" << std::endl << lStreamOut.str() << std::endl;
				lMessageSize = lStreamOut.str().size()+1;
				
				FitnessSimple::Handle lLogFitness = castHandleT<FitnessSimple>(lFitness);
				
				Beagle_LogTraceM(
								 lEvolContext->getSystem().getLogger(),
								 "evaluation", "Beagle::MPIEvaluationOp",
								 std::string("Sending back fitness of value ")+dbl2str(lLogFitness->getValue())
								 );
				
				MPI_Send(&lMessageSize, 1, MPI_INT, lSource, eMessageSize, MPI_COMM_WORLD);
				MPI_Send(const_cast<char*>(lStreamOut.str().data()), lMessageSize, MPI_CHAR, lSource, eFitness, MPI_COMM_WORLD);
			}
		}
	} catch(Exception& inException) {
		std::cerr << "Exception catched in evaluator:" << std::endl << std::flush;
		std::cerr << inException.what() << std::endl << std::flush;
		exit(1);
	}
	catch(std::exception& inException) {
		std::cerr << "Standard exception catched in evaluator:" << std::endl << std::flush;
		std::cerr << inException.what() << std::endl << std::flush;
		exit(1);
	}
}

void Beagle::MPI::Coev::FitnessEvaluationClient::initialize(int& ioArgc, char** ioArgv) {
	mSystem = new System;

	// Get rank
	int lRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &lRank);	
	
	//Calling user defined initialization
	init();
	
	// Initialize the system.
	mSystem->initialize(ioArgc, ioArgv);
	
	//Make individual log filename for each process
	if(mSystem->getRegister().isRegistered("lg.file.name")) {
		String::Handle lName = castHandleT<String>(mSystem->getRegister().getEntry("lg.file.name"));
		
		std::istringstream lNameStream(lName->getWrappedValue());
		PACC::Tokenizer lTokenizer(lNameStream);
		lTokenizer.setDelimiters(".","");
		
		std::string lPrefix = lTokenizer.getNextToken();
		std::string lSuffix = lTokenizer.getNextToken();
		
		lName->setWrappedValue(lPrefix+"-"+int2str(lRank)+"."+lSuffix);
	} 
	
	// Calling system post-initialization hook.
	mSystem->postInit();
	
	//Calling user defined post-initialization.
	postInit();
}

void Beagle::MPI::Coev::FitnessEvaluationClient::initialize(string inConfigFilename) {
	mSystem = new System;
	
	// Get rank
	int lRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &lRank);
	
	//Calling user defined initialization
	init();
	
	// Initialize the system.
	mSystem->initialize(inConfigFilename);
	
	//Make individual log filename for each process
	if(mSystem->getRegister().isRegistered("lg.file.name")) {
		String::Handle lName = castHandleT<String>(mSystem->getRegister().getEntry("lg.file.name"));
		
		std::istringstream lNameStream(lName->getWrappedValue());
		PACC::Tokenizer lTokenizer(lNameStream);
		lTokenizer.setDelimiters(".","");
		
		std::string lPrefix = lTokenizer.getNextToken();
		std::string lSuffix = lTokenizer.getNextToken();
		
		lName->setWrappedValue(lPrefix+"-"+int2str(lRank)+"."+lSuffix);
	} 
	
	// Calling system post-initialization hook.
	mSystem->postInit();
	
	//Calling user defined post-initialization.
	postInit();
}