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
 *  \file   beagle/src/MPIEvaluationOp.cpp
 *  \brief  Source code of class MPIEvaluationOp.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.19 $
 *  $Date: 2010/08/23 16:08:01 $
 */

#include "beagle/Beagle.hpp"
#include "MPI_EvaluationOp.hpp"
#include <mpi.h>
#include <algorithm>
#include <string>
#include <sstream>

#include <XML.hpp>

#include <beagle/GA.hpp>
#include "CommunicationMPI.h"
#include "VectorUtil.h"

using namespace Beagle;


/*!
 *  \brief Construct a new evaluation operator.
 *  \param inName Name of the operator.
 */
Beagle::MPI::EvaluationOp::EvaluationOp(std::string inName) :
Beagle::EvaluationOp(inName)
{ }


/*!
 *  \brief Apply the evaluation operation on a breeding pool, returning a evaluated bred individual.
 *  \param inBreedingPool Breeding pool to use for the breeding operation.
 *  \param inChild Node handle associated to child node in the breeder tree.
 *  \param ioContext Evolutionary context of the breeding operation.
 *  \return Evaluated bred individual.
 */
Individual::Handle Beagle::MPI::EvaluationOp::breed(Individual::Bag& inBreedingPool,
													BreederNode::Handle inChild,
													Context& ioContext)
{
	Beagle_NonNullPointerAssertM(inChild);
	
	Deme& lDeme = *ioContext.getDemeHandle();
	if(lDeme.getStats()->isValid()) {
		ioContext.setProcessedDeme(0);
		if((ioContext.getGeneration()!=0) && (lDeme.getStats()->existItem("total-processed"))) {
			ioContext.setTotalProcessedDeme((unsigned int)lDeme.getStats()->getItem("total-processed"));
		}
		else ioContext.setTotalProcessedDeme(0);
		lDeme.getStats()->setInvalid();
		
		if(ioContext.getDemeIndex()==0) {
			Stats& lVivaStats = *ioContext.getVivarium().getStats();
			ioContext.setProcessedVivarium(0);
			if((ioContext.getGeneration()!=0) && (lVivaStats.existItem("total-processed"))) {
				ioContext.setTotalProcessedVivarium((unsigned int)lVivaStats.getItem("total-processed"));
			}
			else ioContext.setTotalProcessedVivarium(0);
			lVivaStats.setInvalid();
		}
	}
	
	Beagle_NonNullPointerAssertM(inChild);
	Beagle_NonNullPointerAssertM(inChild->getBreederOp());
	Individual::Handle lBredIndividual =
    inChild->getBreederOp()->breed(inBreedingPool, inChild->getFirstChild(), ioContext);
    
	if((lBredIndividual->getFitness()==NULL) || (lBredIndividual->getFitness()->isValid()==false)) {
		Beagle_LogVerboseM(
						   ioContext.getSystem().getLogger(),
						   "evaluation", "Beagle::MPIEvaluationOp",
						   "Evaluating the fitness of a new bred individual"
						   );

		individualEvaluation(*lBredIndividual, ioContext);
//		lBredIndividual->setFitness(evaluate(*lBredIndividual, ioContext));
//		lBredIndividual->getFitness()->setValid();
		
		ioContext.setProcessedDeme(ioContext.getProcessedDeme()+1);
		ioContext.setTotalProcessedDeme(ioContext.getTotalProcessedDeme()+1);
		ioContext.setProcessedVivarium(ioContext.getProcessedVivarium()+1);
		ioContext.setTotalProcessedVivarium(ioContext.getTotalProcessedVivarium()+1);
		
		Beagle_LogVerboseM(
						   ioContext.getSystem().getLogger(),
						   "evaluation", "Beagle::MPIEvaluationOp",
						   std::string("The individual fitness value is: ")+
						   lBredIndividual->getFitness()->serialize()
						   );
		
		if(mDemeHOFSize->getWrappedValue() > 0) {
			Beagle_LogVerboseM(
							   ioContext.getSystem().getLogger(),
							   "evaluation", "Beagle::MPIEvaluationOp",
							   "Updating the deme hall-of-fame"
							   );
			lDeme.getHallOfFame().updateWithIndividual(mDemeHOFSize->getWrappedValue(),
													   *lBredIndividual, ioContext);
		}
		if(mVivaHOFSize->getWrappedValue() > 0) {
			Beagle_LogVerboseM(
							   ioContext.getSystem().getLogger(),
							   "evaluation", "Beagle::MPIEvaluationOp",
							   "Updating the vivarium hall-of-fame"
							   );
			ioContext.getVivarium().getHallOfFame().updateWithIndividual(mVivaHOFSize->getWrappedValue(),
																		 *lBredIndividual, ioContext);
		}
	}
	
	return lBredIndividual;
}


/*!
 *  \return Return selection probability of breeder operator.
 *  \param inChild Child node in the breeder tree.
 */
float Beagle::MPI::EvaluationOp::getBreedingProba(BreederNode::Handle inChild)
{
	Beagle_NonNullPointerAssertM(inChild);
	Beagle_NonNullPointerAssertM(inChild->getBreederOp());
	return inChild->getBreederOp()->getBreedingProba(inChild->getFirstChild());
}


/*!
 *  \brief Initialize the parameters evaluation operator.
 *  \param ioSystem System to use to initialize the operator.
 */
void Beagle::MPI::EvaluationOp::initialize(System& ioSystem)
{
	//Get MPI information
	MPI_Comm_rank(MPI_COMM_WORLD, &mRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mProcessSize);
	
	BreederOp::initialize(ioSystem);
	
	if(ioSystem.getRegister().isRegistered("ec.hof.vivasize")) {
		mVivaHOFSize =
		castHandleT<UInt>(ioSystem.getRegister().getEntry("ec.hof.vivasize"));
	} else {
		mVivaHOFSize = new UInt(1);
		std::string lLongDescript = "Number of individuals kept in vivarium's hall-of-fame ";
		lLongDescript += "(best individuals so far). Note that a hall-of-fame contains only ";
		lLongDescript += "copies of the best individuals so far and is not used by the evolution ";
		lLongDescript += "process.";
		Register::Description lDescription(
										   "Vivarium's hall-of-fame size",
										   "UInt",
										   "1",
										   lLongDescript
										   );
		ioSystem.getRegister().addEntry("ec.hof.vivasize", mVivaHOFSize, lDescription);
	}
	
	if(ioSystem.getRegister().isRegistered("ec.hof.demesize")) {
		mDemeHOFSize =
		castHandleT<UInt>(ioSystem.getRegister().getEntry("ec.hof.demesize"));
	} else {
		mDemeHOFSize = new UInt(0);
		std::string lLongDescript = "Number of individuals kept in each deme's hall-of-fame ";
		lLongDescript += "(best individuals so far). Note that a hall-of-fame contains only ";
		lLongDescript += "copies of the best individuals so far and is not used by the evolution ";
		lLongDescript += "process.";
		Register::Description lDescription(
										   "Demes' hall-of-fame size",
										   "UInt",
										   "0",
										   lLongDescript
										   );
		ioSystem.getRegister().addEntry("ec.hof.demesize", mDemeHOFSize, lDescription);
	}
}


/*!
 *  \brief Apply the evaluation process on the invalid individuals of the deme.
 *  \param ioDeme Deme to process.
 *  \param ioContext Context of the evolution.
 */
void Beagle::MPI::EvaluationOp::operate(Deme& ioDeme, Context& ioContext)
{
	if(mProcessSize == 1)
		Beagle::EvaluationOp::operate(ioDeme,ioContext);
	else {
		if(mRank == 0) { 
			evolverOperate(ioDeme, ioContext);
		}
		else {
			evaluatorOperate(ioDeme, ioContext);
		}
	}
}

void Beagle::MPI::EvaluationOp::individualEvaluation(Individual& ioIndividal, Context& ioContext) {
	Fitness::Handle lFitness = evaluate(ioIndividal, ioContext);
	//Assign the fitness
	ioIndividal.setFitness(lFitness);
	ioIndividal.getFitness()->setValid();
	return;
}


void Beagle::MPI::EvaluationOp::distributeDemeEvaluation(Deme& ioDeme, Context& ioContext) {
	try{
		std::vector<int> lProcess(mProcessSize, -1);
		lProcess[0] = -2; //Master should not be pick
		int lCurrentIndividual = 0;
		std::ostringstream lStreamOut;

		PACC::XML::Streamer lXMLStream(lStreamOut);

		//char lSizeMessage[256];
		int lMessageSize;
		MPI_Status lStatus;
		
		int lFlag;
		unsigned int lSource = 1;
		unsigned int lProcessIdx = 0;
		unsigned int lRecvIndividualIdx = 0;

		unsigned int lNbReceived = 0;		
		unsigned int lNbSent = 0;
		bool lAllSent = false;
		
		while( (lNbReceived < lNbSent) || !lAllSent ) {
			if(!lAllSent) {
				lProcessIdx = find(lProcess, -1, 0, lProcess.size());
				if( lProcessIdx != lProcess.size() ) {
					if((ioDeme[lCurrentIndividual]->getFitness() == NULL) ||
					   (ioDeme[lCurrentIndividual]->getFitness()->isValid() == false)) {
						
						//There is a process idle
						Beagle_LogVerboseM(   
										   ioContext.getSystem().getLogger(),
										   "evaluation", "Beagle::MPIEvaluationOp",
										   std::string("Evaluating the fitness of the ")+uint2ordinal(lCurrentIndividual+1)+
										   " individual"
										   );
						
						ioContext.setIndividualIndex(lCurrentIndividual);
						ioContext.setIndividualHandle(ioDeme[lCurrentIndividual]);
						
						//Send the individual to be evaluated
						Beagle_LogTraceM(
										   ioContext.getSystem().getLogger(),
										   "evaluation", "Beagle::MPIEvaluationOp",
										   std::string("Sending the ") + uint2ordinal(lCurrentIndividual+1) + std::string(" individual to ")+
										   uint2ordinal(lProcessIdx) + std::string(" evaluator")
										   );
						
						lStreamOut.str("");
						ioDeme[lCurrentIndividual]->write(lXMLStream);
						lMessageSize = lStreamOut.str().size()+1;
						MPI_Send(&lMessageSize, 1, MPI_INT, lProcessIdx, eMessageSize, MPI_COMM_WORLD);
						MPI_Send(const_cast<char*>(lStreamOut.str().data()), lMessageSize, MPI_CHAR, lProcessIdx, eIndividual, MPI_COMM_WORLD);
						//std::cout << "Sending individual : " << lStreamOut.str().data() << std::endl;
						unsigned int lGeneration = ioContext.getGeneration();
						MPI_Send(&lGeneration, 1, MPI_INT, lProcessIdx, eIndividual, MPI_COMM_WORLD);
						lProcess[lProcessIdx] = lCurrentIndividual;	
						++lNbSent;
					}
					++lCurrentIndividual;
					if(lCurrentIndividual >= ioDeme.size()) {
						lAllSent = true;
					}
				}
			}
			
			//Look if any cruncher sent a fitness back
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lFlag, &lStatus);
			if (lFlag) {
				//Receive the evaluated fitness
				lSource = lStatus.MPI_SOURCE;
				MPI_Recv(&lMessageSize, 1, MPI_INT, lSource, eMessageSize, MPI_COMM_WORLD, &lStatus);
				char *lMessage = new char[lMessageSize];
				MPI_Recv(lMessage, lMessageSize, MPI_CHAR, lSource, eFitness, MPI_COMM_WORLD, &lStatus);
				lRecvIndividualIdx = lProcess[lSource];
				lProcess[lSource] = -1;
				++lNbReceived;
				
				Beagle_LogTraceM(
								   ioContext.getSystem().getLogger(),
								   "evaluation", "Beagle::MPIEvaluationOp",
								   std::string("Receiving the fitness of the ") + uint2ordinal(lRecvIndividualIdx+1) + 
								   std::string(" individual from ")+uint2ordinal(lSource) + std::string(" evaluator")
								   );
				
				//Read the received fitness
				std::istringstream lStreamIn(lMessage);
				PACC::XML::Document lXMLParser;

				lXMLParser.parse(lStreamIn);
				
				Fitness::Handle lFitness = castHandleT<Fitness>(ioDeme[lRecvIndividualIdx]->getFitnessAlloc()->allocate());
				PACC::XML::ConstIterator lFitnessRootNode = lXMLParser.getFirstRoot(); 

				lFitness->read(lFitnessRootNode);
				
				//Free message space
				delete [] lMessage;
				
				//Assign the fitness
				ioDeme[lRecvIndividualIdx]->setFitness(lFitness);
				ioDeme[lRecvIndividualIdx]->getFitness()->setValid();
				
				//Update stats
				ioContext.setProcessedDeme(ioContext.getProcessedDeme()+1);
				ioContext.setTotalProcessedDeme(ioContext.getTotalProcessedDeme()+1);
				ioContext.setProcessedVivarium(ioContext.getProcessedVivarium()+1);
				ioContext.setTotalProcessedVivarium(ioContext.getTotalProcessedVivarium()+1);  
				
				Beagle_LogDebugM(
								 ioContext.getSystem().getLogger(),
								 "evaluation", "Beagle::MPIEvaluationOp",
								 std::string("Received fitness of individual: ")+
								 ioDeme[lRecvIndividualIdx]->serialize()
								 );
				
				Beagle_LogDebugM(
								   ioContext.getSystem().getLogger(),
								   "evaluation", "Beagle::MPIEvaluationOp",
								   std::string("The individual\'s fitness is: ")+
								   ioDeme[lRecvIndividualIdx]->getFitness()->serialize()
								   );
			}
		}
	} catch(Exception& inException) {
		std::cerr << "Exception catched in evolver:" << std::endl << std::flush;
		std::cerr << inException.what() << std::endl << std::flush;
		exit(1);
	}
	catch(std::exception& inException) {
		std::cerr << "Standard exception catched in evolver:" << std::endl << std::flush;
		std::cerr << inException.what() << std::endl << std::flush;
		exit(1);
	}
}

void Beagle::MPI::EvaluationOp::evolverOperate(Deme& ioDeme, Context& ioContext) {
	Beagle_LogTraceM(
					 ioContext.getSystem().getLogger(),
					 "evaluation", "Beagle::MPIEvaluationOp",
					 std::string("Evaluating the individuals fitness of the ")+
					 uint2ordinal(ioContext.getDemeIndex()+1)+" deme"
					 );
	
	if(ioDeme.size() == 0)
		return;
	
	Individual::Handle lOldIndividualHandle = ioContext.getIndividualHandle();
	unsigned int lOldIndividualIndex = ioContext.getIndividualIndex();    
	
	ioContext.setProcessedDeme(0);
	if((ioContext.getGeneration()!=0) && (ioDeme.getStats()->existItem("total-processed"))) {
		ioContext.setTotalProcessedDeme((unsigned int)ioDeme.getStats()->getItem("total-processed"));
	}
	else ioContext.setTotalProcessedDeme(0);
	ioDeme.getStats()->setInvalid();
	
	if(ioContext.getDemeIndex()==0) {
		Stats& lVivaStats = *ioContext.getVivarium().getStats();
		ioContext.setProcessedVivarium(0);
		if((ioContext.getGeneration()!=0) && (lVivaStats.existItem("total-processed"))) {
			ioContext.setTotalProcessedVivarium((unsigned int)lVivaStats.getItem("total-processed"));
		}
		else ioContext.setTotalProcessedVivarium(0);
		lVivaStats.setInvalid();
	}
	
	
	distributeDemeEvaluation(ioDeme, ioContext);
	
	ioContext.setIndividualIndex(lOldIndividualIndex);
	ioContext.setIndividualHandle(lOldIndividualHandle);
	
	if(mDemeHOFSize->getWrappedValue() > 0) {
		Beagle_LogDetailedM(
							ioContext.getSystem().getLogger(),
							"evaluation", "Beagle::MPIEvaluationOp",
							"Updating the deme's hall-of-fame"
							);
		ioDeme.getHallOfFame().updateWithDeme(mDemeHOFSize->getWrappedValue(), ioDeme, ioContext);
		ioDeme.getHallOfFame().log(Logger::eVerbose, ioContext);
	}
	
	if(mVivaHOFSize->getWrappedValue() > 0) {
		Beagle_LogDetailedM(
							ioContext.getSystem().getLogger(),
							"evaluation", "Beagle::MPIEvaluationOp",
							"Updating the vivarium's hall-of-fame"
							);
		ioContext.getVivarium().getHallOfFame().updateWithDeme(mVivaHOFSize->getWrappedValue(),
															   ioDeme, ioContext);
		ioContext.getVivarium().getHallOfFame().log(Logger::eVerbose, ioContext);
	}
}

void Beagle::MPI::EvaluationOp::evaluatorOperate(Deme& ioDeme, Context& ioContext) {
	try {
		//char lMessage[4096];
		int lMessageSize;
		MPI_Status lStatus;
		int lSource;

		bool lDone = false;
		while(!lDone) {
			//Receive an individual to evaluate
			MPI_Recv(&lMessageSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
			lSource = lStatus.MPI_SOURCE;
			if(lStatus.MPI_TAG == eEvolutionEnd) {
				Beagle_LogDetailedM(
								   ioContext.getSystem().getLogger(),
								   "evaluation", "Beagle::MPIEvaluationOp",
								   std::string("End of evolution received from process ") + int2str(lSource)
								   );
				lDone = true;
			} else {
				char *lMessage = new char[lMessageSize];
				unsigned int lGeneration;
				MPI_Recv(lMessage, lMessageSize, MPI_CHAR, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
				MPI_Recv(&lGeneration, 1, MPI_INT, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
				ioContext.setGeneration(lGeneration);
				Beagle_LogTraceM(
								   ioContext.getSystem().getLogger(),
								   "evaluation", "Beagle::MPIEvaluationOp",
								   std::string("Evaluating individual send from process ") + int2str(lSource)
								   );
				
				//Parse the received individual
				std::istringstream lStreamIn(lMessage);

				PACC::XML::Document lXMLParser;
				lXMLParser.parse(lStreamIn);

				PACC::XML::ConstIterator lIndividualRootNode = lXMLParser.getFirstRoot(); 

				
				//Read the received individual
				ioContext.getDeme().resize(0);
				Individual::Handle lIndividual = castHandleT<Individual>(ioContext.getDeme().getTypeAlloc()->allocate());
				lIndividual->readWithContext(lIndividualRootNode,ioContext);
				ioContext.setIndividualHandle(lIndividual);
				ioContext.setIndividualIndex(0);
				
//				Beagle_LogDebugM(
//								 ioContext.getSystem().getLogger(),
//								 "evaluation", "Beagle::MPIEvaluationOp",
//								 std::string("Individual received: ") + lIndividual->serialize()
//								 );
				
				//Free message string
				delete [] lMessage;
				
				//Evaluated the fitness of the received individual
				Fitness::Handle lFitness = evaluate(*lIndividual, ioContext);
			
				//Send back the fitness
				std::ostringstream lStreamOut;
				PACC::XML::Streamer lXMLStream(lStreamOut);

				lFitness->write(lXMLStream);
				//std::cout << "Sending fitness of size " << lStreamOut.str().size()+1 << ":" << std::endl << lStreamOut.str() << std::endl;
				lMessageSize = lStreamOut.str().size()+1;
				
				Beagle_LogTraceM(
									ioContext.getSystem().getLogger(),
									"evaluation", "Beagle::MPIEvaluationOp",
									std::string("Sending back fitness")
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

/*!
 *  \brief Test the fitness of a given individual.
 *  \param inIndividual Handle to the individual to test.
 *  \param ioSystem Handle to the system to use to test the individual.
 *  \par Note:
 *    This method is provided as a mean to test some individuals after an evolution.
 */
Fitness::Handle Beagle::MPI::EvaluationOp::test(Individual::Handle inIndividual, System::Handle ioSystem)
{
	Beagle_LogInfoM(
					ioSystem->getLogger(),
					"evaluation", "Beagle::MPIEvaluationOp",
					std::string("Testing the following individual: ")+inIndividual->serialize()
					);
    
	Context::Alloc::Handle lContextAlloc =
    castHandleT<Context::Alloc>(ioSystem->getContextAllocatorHandle());
	Context::Handle lContext = castHandleT<Context>(lContextAlloc->allocate());
	lContext->setSystemHandle(ioSystem);
	lContext->setIndividualHandle(inIndividual);
	Fitness::Handle lFitness = evaluate(*inIndividual, *lContext);
	
	Beagle_LogInfoM(
					ioSystem->getLogger(),
					"evaluation", "Beagle::MPIEvaluationOp",
					std::string("New fitness of the individual: ")+lFitness->serialize()
					);
    
	return lFitness;
}



