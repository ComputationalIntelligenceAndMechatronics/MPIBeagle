/*
 *  MPI_Coev_EvaluationOp.cpp
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
 *  This file was created by Jean-Francois Dupuis on 13/01/10.
 */

#include "MPI_Coev_EvaluationOp.hpp"

#include "beagle/Coev.hpp"
#include "beagle/Beagle.hpp"
#include <beagle/GA.hpp>
#include "MPI_EvaluationOp.hpp"
#include <mpi.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <XML.hpp>

#include "CommunicationMPI.h"
#include "VectorUtil.h"

using namespace Beagle;


/*
 *  Initialize static evaluation operator evaluation condition.
 */
PACC::Threading::Condition Beagle::MPI::Coev::EvaluationOp::smCondition;

/*
 *  Initialize static evaluation operator sets stack.
 */
Coev::EvaluationOp::EvalSetVector Beagle::MPI::Coev::EvaluationOp::smEvalSets;

/*
 *  Initialize static evaluation operator trigger.
 */
unsigned int Beagle::MPI::Coev::EvaluationOp::smTrigger = 0;

/*!
 *  \brief Construct a co-evolutionary evaluation operator.
 *  \param inTrigger Number of sets to accumulate before triggering evaluation procedure.
 *  \param inName Name of the co-evolutionary evaluation operator.
 */
Beagle::MPI::Coev::EvaluationOp::EvaluationOp(unsigned int inTrigger, Beagle::string inName) :
Beagle::MPI::EvaluationOp(inName)
{
	Beagle_StackTraceBeginM();
	smCondition.lock();
	if(smTrigger == 0) smTrigger = inTrigger;
	else if(inTrigger != smTrigger) {
		std::ostringstream lOSS;
		lOSS << "trigger value given as argument to constructor of Coev::EvaluationOp (";
		lOSS << inTrigger << ") is different from the actual non-zero value of the trigger (";
		lOSS << smTrigger << ")!";
		smCondition.unlock();
		throw Beagle_RunTimeExceptionM(lOSS.str().c_str());
	}
	smCondition.unlock();
	Beagle_StackTraceEndM("Coev::EvaluationOp::EvaluationOp(unsigned int inTrigger, std::string inName)");
}


/*!
 *  \brief Add evaluation set into shared structure.
 *  \param inEvalSet Evaluation set to add.
 *  \param inBlocking If true, the add set operation block until fitness evaluation is done.
 */
void Beagle::MPI::Coev::EvaluationOp::addSet(Beagle::Coev::EvaluationOp::EvalSet& inEvalSet, bool inBlocking)
{
	Beagle_StackTraceBeginM();
	smCondition.lock();
	if(smTrigger == 0) {
		smCondition.unlock();
		throw Beagle_RunTimeExceptionM("co-evolution trigger value is zero!");
	}
	if(smEvalSets.size() >= smTrigger) {
		std::ostringstream lOSS;
		lOSS << "number of evaluation sets in co-evolution evaluation operator (";
		lOSS << smEvalSets.size() << ") is equal or bigger than the trigger value (";
		lOSS << smTrigger << ")!";
		smCondition.unlock();
		throw Beagle_RunTimeExceptionM(lOSS.str().c_str());
	}
	// Add evaluation set
	smEvalSets.push_back(inEvalSet);
	// If all needed evaluation sets are added
	if(smEvalSets.size() == smTrigger) {  
		evaluateSets(smEvalSets);
		smEvalSets.clear();
		smCondition.broadcast();
	}
	// Othewize, wait to get all needed evaluation sets
	else if(inBlocking) {
		smCondition.wait();
	}
	smCondition.unlock();
	Beagle_StackTraceEndM("void Coev::EvaluationOp::addSet(EvalSet& inEvalSet, bool inBlocking)");
}


/*!
 *  \brief Set fitness value of an individual.
 *  \param inFitness New fitness value of the individual.
 *  \param ioIndividual Individual which fitness value is set.
 *  \param ioContext Evolutionary context.
 */
void Beagle::MPI::Coev::EvaluationOp::assignFitness(Fitness::Handle inFitness,
													Individual& ioIndividual,
													Context& ioContext) const
{
	Beagle_StackTraceBeginM();
	ioIndividual.setFitness(inFitness);
	inFitness->setValid();
	ioContext.setProcessedDeme(ioContext.getProcessedDeme()+1);
	ioContext.setTotalProcessedDeme(ioContext.getTotalProcessedDeme()+1);
	ioContext.setProcessedVivarium(ioContext.getProcessedVivarium()+1);
	ioContext.setTotalProcessedVivarium(ioContext.getTotalProcessedVivarium()+1);
	Beagle_StackTraceEndM("void Coev::EvaluationOp::assignFitness(Fitness::Handle inFitness, Individual& ioIndividual, Context& ioContext) const");
}


/*!
 *  \brief Apply the evaluation operation on a breeding pool, returning a evaluated bred individual.
 *  \param inBreedingPool Breeding pool to use for the breeding operation.
 *  \param inChild Node handle associated to child node in the breeder tree.
 *  \param ioContext Evolutionary context of the breeding operation.
 *  \return Evaluated bred individual.
 */
Individual::Handle Beagle::MPI::Coev::EvaluationOp::breed(Individual::Bag& inBreedingPool,
														  BreederNode::Handle inChild,
														  Context& ioContext)
{
	Beagle_StackTraceBeginM();
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
						   "evaluation", "Beagle::Coev::EvaluationOp",
						   "Evaluating the fitness of a new bred individual"
						   );
		
		Individual::Bag lNonEvalIndiv;
		lNonEvalIndiv.push_back(lBredIndividual);
		Context::Handle lContext = &ioContext;
		makeSets(lNonEvalIndiv, lContext);
		
		Beagle_LogVerboseM(
						   ioContext.getSystem().getLogger(),
						   "evaluation", "Beagle::Coev::EvaluationOp",
						   std::string("The individual fitness value is: ")+
						   lBredIndividual->getFitness()->serialize()
						   );
		
		if(mDemeHOFSize->getWrappedValue() > 0) {
			Beagle_LogVerboseM(
							   ioContext.getSystem().getLogger(),
							   "evaluation", "Beagle::Coev::EvaluationOp",
							   "Updating the deme hall-of-fame"
							   );
			lDeme.getHallOfFame().updateWithIndividual(mDemeHOFSize->getWrappedValue(),
													   *lBredIndividual, ioContext);
		}
		if(mVivaHOFSize->getWrappedValue() > 0) {
			Beagle_LogVerboseM(
							   ioContext.getSystem().getLogger(),
							   "evaluation", "Beagle::Coev::EvaluationOp",
							   "Updating the vivarium hall-of-fame"
							   );
			ioContext.getVivarium().getHallOfFame().updateWithIndividual(mVivaHOFSize->getWrappedValue(),
																		 *lBredIndividual, ioContext);
		}
	}
	
	return lBredIndividual;
	Beagle_StackTraceEndM("Individual::Handle Coev::EvaluationOp::breed(Individual::Bag& inBreedingPool, BreederNode::Handle inChild, Context& ioContext)");
}


/*!
 *  \brief Evaluating an individual in co-evolution is not that simple. Define makeSets and
 *     evaluateSets methods instead.
 */
Fitness::Handle Beagle::MPI::Coev::EvaluationOp::evaluate(Individual& inIndividual, Context& ioContext)
{
	Beagle_StackTraceBeginM();
	throw Beagle_UndefinedMethodInternalExceptionM("evaluate","Coev::EvaluationOp",getName());
	Beagle_StackTraceEndM("Fitness::Handle Coev::EvaluationOp::evaluate(Individual& inIndividual, Context& ioContext)");
}

/*!
 *  \brief Initialize the parameters evaluation operator.
 *  \param ioSystem System to use to initialize the operator.
 */
void Beagle::MPI::Coev::EvaluationOp::initialize(System& ioSystem)
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
void Beagle::MPI::Coev::EvaluationOp::operate(Deme& ioDeme, Context& ioContext)
{
	//	if(mProcessSize == 1)
	//		Beagle::EvaluationOp::operate(ioDeme,ioContext);
	//	else {
	if(mRank == 0) { 
		evolverOperate(ioDeme, ioContext);
	}
	else {
		evaluatorOperate(ioDeme, ioContext);
	}
	//	}
}

/*!
 *  \brief Apply co-evolutionary evaluation operation on the deme.
 *  \param ioDeme Deme to evaluate fitness.
 *  \param ioContext Evolutionary context.
 */
void Beagle::MPI::Coev::EvaluationOp::evolverOperate(Deme& ioDeme, Context& ioContext)
{
	Beagle_StackTraceBeginM();
	
	Beagle_LogTraceM(
					 ioContext.getSystem().getLogger(),
					 "evaluation", "Beagle::Coev::EvaluationOp",
					 std::string("Evaluating the individuals fitness of the ")+
					 uint2ordinal(ioContext.getDemeIndex()+1)+" deme in co-evolution mode"
					 );
	
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
	
	Beagle_LogVerboseM(
					   ioContext.getSystem().getLogger(),
					   "evaluation", "Beagle::Coev::EvaluationOp",
					   std::string("Calling co-evolution evaluation hook for the ")+
					   uint2ordinal(ioContext.getDemeIndex()+1)+" deme"
					   );
	
	Context::Handle lContext = &ioContext;
	makeSets(ioDeme, lContext);
	
	if(mDemeHOFSize->getWrappedValue() > 0) {
		Beagle_LogDetailedM(
							ioContext.getSystem().getLogger(),
							"evaluation", "Beagle::Coev::EvaluationOp",
							"Updating the deme's hall-of-fame"
							);
		ioDeme.getHallOfFame().updateWithDeme(mDemeHOFSize->getWrappedValue(), ioDeme, ioContext);
		ioDeme.getHallOfFame().log(Logger::eVerbose, ioContext);
	}
	
	if(mVivaHOFSize->getWrappedValue() > 0) {
		Beagle_LogDetailedM(
							ioContext.getSystem().getLogger(),
							"evaluation", "Beagle::Coev::EvaluationOp",
							"Updating the vivarium's hall-of-fame"
							);
		ioContext.getVivarium().getHallOfFame().updateWithDeme(mVivaHOFSize->getWrappedValue(),
															   ioDeme, ioContext);
		ioContext.getVivarium().getHallOfFame().log(Logger::eVerbose, ioContext);
	}
	Beagle_StackTraceEndM("void Coev::EvaluationOp::operate(Deme& ioDeme, Context& ioContext)");
}

//void Beagle::MPI::Coev::EvaluationOp::evaluatorOperate(Deme& ioDeme, Context& ioContext) {
//	try {
//		//char lMessage[4096];
//		int lMessageSize;
//		int lNbIndividuals = 0;
//		MPI_Status lStatus;
//		int lSource;
//		
//		bool lDone = false;
//		while(!lDone) {
//			//Receive an individual to evaluate
//			MPI_Recv(&lNbIndividuals, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
//			lSource = lStatus.MPI_SOURCE;
//			if(lStatus.MPI_TAG == eEvolutionEnd) {
//				Beagle_LogDetailedM(
//									ioContext.getSystem().getLogger(),
//									"evaluation", "Beagle::MPIEvaluationOp",
//									std::string("End of evolution received from process ") + int2str(lSource)
//									);
//				lDone = true;
//			} else {
//				Individual::Bag lIndividuals;
//				for(unsigned int i = 0; i < lNbIndividuals; ++i) {
//					MPI_Recv(&lMessageSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
//					char *lMessage = new char[lMessageSize];
//					MPI_Recv(lMessage, lMessageSize, MPI_CHAR, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
//					
//					//Parse the received individual
//					std::istringstream lStreamIn(lMessage);
//					PACC::XML::Document lXMLParser;
//					lXMLParser.parse(lStreamIn);
//					PACC::XML::ConstIterator lIndividualRootNode = lXMLParser.getFirstRoot(); 
//					
//					//Read the received individual
//					ioContext.getDeme().resize(0);
//					Individual::Handle lIndividual = castHandleT<Individual>(ioContext.getDeme().getTypeAlloc()->allocate());
//					lIndividual->readWithContext(lIndividualRootNode,ioContext);
//				
//					//Free message std::string
//					delete [] lMessage;
//					
//					lIndividuals.push_back(lIndividual);
//				}
//				unsigned int lGeneration;
//				MPI_Recv(&lGeneration, 1, MPI_INT, lSource, MPI_ANY_TAG, MPI_COMM_WORLD, &lStatus);
//				ioContext.setGeneration(lGeneration);
//				
//				Beagle_LogTraceM(
//								 ioContext.getSystem().getLogger(),
//								 "evaluation", "Beagle::MPIEvaluationOp",
//								 std::string("Evaluating individual send from process ") + int2str(lSource)
//								 );
//				
//				
//				//Evaluated the fitness of the received individual
//				Fitness::Handle lFitness =  evaluateIndividuals(lIndividuals,ioContext);
//				//Fitness::Handle lFitness = evaluate(*lIndividual, ioContext);
//				
//				
//				//Send back the fitness
//				std::ostringstream lStreamOut;
//				PACC::XML::Streamer lXMLStream(lStreamOut);
//				
//				lFitness->write(lXMLStream);
//				//std::cout << "Sending fitness of size " << lStreamOut.str().size()+1 << ":" << std::endl << lStreamOut.str() << std::endl;
//				lMessageSize = lStreamOut.str().size()+1;
//				
//				Beagle_LogTraceM(
//								 ioContext.getSystem().getLogger(),
//								 "evaluation", "Beagle::MPIEvaluationOp",
//								 std::string("Sending back fitness")
//								 );
//				
//				MPI_Send(&lMessageSize, 1, MPI_INT, lSource, eMessageSize, MPI_COMM_WORLD);
//				MPI_Send(const_cast<char*>(lStreamOut.str().data()), lMessageSize, MPI_CHAR, lSource, eFitness, MPI_COMM_WORLD);
//			}
//		}
//	} catch(Exception& inException) {
//		std::cerr << "Exception catched in evaluator:" << std::endl << std::flush;
//		std::cerr << inException.what() << std::endl << std::flush;
//		exit(1);
//	}
//	catch(std::exception& inException) {
//		std::cerr << "Standard exception catched in evaluator:" << std::endl << std::flush;
//		std::cerr << inException.what() << std::endl << std::flush;
//		exit(1);
//	}
//}



/*! \brief Distribute the individual pairs to be evaluated
 *  This method receives a vector of each pairs of individual to be
 *	evaluated together. Each pair are received by the evaluator to
 *	compute the fitness.
 *	\param inAssignmentVector	Indicate which individual of the group to assign new fitness. 0 meaning all individual.
 */
void Beagle::MPI::Coev::EvaluationOp::distributeIndividuals(vector<Individual::Bag>& inIndividuals, Context& ioContext, std::vector<int>& inAssignmentVector) {
	try{
		Beagle_LogVerboseM(   
						   ioContext.getSystem().getLogger(),
						   "evaluation", "Beagle::MPI::Coev::EvaluationOp",
						   std::string("Distributing the ")+uint2str(inIndividuals.size())+std::string(" individuals group for evaluation.")
						   );
		
		std::vector<int> lProcess(mProcessSize, -1);
		lProcess[0] = -2; //Master should not be pick
		int lCurrentIndGroup = 0;
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
					
					//There is a process idle
					Beagle_LogVerboseM(   
									   ioContext.getSystem().getLogger(),
									   "evaluation", "Beagle::MPI::Coev::EvaluationOp",
									   std::string("Evaluating the fitness of the ")+uint2ordinal(lCurrentIndGroup+1)+
									   " individual"
									   );
					
					//						ioContext.setIndividualIndex(lCurrentIndGroup);
					//						ioContext.setIndividualHandle(inIndividuals[lCurrentIndGroup]);
					
					//Send the individual to be evaluated
					Beagle_LogTraceM(
									 ioContext.getSystem().getLogger(),
									 "evaluation", "Beagle::MPI::Coev::EvaluationOp",
									 std::string("Sending the ") + uint2ordinal(lCurrentIndGroup+1) + std::string(" individual group to ")+
									 uint2ordinal(lProcessIdx) + std::string(" evaluator")
									 );
					unsigned int lNbIndividual = inIndividuals[lCurrentIndGroup].size();
					MPI_Send(&lNbIndividual, 1, MPI_INT, lProcessIdx, eNbIndividual, MPI_COMM_WORLD);
					for(unsigned int i = 0; i < inIndividuals[lCurrentIndGroup].size(); ++i) {
						inIndividuals[lCurrentIndGroup][i]->getFitness()->setInvalid();
						lStreamOut.str("");
						inIndividuals[lCurrentIndGroup][i]->write(lXMLStream);
						lMessageSize = lStreamOut.str().size()+1;
						MPI_Send(&lMessageSize, 1, MPI_INT, lProcessIdx, eMessageSize, MPI_COMM_WORLD);
						MPI_Send(const_cast<char*>(lStreamOut.str().data()), lMessageSize, MPI_CHAR, lProcessIdx, eIndividual, MPI_COMM_WORLD);
					}
					unsigned int lGeneration = ioContext.getGeneration();
					MPI_Send(&lGeneration, 1, MPI_INT, lProcessIdx, eIndividual, MPI_COMM_WORLD);
					lProcess[lProcessIdx] = lCurrentIndGroup;	
					++lNbSent;
					++lCurrentIndGroup;
	
					if(lCurrentIndGroup >= inIndividuals.size()) {
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
				
				Fitness::Handle lFitness = castHandleT<Fitness>(inIndividuals[lRecvIndividualIdx][0]->getFitnessAlloc()->allocate());
				PACC::XML::ConstIterator lFitnessRootNode = lXMLParser.getFirstRoot(); 
				
				lFitness->read(lFitnessRootNode);
				
				//Free message space
				delete [] lMessage;
				
				//Assign the fitness
				if(inAssignmentVector[lRecvIndividualIdx] == 0) {
					//Assign the computed value to all individual in the group
					for(unsigned int i = 0; i < inIndividuals[lRecvIndividualIdx].size(); ++i) {
						inIndividuals[lRecvIndividualIdx][i]->setFitness(lFitness);
						inIndividuals[lRecvIndividualIdx][i]->getFitness()->setValid();
					}
				} else {
					//Assign the fitness the correct individual
					inIndividuals[lRecvIndividualIdx][inAssignmentVector[lRecvIndividualIdx]-1]->setFitness(lFitness);
					inIndividuals[lRecvIndividualIdx][inAssignmentVector[lRecvIndividualIdx]-1]->getFitness()->setValid();
				}
				
				//Update stats
				ioContext.setProcessedDeme(ioContext.getProcessedDeme()+1);
				ioContext.setTotalProcessedDeme(ioContext.getTotalProcessedDeme()+1);
				ioContext.setProcessedVivarium(ioContext.getProcessedVivarium()+1);
				ioContext.setTotalProcessedVivarium(ioContext.getTotalProcessedVivarium()+1);  
				
				Beagle_LogTraceM(
								 ioContext.getSystem().getLogger(),
								 "evaluation", "Beagle::MPIEvaluationOp",
								 std::string("The individual\'s fitness is: ")+
								 lFitness->serialize()
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

/*!
 *  \brief Test the fitness of a given individual.
 *  \param inIndividual Handle to the individual to test.
 *  \param ioSystem Handle to the system to use to test the individual.
 *  \par Note:
 *    This method is provided as a mean to test some individuals after an evolution.
 */
Fitness::Handle Beagle::MPI::Coev::EvaluationOp::test(Individual::Handle inIndividual, System::Handle ioSystem)
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

