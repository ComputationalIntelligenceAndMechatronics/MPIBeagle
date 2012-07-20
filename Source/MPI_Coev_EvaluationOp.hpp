/*
 *  EvaluationOp.h
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

#ifndef MPI_Coev_EvaluationOp_H
#define MPI_Coev_EvaluationOp_H

#include <string>

#include "beagle/config.hpp"
#include "beagle/macros.hpp"
#include "beagle/Object.hpp"
#include "beagle/AbstractAllocT.hpp"
#include "beagle/PointerT.hpp"
#include "beagle/ContainerT.hpp"
#include "beagle/EvaluationOp.hpp"
#include "beagle/GP/Individual.hpp"
#include "beagle/GP/Context.hpp"
#include "beagle/GP/Datum.hpp"

#include "MPI_EvaluationOp.hpp"
#include "beagle/Coev/EvaluationOp.hpp"

namespace Beagle {
namespace MPI {
namespace Coev {

class EvaluationOp : public Beagle::MPI::EvaluationOp {
	
public:
	
	//! Coev::EvaluationOp allocator type.
	typedef AbstractAllocT<EvaluationOp,Beagle::MPI::EvaluationOp::Alloc>
	Alloc;
	//! Coev::EvaluationOp handle type.
	typedef PointerT<EvaluationOp,Beagle::MPI::EvaluationOp::Handle>
	Handle;
	//! Coev::EvaluationOp bag type.
	typedef ContainerT<EvaluationOp,Beagle::MPI::EvaluationOp::Bag>
	Bag;
	
	explicit EvaluationOp(unsigned int inTrigger=1, std::string inName="MPI-Coev-EvaluationOp");
	virtual ~EvaluationOp() { }
	
	//! Vector of evaluation set.
	typedef std::vector< Beagle::Coev::EvaluationOp::EvalSet,BEAGLE_STLALLOCATOR<Beagle::Coev::EvaluationOp::EvalSet> > EvalSetVector;

	
	/*!
	 *  \brief Evaluate fitness of a bunch of individual in an evaluation sets.
	 *  \param ioSets Sets to evaluate fitness.
	 *
	 *  You cannot make assumption about the order of the evaluation sets are in
	 *  the ioSets structure. These can change from evaluation to evaluation. Use
	 *  struct EvalSet member inID to identify the evaluation sets you put into
	 *  shared storage by call to add set. This ID is for convenience and is not used
	 *  internally to identify sets.
	 *
	 *  Assign fitness value to individuals using method assignFitness, otherwise
	 *  statistics value will be erroneous.
	 */
	virtual void evaluateSets(EvalSetVector& ioSets) =0;
	
	/*!
	 *  \brief Make evaluation sets from given deme.
	 *  \param ioIndivBag Bag of individuals to evaluate.
	 *  \param ioContext Evolutionary context.
	 *
	 *  This method is call one time at each generation, for each evolving
	 *  thread/population. This method is problem-specific and consists
	 *  to make evaluation sets that are mated with other thread/population
	 *  evaluation sets for co-evolutionary fitness evaluation. Sets are added
	 *  into a shared structure with a call to method addSet.
	 */
	virtual void makeSets(Beagle::Individual::Bag& ioIndivBag, Beagle::Context::Handle ioContext) =0;
	
	virtual void addSet(Beagle::Coev::EvaluationOp::EvalSet& inEvalSet, bool inBlocking=true);
	
	
	virtual void assignFitness(Beagle::Fitness::Handle inFitness,
							   Beagle::Individual& ioIndividual,
							   Beagle::Context& ioContext) const;
	
	virtual Individual::Handle breed(Individual::Bag& inBreedingPool,
									 BreederNode::Handle inChild,
									 Context& ioContext);
	
	virtual void               initialize(System& ioSystem);
	virtual void               operate(Deme& ioDeme, Context& ioContext);
	virtual Fitness::Handle    test(Individual::Handle inIndividual, System::Handle ioSystem);
	
protected:
	void evolverOperate(Deme& ioDeme, Context& ioContext);
	void distributeIndividuals(std::vector<Individual::Bag>& inIndividuals, Context& ioContext, std::vector<int>& inAssignmentVector);
	
	UInt::Handle mVivaHOFSize;
	UInt::Handle mDemeHOFSize;
	
	int mRank;         //!< MPI rank for this process
	int mProcessSize;  //!< Number of process running 
	
	static PACC::Threading::Condition smCondition;      //!< Condition of co-evaluation
	static EvalSetVector              smEvalSets;       //!< Shared storage of evaluation sets
	static unsigned int               smTrigger;        //!< Number of sets needed to start an evaluation
	
private:
	
	virtual Fitness::Handle evaluate(Individual& inIndividual, Context& ioContext);
	
};
	
}
}
}

#endif
