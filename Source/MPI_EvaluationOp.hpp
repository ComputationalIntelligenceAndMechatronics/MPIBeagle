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
 *  \file   beagle/MPIEvaluationOp.hpp
 *  \brief  Definition of the class MPIEvaluationOp.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.4 $
 *  $Date: 2008/10/15 19:46:30 $
 */

#ifndef Beagle_MPI_EvaluationOp_hpp
#define Beagle_MPI_EvaluationOp_hpp

#include <string>
#include <iostream>

#include "beagle/EvaluationOp.hpp"
#include "beagle/config.hpp"
#include "beagle/macros.hpp"
#include "beagle/Object.hpp"
#include "beagle/AbstractAllocT.hpp"
#include "beagle/PointerT.hpp"
#include "beagle/ContainerT.hpp"
#include "beagle/Operator.hpp"
#include "beagle/UInt.hpp"
#include "beagle/System.hpp"
#include "beagle/Context.hpp"
#include "beagle/Logger.hpp"
#include "beagle/BreederOp.hpp"

namespace Beagle {
namespace MPI {
/*!
 *  \class MPIEvaluationOp beagle/MPIEvaluationOp.hpp "beagle/MPIEvaluationOp.hpp"
 *  \brief Abstract evaluation operator class.
 *  \ingroup ECF
 *  \ingroup Op
 *  \ingroup FitStats
 */
class EvaluationOp : public Beagle::EvaluationOp {
	
public:
	
	//! MPIEvaluationOp allocator type.
	typedef AbstractAllocT<EvaluationOp,Beagle::EvaluationOp::Alloc>
	Alloc;
	//! MPIEvaluationOp handle type.
	typedef PointerT<EvaluationOp,Beagle::EvaluationOp::Handle>
	Handle;
	//! MPIEvaluationOp bag type.
	typedef ContainerT<EvaluationOp,Beagle::EvaluationOp::Bag>
	Bag;
	
	explicit EvaluationOp(std::string inName="MPI-EvaluationOp");
	virtual ~EvaluationOp() { }
	
	/*!
	 *  \brief Evaluate the fitness of the given individual.
	 *  \param inIndividual Current individual to evaluate.
	 *  \param ioContext Evolutionary context.
	 *  \return Handle to the fitness value of the individual.
	 */
	virtual Fitness::Handle evaluate(Individual& inIndividual, Context& ioContext) = 0;
	
	virtual Individual::Handle breed(Individual::Bag& inBreedingPool,
									 BreederNode::Handle inChild,
									 Context& ioContext);
	virtual float              getBreedingProba(BreederNode::Handle inChild);
	virtual void               initialize(System& ioSystem);
	virtual void               operate(Deme& ioDeme, Context& ioContext);
	virtual Fitness::Handle    test(Individual::Handle inIndividual, System::Handle ioSystem);
	
protected:
	void evolverOperate(Deme& ioDeme, Context& ioContext);
	void evaluatorOperate(Deme& ioDeme, Context& ioContext);
	void distributeDemeEvaluation(Deme& ioDeme, Context& ioContext);
	void individualEvaluation(Individual& ioIndividal, Context& ioContext);
	
	UInt::Handle mVivaHOFSize;
	UInt::Handle mDemeHOFSize;
	
	int mRank;         //!< MPI rank for this process
	int mProcessSize;  //!< Number of process running 
	
};
}
}

#endif // Beagle_EvaluationOp_hpp
