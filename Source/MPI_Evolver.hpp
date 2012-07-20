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
 *  \file   beagle/MPI/Evolver.hpp
 *  \brief  Definition of the class MPIEvolver.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.4 $
 *  $Date: 2010/01/19 13:28:36 $
 */

#ifndef Beagle_MPI_Evolver_hpp
#define Beagle_MPI_Evolver_hpp

#include "beagle/config.hpp"
#include "beagle/macros.hpp"
#include "beagle/Object.hpp"
#include "beagle/Pointer.hpp"
#include "beagle/PointerT.hpp"
#include "beagle/Allocator.hpp"
#include "beagle/AllocatorT.hpp"
#include "beagle/Container.hpp"
#include "beagle/ContainerT.hpp"
#include "beagle/AbstractAllocT.hpp"
#include "beagle/WrapperT.hpp"
#include "beagle/System.hpp"
#include "beagle/Operator.hpp"
#include "beagle/Map.hpp"
#include "beagle/OperatorMap.hpp"
#include "beagle/Vivarium.hpp"
#include "beagle/ConfigurationDumper.hpp"

//#include "beagle/IntegerVector.hpp"
#include "beagle/Int.hpp"

namespace Beagle {
namespace MPI {
/*!
 *  \class MPIEvolver beagle/MPIEvolver.hpp "beagle/MPIEvolver.hpp"
 *  \brief Beagle's basic evolver class.
 *  \ingroup ECF
 *  \ingroup Op
 */
class Evolver : public Beagle::Evolver {
	
public:
	
	//! MPIEvolver allocator type.
	typedef AllocatorT<Evolver,Beagle::Evolver::Alloc>
	Alloc;
	//! MPIEvolver handle type.
	typedef PointerT<Evolver,Beagle::Evolver::Handle>
	Handle;
	//! MPIEvolver bag type.
	typedef ContainerT<Evolver,Beagle::Evolver::Bag>
	Bag;
	
	Evolver();
	Evolver(Beagle::EvaluationOp::Handle inEvalOp);
	virtual ~Evolver();

	virtual void           initialize(System::Handle ioSystem, int& ioArgc, char** ioArgv);
	virtual void           initialize(System::Handle ioSystem, string inConfigFilename);
	virtual void           evolve(Vivarium::Handle ioVivarium);

protected:
	
	void evolver(Vivarium::Handle ioVivarium);
	void evaluater(Vivarium::Handle ioVivarium);
	void stopEvaluater();
	
	int mRank;			       //!< MPI rank for this process
	Int::Handle mProcessSize;  //!< Number of process running 
	Beagle::EvaluationOp::Handle mEvaluator; //!< Evaluation operator used to evaluate individual
};
}	
}

#endif // Beagle_MPI_Evolver_hpp
