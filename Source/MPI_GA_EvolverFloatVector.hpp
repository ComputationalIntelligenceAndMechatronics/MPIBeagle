/*
 *  MPI_GA_EvolverFloatVector.h
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

#ifndef MPI_GA_EvolverFloatVector_H
#define MPI_GA_EvolverFloatVector_H

class MPI_GA_EvolverFloatVector {
public:
	MPI_GA_EvolverFloatVector() {}
	~MPI_GA_EvolverFloatVector() {}
};

#include <string>

#include "beagle/config.hpp"
#include "beagle/macros.hpp"
#include "beagle/Object.hpp"
#include "beagle/Pointer.hpp"
#include "beagle/PointerT.hpp"
#include "beagle/Allocator.hpp"
#include "beagle/AllocatorT.hpp"
#include "beagle/ContainerT.hpp"
#include "beagle/Evolver.hpp"
#include "beagle/EvaluationOp.hpp"
#include "beagle/TerminationOp.hpp"


#include "MPI_Evolver.hpp"
#include "MPI_EvaluationOp.hpp"

namespace Beagle {
namespace MPI {
namespace GA {

/*!
 *  \class EvolverFloatVector beagle/GA/EvolverFloatVector.hpp "beagle/GA/EvolverFloatVector.hpp"
 *  \brief Bit string GA evolver class.
 *  \ingroup GAF
 *  \ingroup GABS
 */
class EvolverFloatVector : public Beagle::MPI::Evolver {
	
public:
	
	//! GA::EvolverFloatVector allocator type.
	typedef AllocatorT<EvolverFloatVector,Beagle::MPI::Evolver::Alloc>
	Alloc;
	//! GA::EvolverFloatVector handle type.
	typedef PointerT<EvolverFloatVector,Beagle::MPI::Evolver::Handle>
	Handle;
	//! GA::EvolverFloatVector bag type.
	typedef ContainerT<EvolverFloatVector,Beagle::MPI::Evolver::Bag>
	Bag;
	
	explicit EvolverFloatVector(unsigned int inInitSize=1);
	explicit EvolverFloatVector(MPI::EvaluationOp::Handle inEvalOp, unsigned int inInitSize=1);
	explicit EvolverFloatVector(UIntArray inInitSize);
	explicit EvolverFloatVector(MPI::EvaluationOp::Handle inEvalOp, UIntArray inInitSize);
	virtual ~EvolverFloatVector() { }
	
	virtual void initialize(System::Handle ioSystem, int ioArgc, char** ioArgv);
	virtual void initialize(System::Handle ioSystem, std::string inConfigFilename);
};

}
}
}

#endif
