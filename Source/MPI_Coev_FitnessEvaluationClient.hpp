/*
 *  MPI_Coev_FitnessEvaluationClient.h
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

#ifndef MPI_Coev_FitnessEvaluationClient_H
#define MPI_Coev_FitnessEvaluationClient_H

#include <beagle/Allocator.hpp>
#include <beagle/Fitness.hpp>
#include <beagle/IndividualBag.hpp>
#include <beagle/System.hpp>
#include <vector>

namespace Beagle {
namespace MPI {
namespace Coev {

class FitnessEvaluationClient {
public:
	FitnessEvaluationClient() {}
	~FitnessEvaluationClient() {}
	
	void operate(Beagle::Allocator::Handle inContextAllocator, 
				 std::vector<Beagle::Genotype::Alloc::Handle>& inGenotypeAlloc, 
				 Beagle::Fitness::Alloc::Handle inFitnessAlloc);
	
	virtual void initialize(int& ioArgc, char** ioArgv);
	virtual void initialize(std::string inConfigFilename);
	
protected:
	virtual void init() = 0;
	virtual void postInit() = 0;
	virtual Beagle::Fitness::Handle evaluate(Beagle::Individual::Bag& inIndividuals, Beagle::Context& ioContext) = 0;
	Beagle::System::Handle mSystem;
};
	
}
}
}
#endif
