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
 *  \file   beagle/GA/EvolverBitString.hpp
 *  \brief  Definition of the class GA::EvolverBitString.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.4 $
 *  $Date: 2010/01/19 13:28:36 $
 */

#ifndef Beagle_MPI_GA_EvolverBitString_hpp
#define Beagle_MPI_GA_EvolverBitString_hpp

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
 *  \class EvolverBitString beagle/GA/EvolverBitString.hpp "beagle/GA/EvolverBitString.hpp"
 *  \brief Bit string GA evolver class.
 *  \ingroup GAF
 *  \ingroup GABS
 */
class EvolverBitString : public Beagle::MPI::Evolver {
	
public:
	
	//! GA::EvolverBitString allocator type.
	typedef AllocatorT<EvolverBitString,Beagle::MPI::Evolver::Alloc>
	Alloc;
	//! GA::EvolverBitString handle type.
	typedef PointerT<EvolverBitString,Beagle::MPI::Evolver::Handle>
	Handle;
	//! GA::EvolverBitString bag type.
	typedef ContainerT<EvolverBitString,Beagle::MPI::Evolver::Bag>
	Bag;
	
	explicit EvolverBitString(unsigned int inInitSize=1);
	explicit EvolverBitString(MPI::EvaluationOp::Handle inEvalOp, unsigned int inInitSize=1);         

	explicit EvolverBitString(UIntArray inInitSize);
	explicit EvolverBitString(MPI::EvaluationOp::Handle inEvalOp, UIntArray inInitSize);

	virtual ~EvolverBitString() { }
	

};

}
}
}

#endif // Beagle_GA_EvolverBitString_hpp
