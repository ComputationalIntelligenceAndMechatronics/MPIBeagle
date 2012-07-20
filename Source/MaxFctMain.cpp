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
 *  \file   MaxFctMain.cpp
 *  \brief  Implementation of the main routine for the maximisation problem.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.3 $
 *  $Date: 2009/11/26 14:41:50 $
 */

#include "beagle/GA.hpp"
#include "MaxFctEvalOp.hpp"
#include "MPI_GA_EvolverBitString.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <numeric>

using namespace std;
using namespace Beagle;

/*!
 *  \brief Main routine for the function maximisation problem.
 *  \param argc Number of arguments on the command-line.
 *  \param argv Arguments on the command-line.
 *  \return Return value of the program.
 *  \ingroup MaxFct
 */
int main(int argc, char** argv) {
	try {
		// 1. Encoding values.

		UIntArray lEncodingOfX;

		lEncodingOfX.resize(5,25);
		const unsigned int lNumberOfBits = std::accumulate(lEncodingOfX.begin(), lEncodingOfX.end(), 0);
		// 2. Build the system.
		System::Handle lSystem = new System;
		// 3. Build evaluation operator.
		MaxFctEvalOp::Handle lEvalOp = new MaxFctEvalOp(lEncodingOfX);
		// 4. Instanciate the evolver and the vivarium for bit string GA population.
		GA::BitString::Alloc::Handle lBSAlloc = new GA::BitString::Alloc;
		Vivarium::Handle lVivarium = new Vivarium(lBSAlloc);
		Beagle::MPI::GA::EvolverBitString::Handle lEvolver = new Beagle::MPI::GA::EvolverBitString(lEvalOp, lNumberOfBits);
		// 5. Initialize the evolver and evolve the vivarium.
		lEvolver->initialize(lSystem, argc, argv);
		lEvolver->evolve(lVivarium);
	}
	catch(Exception& inException) {
		//inException.terminate(cerr);
		cerr << "Exception catched:" << endl << flush;
		cerr << inException.what() << endl << flush;
		return 1;
	}
	catch(std::exception& inException) {
		cerr << "Standard exception catched:" << endl << flush;
		cerr << inException.what() << endl << flush;
		return 1;
	}
	return 0;
}
