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
 *  \file   MaxFctEvalOp.hpp
 *  \brief  Definition of the type MaxFctEvalOp.
 *  \author Christian Gagne
 *  \author Marc Parizeau
 *  $Revision: 1.3 $
 *  $Date: 2009/11/26 14:41:50 $
 */

/*!
 *  \defgroup MaxFct Function Maximization Example
 *  \brief Function maximization (maxfct): A simple GA example with Open BEAGLE.
 *
 *  \par Objective
 *  Find the maximum of the following 5D function:
 *  \f$f(x) = \frac{161.8}{(u_N^2 \sum_{k=0}^{N-1}(x_k^2 + u_k^2))}\f$
 *  with \f$x = <x_0, x_1, ..., x_{N-1}>\f$, \f$u_{k+1} = x_k + u_k\f$, \f$x_k\f$
 *  in \f$[-200,200]\f$ for all \f$k\f$, \f$N = 5\f$ and \f$u_0 = 10\f$.
 *
 *  \par Representation
 *  Bit strings of 125 bits, constructed from the function's five arguments \f$x_i\f$,
 *  each encoded with 25 bits on the interval \f$[-200,200]\f$.
 *
 *  \par Fitness
 *  Value of \f$f(x_1,x_2,x_3,x_4,x_5)\f$.
 *
 */
 
#ifndef MaxFctEvalOp_hpp
#define MaxFctEvalOp_hpp

#include <vector>
#include "beagle/GA.hpp"
#include "MPI_EvaluationOp.hpp"


/*!
 *  \class MaxFctEvalOp MaxFctEvalOp.hpp "MaxFctEvalOp.hpp"
 *  \brief The individual evaluation class operator for the problem of function maximisation.
 *  \ingroup MaxFct
 */
class MaxFctEvalOp : public Beagle::MPI::EvaluationOp {

public:

  //! MaxFctEvalOp allocator type.
	typedef Beagle::AllocatorT<MaxFctEvalOp,Beagle::MPI::EvaluationOp::Alloc>
          Alloc;
  //!< MaxFctEvalOp handle type.
  typedef Beagle::PointerT<MaxFctEvalOp,Beagle::MPI::EvaluationOp::Handle>
          Handle;
  //!< MaxFctEvalOp bag type.
  typedef Beagle::ContainerT<MaxFctEvalOp,Beagle::MPI::EvaluationOp::Bag>
          Bag;


	explicit MaxFctEvalOp(Beagle::UIntArray inEncoding=Beagle::UIntArray(5,25),
						  std::string inName="MaxFctEvalOp");


  virtual Beagle::Fitness::Handle evaluate(Beagle::Individual& inIndividual,
                                           Beagle::Context& ioContext);

protected:
  std::vector<Beagle::GA::BitString::DecodingKey> mDecodingKeys; //!< Decoding keys.

};

#endif // MaxFctEvalOp_hpp
