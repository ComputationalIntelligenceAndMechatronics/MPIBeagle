/*
 *  VectorUtil.h
 *  Copyright 2008 Jean-Francois Dupuis.
 *  
 *  This is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with it.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  This file was created by Jean-Francois Dupuis on 08/07/2008.
 */

#pragma once

#include <vector>

using namespace std;

template <class T>

/*! \brief Find value in vector based on subscripts
 *  \param  inVector Vector in which the value are searched
 *  \param  inStart Start subscript of the search
 *  \param  inEnd End subscript of the search
 *  \return The subscript of the first value found in vector. Return inVector.size() if not found.
 */
unsigned int find(vector<T> inVector, T inValue, unsigned int inStart, unsigned int inEnd) {
	for(unsigned int i = inStart; i < inEnd; ++i) {
		if(inVector[i] == inValue)
			return i;
	}
	//Value not found
	return inVector.size();
}
