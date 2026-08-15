#pragma once

#include <vector>
#include <volt/math/lin_alg.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/neighbor_bond_array.h>

namespace Volt{

template<typename iterator>
void bitmapSort(iterator begin, iterator end, int max){

	int bitarray = 0;
	for(iterator pin = begin; pin != end; ++pin){
		bitarray |= 1 << (*pin);
	}

	iterator pout = begin;
	for(int i = max - 1; i >= 0; i--){
		if(bitarray & (1 << i)){
			*pout++ = i;
		}
	}

}

struct CoordinationStructure{
    int numNeighbors;
    std::vector<Vector3> latticeVectors;
    NeighborBondArray neighborArray;
    int cnaSignatures[MAX_NEIGHBORS];
    int commonNeighbors[MAX_NEIGHBORS][2];
};

}
