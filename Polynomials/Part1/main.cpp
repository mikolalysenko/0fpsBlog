// Compile with:
//
//    g++ --std=c++0x main.cpp
//

#include <iostream>
#include "symmetric_tensor.h"

using namespace std;

int main(int argc, char** argv) {

	SymmetricIndexIterator<3,4> iter, end(-1);


	//Test basic iteration
	cout << "Array\tTensor\tDegree\tComputed Array" << endl;
	while(iter != end) {
	
		auto tsr = iter.tensor_index();
		auto deg = iter.degree_index();
		auto arr = iter.array_index();
	
		//Print out array index
		cout << arr << "\t";
	
		//Print out tensor index
		cout << "(";
		for(int i=0; i<tsr.size(); ++i) {
			if(i > 0) cout << ',';
			cout << tsr[i];
		}
		cout << ")\t";
		
		//Print out degree index
		cout << "(";
		for(int i=0; i<deg.size(); ++i) {
			if(i > 0) cout << ',';
			cout << deg[i];
		}
		cout << ")\t";
		
		//Print out computed index
		cout << iter.degree_to_array(deg) << endl;
		
		++iter;
	}
	
	return 0;
}
