#include <iostream>
#include "symmetric_tensor.h"

using namespace std;

int main(int argc, char** argv) {

    //Create the linear polynomial, 3x + 5y"
    SymmetricTensor<float, 1, 2> linear_poly;
    linear_poly({1,0}) = 3;
    linear_poly({0,1}) = 5;
    cout << "Linear poly = " << linear_poly << endl;

    SymmetricTensor<float, 1, 2> B;
    B({1,0}) = 0;
    B({0,1}) = 1;
    cout << "B = " << B << endl;
    cout << "A*B = " << linear_poly.outer_product(B) << endl;

    //Create a quadratic polynomial, 2x^2 + xy - 3y^2
    SymmetricTensor<float, 2, 2> quadratic_poly;
    quadratic_poly({2,0}) = 2;
    quadratic_poly({1,1}) = 1;
    quadratic_poly({0,2}) = -3;
    cout << "Quadratic poly = " << quadratic_poly << endl;
    
    
    //Compute product
    auto prod = quadratic_poly.outer_product(linear_poly);
    cout << "Product = " << prod << endl;
    
    return 0;
}

