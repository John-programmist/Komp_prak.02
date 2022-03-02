#include "Matrix.h"
#include <cmath>

int main(){
    Vector vector1(VERTICAL, 2);
    
    SquareMatrix a(3);
    
    std::cin >> a;
    
    std::cout << a << std::endl;
    
    std::cout << a.det() << std::endl;
    
    
    
    
    
    
    
    return EXIT_SUCCESS;
}
