#include <iostream>
#include <cstdlib>

int main() {
    freopen("output.txt", "w", stdout);
    system("g++ -O3 -fopenmp -march=native ep-wocd.cpp -o ep-wocd");
    for(int i = 1; i <= 30; i++) {
        system("./ep-wocd");
    }
    
    return 0;
}