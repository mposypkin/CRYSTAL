/*
 * Computes the energy of the perturbed structure
 */

/* 
 * File:   streched.cpp
 * Author: posypkin
 *
 * Created on October 23, 2016, 10:53 PM
 */

#include <iostream>
#include <sstream>
#include <common/vec.hpp>

/**
 * For debugging: 
 * 0.877452 0 0.940712 0.877452 1.32048 0.948432 0.872253 3.65006 0.949996 0.872253 0.372046 0.948432
 */


/**
 * Array size
 */
static const int N = 12;

/*
 * 
 */
int main(int argc, char** argv) {
    
    std::cout << "Enter string\n";
    std::string s;
    std::getline(std::cin, s);
    double x[N];
    snowgoose::VecUtils::vecRead(s, N, x);
    std::cout << snowgoose::VecUtils::vecPrint(N, x) << "\n";
    return 0;
}

