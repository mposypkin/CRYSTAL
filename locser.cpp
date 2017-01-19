/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   locser.cpp
 * Author: posypkin
 *
 * Created on January 5, 2017, 1:52 PM
 */

#include <methods/coordesc/coordesc.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include "crystproblemfact.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    CrystallProblemFactory cpf(argv[1]);
    COMPI::MPProblem<double>& mpp = *cpf.get();
    int cnt = 0;
    auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
        cnt++;
        //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
        if (cnt > 7000)
            return true;
        else
            return false;
    };
    LOCSEARCH::CoorDesc<double> desc(mpp, stopper);
    const int n = 12;
    double *x = new double[n];
    snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
    double v;
    v = mpp.mObjectives[0]->func(x);
    std::cout << "Initial v = " << v << "\n";
    std::cout << "Initial x = " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    return 0;
}

