/*
 * Computes the energy of the set of atoms
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
#include <methods/coordesc/coordesc.hpp>
#include <ppenergy.hpp>
#include <atoms.hpp>
#include <pairpotentials.hpp>
#include <funccnt.hpp>
#include "crystproblemfact.hpp"

/**
 * For debugging: 
 * 0.877452 0 0.940712 0.877452 1.32048 0.948432 0.872253 3.65006 0.949996 0.872253 0.372046 0.948432
 * Lennard-Jones: 0.8576121447 0.0000000000 0.9902911868 0.8576121494 3.4660753358 0.9902787022 0.8576121467 2.9708735663 0.9902911845 0.8576121502 0.4952391918 0.9902787065
 * Tersoff [0.7125554220,0.0000000000,2.2857142816,1.4536654278,0.0000000567,2.2857142775,0.7125554121,1.1217308149,2.2927564072,1.4533956522,1.1217307869,2.2927564059]
 * 0.7125554220 0.0000000000 2.2857142816 1.4536654278 0.0000000567 2.2857142775 0.7125554121 1.1217308149 2.2927564072 1.4533956522 1.1217307869 2.2927564059
 */


/**
 * Array size
 */
static const int N = 12;

void stretch(double nh, double* x) {
    double h = 0;
    for (int i = 0; i < N; i += 3) {
        h += x[i];
    }
    double a = nh / h;
    for (int i = 0; i < N; i += 3) {
        x[i] *= a;
    }
    h = 0;
    for (int i = 0; i < N; i += 3) {
        h += x[i];
    }

}

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc != 2) {
        std::cout << "usage: " << argv[0] << " " << LENNARD_JONES_POTENTIAL << "|" << TERSOFF_POTENTIAL << "\n";
        exit(-1);
    }

    CrystallProblemFactory cpf(argv[1]);
    COMPI::MPProblem<double>& mpp = *cpf.get();

    const int n = mpp.mVarTypes.size();

    CRYSTAL::EnergyFunc* ef = dynamic_cast<CRYSTAL::EnergyFunc*> (mpp.mObjectives.at(0));
    SG_ASSERT(ef);
    ef->getEnergy().setFixedAtoms(true);

    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*(mpp.mObjectives.at(0)));
    mpp.mObjectives.pop_back();
    mpp.mObjectives.push_back(obj);

    std::cout.precision(10);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);

    std::cout << "Enter  string\n";
    std::string s;
    std::getline(std::cin, s);
    double x[N];
    snowgoose::VecUtils::vecRead(s, N, x);
    std::cout << snowgoose::VecUtils::vecPrint(N, x) << "\n";
    std::cout << "Energy = " << obj->func(x) << "\n";
    if (std::string(argv[1]) == TERSOFF_POTENTIAL) {
        lattice::TersoffEnergy& tsen = (lattice::TersoffEnergy&)ef->getEnergy();
        for (auto i : tsen.getLBounds()) {
            std::cout << i << "\n";
        }
        for (auto i : tsen.getUBounds()) {
            std::cout << i << "\n";
        }
    }

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
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    // std::cout << " with height " << computeHeight(x) << "\n";

    return 0;
}

