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
#include <methods/projcoordesc/projcoordesc.hpp>
#include <ppenergy.hpp>
#include <atoms.hpp>
#include <pairpotentials.hpp>
#include <funccnt.hpp>
#include "crystproblemfact.hpp"

/**
 * For debugging: 
 * 0.877452 0 0.940712 0.877452 1.32048 0.948432 0.872253 3.65006 0.949996 0.872253 0.372046 0.948432
 * Lennard-Jones: 0.8576121447 0.0000000000 0.9902911868 0.8576121494 3.4660753358 0.9902787022 0.8576121467 2.9708735663 0.9902911845 0.8576121502 0.4952391918 0.9902787065
 * Tersoff 
 * [0.7301206930,0.0000000000,2.5285255149,1.4604847987,0.0000758707,2.5285088817,0.7301213132,1.2642696344,2.5285238656,1.4604815266,1.2641880933,2.5285509477]
 * 0.7301206930 0.0000000000 2.5285255149 1.4604847987 0.0000758707 2.5285088817 0.7301213132 1.2642696344 2.5285238656 1.4604815266 1.2641880933 2.5285509477
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

int main(int argc, char** argv) {

    if (argc != 3) {
        std::cout << "usage: " << argv[0] << LENNARD_JONES_POTENTIAL << "|" << TERSOFF_POTENTIAL << " stretch_factor \n";
        exit(-1);
    }

    CrystallProblemFactory cpf(argv[1]);
    COMPI::MPProblem<double>& mpp = *cpf.get();

    const double stretchFactor = atof(argv[2]);

    std::cout.precision(10);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);


    const int n = mpp.mVarTypes.size();

    CRYSTAL::EnergyFunc* ef = dynamic_cast<CRYSTAL::EnergyFunc*> (mpp.mObjectives.at(0));
    SG_ASSERT(ef);
    ef->getEnergy().setFixedAtoms(true);

    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*(mpp.mObjectives.at(0)));
    mpp.mObjectives.pop_back();
    mpp.mObjectives.push_back(obj);

    std::cout << "Enter  string\n";
    std::string s;
    std::getline(std::cin, s);
    double x[N];
    snowgoose::VecUtils::vecRead(s, N, x);
    std::cout << snowgoose::VecUtils::vecPrint(N, x) << "\n";
    std::cout << "Energy = " << obj->func(x) << "\n";


    auto computeHeight = [] (double* x) {
        double h = 0;
        for (int i = 0; i < N; i += 3) {
            h += x[i];
        }
        return h;
    };



    const snowgoose::Box<double>& box = *(mpp.mBox);
    int cnt = 0;
    auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
        cnt++;
        //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
        if (cnt > 2000)
            return true;
        else
            return false;
    };


    const double newh = stretchFactor * computeHeight(x);
    std::cout << "newh = " << newh << "\n";

    auto projectorOld = [&newh, &box] (double* x) {
        stretch(newh, x);
        snowgoose::BoxUtils::project(x, box);
    };


    auto projector = [&box, &projectorOld, &computeHeight] (double* x) {
        projectorOld(x);
        //        std::cout << "height = " << computeHeight(x) << "\n";
        //        std::cout << " at " << snowgoose::VecUtils::vecPrint(N, x) << "\n";
    };

    LOCSEARCH::ProjCoorDesc<double> desc(mpp, stopper, projector);
    double v;
    bool rv = desc.search(x, v);

    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    std::cout << " with height " << computeHeight(x) << "\n";

    return 0;
}

