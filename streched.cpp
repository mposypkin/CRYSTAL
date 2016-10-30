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
#include "energyfunc.hpp"

/**
 * For debugging: 
 * 0.877452 0 0.940712 0.877452 1.32048 0.948432 0.872253 3.65006 0.949996 0.872253 0.372046 0.948432
 * 0.8576121447 0.0000000000 0.9902911868 0.8576121494 3.4660753358 0.9902787022 0.8576121467 2.9708735663 0.9902911845 0.8576121502 0.4952391918 0.9902787065
 */


/**
 * Array size
 */
static const int N = 12;

/**
 * Setup optimization problem for Lennard-Jones clusters
 * @param mpp optimization problem to fill in
 */
void setupLJProblem(COMPI::MPProblem<double>& mpp) {

    lattice::LatticeData * data = new lattice::LatticeData();

    /**
     * Setup lattice data
     */
    data->mNumLayers = 4;
    data->mLength = 16;
    data->mRadius = 3;
    data->mLayersAtoms[0] = lattice::CARBON;
    data->mLayersAtoms[1] = lattice::CARBON;
    data->mLayersAtoms[2] = lattice::CARBON;
    data->mLayersAtoms[3] = lattice::CARBON;


    /**
     * Setup potential
     */
    double qcut = data->mRadius * data->mRadius;
    double d = 0.15;
    double qmin = (data->mRadius - d) * (data->mRadius - d);
    double qdelta = qcut - qmin;
    // Lennard Jones
    lattice::PotentialCutter* pc = new lattice::PotentialCutter(qcut, qdelta, lattice::ljpotent);

    /**
     * Setup energy
     */
    lattice::LatticeUtils* lut = new lattice::LatticeUtils(*data);
    lattice::PairPotentialEnergy* enrg = new lattice::PairPotentialEnergy(*lut, *pc);
    enrg->setFixedAtoms(true);
    CRYSTAL::EnergyFunc* fe = new CRYSTAL::EnergyFunc(*enrg);

    mpp.mObjectives.push_back(fe);

    int n = data->mNumLayers * 3;
    for (int i = 0; i < n; i++) {
        int v = COMPI::MPProblem<double>::VariableTypes::GENERIC;
        mpp.mVarTypes.push_back(v);
    }
    mpp.mBox = new snowgoose::Box<double>(n);
    double A = 0;
    double B = 4;

    double a[] = {0.4, 0, 0.4, 0.4, 0, 0.4, 0.4, 0, 0.4, 0.4, 0, 0.4};
    double b[] = {1, 0, 4, 1, 4, 4, 1, 4, 4, 1, 4, 4};

    for (int i = 0; i < n; i++) {
        mpp.mBox->mA[i] = a[i];
        mpp.mBox->mB[i] = b[i];
    }
    double x[] = {1.1, 0, 1.1, 1.1, 0.4, 0.7, 1.5, 0, 0.1, 1.2, 0.4, 0.8};
    //std::cout << "f(x) = " << enrg->energy(x) << "\n";
}

void stretch(int nh, double* x) {
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
    COMPI::MPProblem<double> mpp;
    setupLJProblem(mpp);
    const int n = mpp.mVarTypes.size();
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
    std::cout << "Energy = " << obj->func(x) << "\n";
    const snowgoose::Box<double>& box = *(mpp.mBox);
    int cnt = 0;
    auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
        cnt++;
        std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
        if (cnt > 2000)
            return true;
        else
            return false;
    };
    auto projector = [&box] (double* x) {
        stretch(3, x);
        snowgoose::BoxUtils::project(x, box);
    };
    LOCSEARCH::ProjCoorDesc<double> desc(mpp, stopper, projector);
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";

    return 0;
}

