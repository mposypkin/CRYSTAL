/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <stdlib.h>
#include <iostream>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include <methods/gfsdesc/gfsdesc.hpp>
#include <methods/coordesc/coordesc.hpp>
#include <ppenergy.hpp>
#include <atoms.hpp>
#include <pairpotentials.hpp>
#include "energyfunc.hpp"


class CoorStopper : public LOCSEARCH::CoorDesc<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gran, double fval, int n) {
        mCnt++;
        //std::cout << "n = " << n << " fval = " << fval << "\n";
        return false;
    }

    int mCnt = 0;
};

class GFSStopper : public LOCSEARCH::GFSDesc<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gmin, double fval, int n) {
        mCnt++;
        if (gmin > -1e-5)
            return true;
        else if(n > 1000) 
            return true;
        else
            return false;
    }

    int mCnt = 0;
};

class QuadStopper : public LOCSEARCH::QuadLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn) {
        mCnt++;
#if 0        
        std::cout << "s = " << s << ", k = " << k << "\n";
        std::cout << "vo = " << vo << ", vn = " << vn << "\n";
#endif
        if (s < 1e-3)
            return true;
        else if (k > 16)
            return true;
        else
            return false;
    }

    int mCnt = 0;
};

class DichStopper : public LOCSEARCH::DichotLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn) {
        mCnt++;
#if 0        
        std::cout << "s = " << s << ", k = " << k << "\n";
        std::cout << "vo = " << vo << ", vn = " << vn << "\n";
#endif
        if (s < 1e-3)
            return true;
        else if (k > 512)
            return true;
        else
            return false;
    }

    int mCnt = 0;
};

void getPoint(const snowgoose::Box<double>& box, double* x) {
    int n = box.mDim;
    for (int i = 0; i < n; i++) {
        x[i] = (double) rand() / (double) RAND_MAX;
    }
}

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
    
    double a[] = {0.8, 0, 0.4,  0.8, 0, 0.4,  0.8, 0, 0.4,  0.8, 0, 0.4};
    double b[] = {1, 0, 4,  1, 4, 4,  1, 4, 4,  1, 4, 4};
    
    for (int i = 0; i < n; i++) {
        mpp.mBox->mA[i] = a[i];
        mpp.mBox->mB[i] = b[i];
    }
    double x[] = {1.1, 0, 1.1, 1.1, 0.4, 0.7, 1.5, 0, 0.1, 1.2, 0.4, 0.8};
    std::cout << "f(x) = " << enrg->energy(x) << "\n";
}

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 100;

    COMPI::MPProblem<double> mpp;
    setupLJProblem(mpp);
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*(mpp.mObjectives.at(0)));
    mpp.mObjectives.pop_back();
    mpp.mObjectives.push_back(obj);

#if 0       
    DichStopper lstp;
    LOCSEARCH::DichotLS<double> ls(mpp, lstp);
    ls.getOptions().mSInit = .1;
    ls.getOptions().mAccelerate = 1.1;
    ls.getOptions().mSlowDown = 0.5;
#else
    QuadStopper lstp;
    LOCSEARCH::QuadLS<double> ls(mpp, lstp);

#endif    

    GFSStopper stp;

#if 1      
    LOCSEARCH::GFSDesc<double> desc(mpp, stp, &ls);
#else
    LOCSEARCH::GFSDesc<double> desc(mpp, stp);
#endif    

#if 0           
    desc.getOptions().mOnlyCoordinateDescent = true;
#endif
    
#if 1
    CoorStopper cstp;
    LOCSEARCH::CoorDesc<double> cdesc(mpp, cstp);
#endif    

    desc.getOptions().mHInit = .01;

    double x[] = {1, 0, 1, 1, 0.5, 1, 1, 0, 1, 1, 0.5, 1};

    int niters = 100;
    for (int i = 0; i < niters; i++) {
        double v;
        getPoint(*(mpp.mBox), (double*)x);
        stp.mCnt = 0;
        obj->reset();
        bool rv = cdesc.search(x, v);
        std::cout << "In " << cstp.mCnt << " iterations found v = " << v << "\n";
        //std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    }
    return 0;
}

