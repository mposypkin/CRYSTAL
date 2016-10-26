/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits>

//#include <common/fileutils.hpp>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include <methods/gfsdesc/gfsdesc.hpp>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <methods/hookejeeves/coorhjexplorer.hpp>
#include <methods/hookejeeves/rndhjexplorer.hpp>
#include <methods/hookejeeves/hookjeeves.hpp>
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

class VarCoorStopper : public LOCSEARCH::VarCoorDesc<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, const double* gran, double fval, int n) {
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
        else if (n > 1000)
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

class HJStopper : public LOCSEARCH::MHookeJeeves<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double fval, int n) {
        //        std::cout << "n = " << n << "fdiff = " << fdiff << "\n";
        mCnt++;
        return false;
    }

    int mCnt = 0;
};

void getPoint(const snowgoose::Box<double>& box, double* x) {
    int n = box.mDim;
    for (int i = 0; i < n; i++) {
        x[i] = box.mA[i] + (box.mB[i] - box.mA[i]) * (double) rand() / (double) RAND_MAX;
    }
}

void setupPoints(const snowgoose::Box<double>& box, int np, int n, double* x) {
    srand(1);
    for (int i = 0; i < np; i++) {
        double* y = x + i * n;
        getPoint(box, y);
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

    double a[] = {0.8, 0, 0.4, 0.8, 0, 0.4, 0.8, 0, 0.4, 0.8, 0, 0.4};
    double b[] = {1, 0, 4, 1, 4, 4, 1, 4, 4, 1, 4, 4};

    for (int i = 0; i < n; i++) {
        mpp.mBox->mA[i] = a[i];
        mpp.mBox->mB[i] = b[i];
    }
    double x[] = {1.1, 0, 1.1, 1.1, 0.4, 0.7, 1.5, 0, 0.1, 1.2, 0.4, 0.8};
    std::cout << "f(x) = " << enrg->energy(x) << "\n";
}

struct MyLog {

    MyLog() : mItime(time(NULL)) {
        std::cout << "format:\n";
        std::cout << "log-value:time(s):objcalls:value\n";
        std::cout << "log-record:time(s):objcalls:value\n";
    }

    time_t log(const std::string& pstfx, double v, int fcalls) {
        time_t t = time(NULL) - mItime;
        std::cout << "log-" << pstfx << ":";
        std::cout << t << ":";
        std::cout << fcalls << ":";
        std::cout << v << "\n";
        return t;
    }

    time_t mItime;
};

/**
 * 
 */
int main(int argc, char** argv) {
    COMPI::MPProblem<double> mpp;
    setupLJProblem(mpp);
    const int n = mpp.mVarTypes.size();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*(mpp.mObjectives.at(0)));
    mpp.mObjectives.pop_back();
    mpp.mObjectives.push_back(obj);

#if 1       
    DichStopper lstp;
    LOCSEARCH::DichotLS<double> ls(mpp, lstp);
    ls.getOptions().mSInit = .1;
    ls.getOptions().mAccelerate = 1.1;
    ls.getOptions().mSlowDown = 0.5;
#else
    QuadStopper lstp;
    LOCSEARCH::QuadLS<double> ls(mpp, lstp);

#endif    

    GFSStopper gstp;

#if 1      
    LOCSEARCH::GFSDesc<double> gdesc(mpp, gstp, &ls);
#else
    LOCSEARCH::GFSDesc<double> gdesc(mpp, gstp);
#endif    
    gdesc.getOptions().mHInit = .01;

    CoorStopper cstp;
    LOCSEARCH::CoorDesc<double> cdesc(mpp, cstp);

    VarCoorStopper vcstp;
    LOCSEARCH::VarCoorDesc<double> vcdesc(mpp, vcstp);
    vcdesc.getOptions().mHInit = 0.5;
    vcdesc.getOptions().mInc = 3;


    HJStopper hjstp;
#if 1    
    LOCSEARCH::CoorHJExplorer<double> explr(mpp);
    explr.getOptions().mHInit = 0.01;
    explr.getOptions().mHLB = 1e-6;
    explr.getOptions().mResetEveryTime = true;
#endif
#if 0    
    LOCSEARCH::RndHJExplorer<double> explr(mpp);
    explr.getOptions().mMaxTries = 512;
    explr.getOptions().mHInit = 1;
#endif    

#if 1
    LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr);
#else    
    LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr, &ls);
#endif    
    hjdesc.getOptions().mLambda = 1;




    double vbest = std::numeric_limits<double>::max();
    double x[n], xbest[n];
    const int niters = 3;
    MyLog log;

    double avev = 0;
    double avefc = 0;
    double avet = 0;

    double xx[niters * n];
    setupPoints(*(mpp.mBox), niters, n, xx);

    time_t tp = 0;

    for (int i = 0; i < niters; i++) {
        double v;

        snowgoose::VecUtils::vecCopy(n, xx + n * i, x);
        obj->reset();
        int its = 0;
#if 0       
        gstp.mCnt = 0;
        bool rv = gdesc.search(x, v);
        its = gstp.mCnt;
#endif        

#if 1       
        cstp.mCnt = 0;
        bool rv = cdesc.search(x, v);
        its = cstp.mCnt;
#endif        

#if 0
        vcstp.mCnt = 0;
        bool rv = vcdesc.search(x, v);
        its = vcstp.mCnt;
#endif        


#if 0        
        hjstp.mCnt = 0;
        bool rv = hjdesc.search(x, v);
        its = hjstp.mCnt;
#endif        
        auto getave = [](int in, double ave, double nval) {
            double rval = ave * (((double) in) / ((double) in + 1)) + nval / ((double) (in + 1));
            return rval;
        };
        time_t t = log.log("value", v, obj->mCounters.mFuncCalls);
        avev = getave(i, avev, v);
        avefc = getave(i, avefc, (double) obj->mCounters.mFuncCalls);
        avet = getave(i, avet, (double) t - (double) tp);
        tp = t;

        //std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        std::cout << "In " << its << " iterations found v = " << v << "\n";
        std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
        if (v < vbest) {
            vbest = v;
            snowgoose::VecUtils::vecCopy(n, x, xbest);
            log.log("record", v, obj->mCounters.mFuncCalls);
        }
    }

    
    std::cout << "Record = " << vbest << "\n";
    std::cout << "Average value = " << avev << "\n";
    std::cout << "Average objective calls = " << avefc << "\n";
    std::cout << "Average time = " << avet << "\n";

    std::cout << "x = " << snowgoose::VecUtils::vecPrint(n, xbest, 10);
    std::cout << "Energy = " << obj->func(xbest) << "\n";
    
    // TMP
    
    const double tmpx[n] = {0.8576121447,0.0000000000,0.9902911868,0.8576121494,3.4660753358,0.9902787022,0.8576121467,2.9708735663,0.9902911845,0.8576121502,0.4952391918,0.9902787065};
    std::cout << "TMP Energy = " << obj->func(tmpx) << "\n";
    
    std::cout << "Method:\n";

#if 0    
    std::cout << gdesc.about() << "\n";
#endif    

#if 1    
    std::cout << cdesc.about() << "\n";
#endif    

#if 0    
    std::cout << vcdesc.about() << "\n";
#endif    


#if 0    
    std::cout << hjdesc.about() << "\n";
#endif    

    return 0;
}

