/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tsofproblemfactory.hpp
 * Author: posypkin
 *
 * Created on January 8, 2017, 12:32 PM
 */

#ifndef TSOFPROBLEMFACTORY_HPP
#define TSOFPROBLEMFACTORY_HPP

#include <ppenergy.hpp>
#include <atoms.hpp>
#include <tsofenergy.hpp>
#include <carbontersoff.hpp>
#include <common/vec.hpp>
#include <mpproblem.hpp>
#include "energyfunc.hpp"

/**
 * Problem factory for Tersoff potential
 */
class TsofProblemFactory {
public:

    TsofProblemFactory(const lattice::LatticeData& data) : mData(data) {
    }

    COMPI::MPProblem<double>* getTsofProblem() {
        COMPI::MPProblem<double>* mpp = new COMPI::MPProblem<double>();
        
        /**
         * Setup potential
         */
        lattice::TersoffParams tparam;
        //lattice::fillCarbonParametersTersoffOriginal(tparam);
        //lattice::fillCarbonParametersTersoffOptimized(tparam);
        lattice::fillCarbonParametersPowell(tparam);
        lattice::TersoffUtils tutils(tparam);
        lattice::LatticeUtils lut(mData);
        lattice::TersoffEnergy *enrg = new lattice::TersoffEnergy(lut, tutils);


        /**
         * Setup objective
         */
        CRYSTAL::EnergyFunc* fe = new CRYSTAL::EnergyFunc(*enrg);
        mpp->mObjectives.push_back(fe);

        int n = mData.mNumLayers * 3;
        mpp->mVarTypes.assign(n, COMPI::MPProblem<double>::VariableTypes::GENERIC);
        mpp->mBox = new snowgoose::Box<double>(n);

        for (int i = 0; i < mData.mNumLayers; i++) {
            const int j = 3 * i;
            // height
            mpp->mBox->mA[j] = 0.5;
            mpp->mBox->mB[j] = 2;
            // displacement
            if (i == 0) {
                mpp->mBox->mA[j + 1] = 0;
                mpp->mBox->mB[j + 1] = 0;
            } else {
                mpp->mBox->mA[j + 1] = 0;
                mpp->mBox->mB[j + 1] = 4;

            }
            // stride
            mpp->mBox->mA[j + 2] = 0.4;
            mpp->mBox->mB[j + 2] = 4;
        }
        return mpp;
    }

private:

    const lattice::LatticeData& mData;
};

#endif /* TSOFPROBLEMFACTORY_HPP */

