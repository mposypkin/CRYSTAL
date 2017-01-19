/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ljproblemfactory.hpp
 * Author: posypkin
 *
 * Created on January 5, 2017, 9:23 PM
 */

#ifndef LJPROBLEMFACTORY_HPP
#define LJPROBLEMFACTORY_HPP

#include <ppenergy.hpp>
#include <atoms.hpp>
#include <pairpotentials.hpp>
#include <common/vec.hpp>
#include <mpproblem.hpp>
#include "energyfunc.hpp"

/**
 * Problem factory to get instance of the Lennard-Jones problem
 */
class LJProblemFactory {
public:

    LJProblemFactory(const lattice::LatticeData& data) : mData(data) {
    }

    COMPI::MPProblem<double>* getLJProblem() {

        COMPI::MPProblem<double>* mpp = new COMPI::MPProblem<double>();

        /**
         * Setup potential
         */
        const double qcut = mData.mRadius * mData.mRadius;
        const double d = 0.15;
        const double qmin = (mData.mRadius - d) * (mData.mRadius - d);
        const double qdelta = qcut - qmin;
        // Lennard Jones
        lattice::PotentialCutter pc(qcut, qdelta, lattice::ljpotent);

        /**
         * Setup energy
         */
        lattice::LatticeUtils lut(mData);
        lattice::PairPotentialEnergy* enrg = new lattice::PairPotentialEnergy(lut, pc);

        CRYSTAL::EnergyFunc* fe = new CRYSTAL::EnergyFunc(*enrg);
        mpp->mObjectives.push_back(fe);

        const int n = mData.mNumLayers * 3;
        mpp->mVarTypes.assign(n, COMPI::MPProblem<double>::VariableTypes::GENERIC);
        mpp->mBox = new snowgoose::Box<double>(n);

        for (int i = 0; i < mData.mNumLayers; i++) {
            const int j = 3 * i;
            // height
            mpp->mBox->mA[j] = 0.5;
            mpp->mBox->mB[j] = 1;
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

#endif /* LJPROBLEMFACTORY_HPP */

