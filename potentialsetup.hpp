/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pairpotentialsetup.hpp
 * Author: posypkin
 *
 * Created on November 12, 2016, 12:58 PM
 */

#ifndef POTENTIALSETUP_HPP
#define POTENTIALSETUP_HPP

#include <ppenergy.hpp>
#include <atoms.hpp>
#include <pairpotentials.hpp>
#include <funccnt.hpp>
#include <common/vec.hpp>
#include <tsofenergy.hpp>
#include <carbontersoff.hpp>
#include "energyfunc.hpp"

// Potential names

#define LENNARD_JONES_POTENTIAL "lj"

#define TERSOFF_POTENTIAL "tsof"

/**
 * Setup pair potentials
 */
class PotentialSetup {
public:

    /**
     * Setup potential factory
     * @param data lattice data
     */
    PotentialSetup(const lattice::LatticeData& data) : mData(data) {



        /**
         * Setup lattice data
         *
        mData.mNumLayers = 4;
        mData.mLength = 16;
        mData.mRadius = 3;
        snowgoose::VecUtils::vecSet(mData.mNumLayers, lattice::CARBON, mData.mLayersAtoms);
         * */
    }

    /**
     * Choose and setup proper potential
     * @param pname potential name
     * @param data lattice data
     * @param mpp mathematical programming problem
     */
    static void choosePotential(char* pname, const lattice::LatticeData& data, COMPI::MPProblem<double>& mpp) {
        const std::string potent(pname);
        PotentialSetup ps(data);
        if (potent == LENNARD_JONES_POTENTIAL)
            ps.setupLJProblem(mpp);
        else if (potent == TERSOFF_POTENTIAL)
            ps.setupTersoffProblem(mpp);
        else
            SG_ERROR_REPORT("wrong potential specified\n");
    }

    /**
     * Setup Lennard-Jones problem
     * @param mpp mathematical programming problem
     */
    void setupLJProblem(COMPI::MPProblem<double>& mpp) {


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
        mpp.mObjectives.push_back(fe);

        const int n = mData.mNumLayers * 3;
        for (int i = 0; i < n; i++) {
            int v = COMPI::MPProblem<double>::VariableTypes::GENERIC;
            mpp.mVarTypes.push_back(v);
        }
        mpp.mBox = new snowgoose::Box<double>(n);
        
        for (int i = 0; i < mData.mNumLayers; i++) {
            const int j = 3 * i;
            // height
            mpp.mBox->mA[j] = 0.5;
            mpp.mBox->mB[j] = 1;
            // displacement
            if (i == 0) {
                mpp.mBox->mA[j + 1] = 0;
                mpp.mBox->mB[j + 1] = 0;
            } else {
                mpp.mBox->mA[j + 1] = 0;
                mpp.mBox->mB[j + 1] = 4;

            }
            // stride
            mpp.mBox->mA[j + 2] = 0.4;
            mpp.mBox->mB[j + 2] = 4;
        }
        
    }

    /**
     * Setup Tersoff potential problem
     * @param mpp mathematical programming problem
     */
    void setupTersoffProblem(COMPI::MPProblem<double>& mpp) {



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
        mpp.mObjectives.push_back(fe);

        int n = mData.mNumLayers * 3;
        for (int i = 0; i < n; i++) {
            int v = COMPI::MPProblem<double>::VariableTypes::GENERIC;
            mpp.mVarTypes.push_back(v);
        }
        mpp.mBox = new snowgoose::Box<double>(n);
               
        for (int i = 0; i < mData.mNumLayers; i++) {
            const int j = 3 * i;
            // height
            mpp.mBox->mA[j] = 0.5;
            mpp.mBox->mB[j] = 2;
            // displacement
            if (i == 0) {
                mpp.mBox->mA[j + 1] = 0;
                mpp.mBox->mB[j + 1] = 0;
            } else {
                mpp.mBox->mA[j + 1] = 0;
                mpp.mBox->mB[j + 1] = 4;

            }
            // stride
            mpp.mBox->mA[j + 2] = 0.4;
            mpp.mBox->mB[j + 2] = 4;
        }
    }

private:

    const lattice::LatticeData& mData;

};

#endif /* POTENTIALSETUP_HPP */

