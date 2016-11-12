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
     * Choose and setup proper potential
     * @param pname potential name
     * @param mpp mathematical programming problem
     */
    static void choosePotential(char* pname, COMPI::MPProblem<double>& mpp) {
        const std::string potent(pname);
        if (potent == LENNARD_JONES_POTENTIAL)
            PotentialSetup::setupLJProblem(mpp);
        else if (potent == TERSOFF_POTENTIAL)
            PotentialSetup::setupTersoffProblem(mpp);
        else
            SG_ERROR_REPORT("wrong potential specified\n");
    }

    /**
     * Setup Lennard-Jones problem
     * @param mpp mathematical programming problem
     */
    static void setupLJProblem(COMPI::MPProblem<double>& mpp) {

        lattice::LatticeData data;

        /**
         * Setup lattice data
         */
        data.mNumLayers = 4;
        data.mLength = 16;
        data.mRadius = 3;
        snowgoose::VecUtils::vecSet(data.mNumLayers, lattice::CARBON, data.mLayersAtoms);

        /**
         * Setup potential
         */
        const double qcut = data.mRadius * data.mRadius;
        const double d = 0.15;
        const double qmin = (data.mRadius - d) * (data.mRadius - d);
        const double qdelta = qcut - qmin;
        // Lennard Jones
        lattice::PotentialCutter pc(qcut, qdelta, lattice::ljpotent);

        /**
         * Setup energy
         */
        lattice::LatticeUtils lut(data);
        lattice::PairPotentialEnergy* enrg = new lattice::PairPotentialEnergy(lut, pc);

        CRYSTAL::EnergyFunc* fe = new CRYSTAL::EnergyFunc(*enrg);
        mpp.mObjectives.push_back(fe);

        const int n = data.mNumLayers * 3;
        for (int i = 0; i < n; i++) {
            int v = COMPI::MPProblem<double>::VariableTypes::GENERIC;
            mpp.mVarTypes.push_back(v);
        }
        mpp.mBox = new snowgoose::Box<double>(n);
        double A = 0;
        double B = 4;

        double a[] = {0.5, 0, 0.4, 0.5, 0, 0.4, 0.5, 0, 0.4, 0.5, 0, 0.4};
        double b[] = {1, 0, 4, 1, 4, 4, 1, 4, 4, 1, 4, 4};

        for (int i = 0; i < n; i++) {
            mpp.mBox->mA[i] = a[i];
            mpp.mBox->mB[i] = b[i];
        }
        double x[] = {1.1, 0, 1.1, 1.1, 0.4, 0.7, 1.5, 0, 0.1, 1.2, 0.4, 0.8};
        std::cout << "f(x) = " << enrg->energy(x) << "\n";
    }

    /**
     * Setup Tersoff potential problem
     * @param mpp mathematical programming problem
     */
    static void setupTersoffProblem(COMPI::MPProblem<double>& mpp) {

        lattice::LatticeData data;

        /**
         * Setup lattice data
         */
        data.mNumLayers = 4;
        data.mLength = 16;
        data.mRadius = 3;
        snowgoose::VecUtils::vecSet(data.mNumLayers, lattice::CARBON, data.mLayersAtoms);


        /**
         * Setup potential
         */
        lattice::TersoffParams tparam;
        lattice::fillCarbonParametersTersoffOriginal(tparam);
        lattice::TersoffUtils tutils(tparam);
        lattice::LatticeUtils lut(data);
        lattice::TersoffEnergy *enrg = new lattice::TersoffEnergy(lut, tutils);


        /**
         * Setup objective
         */
        CRYSTAL::EnergyFunc* fe = new CRYSTAL::EnergyFunc(*enrg);
        mpp.mObjectives.push_back(fe);

        int n = data.mNumLayers * 3;
        for (int i = 0; i < n; i++) {
            int v = COMPI::MPProblem<double>::VariableTypes::GENERIC;
            mpp.mVarTypes.push_back(v);
        }
        mpp.mBox = new snowgoose::Box<double>(n);
        double A = 0;
        double B = 4;

        double a[] = {0.5, 0, 0.4, 0.5, 0, 0.4, 0.5, 0, 0.4, 0.5, 0, 0.4};
        double b[] = {2, 0, 4, 2, 4, 4, 2, 4, 4, 2, 4, 4};

        for (int i = 0; i < n; i++) {
            mpp.mBox->mA[i] = a[i];
            mpp.mBox->mB[i] = b[i];
        }
    }

};

#endif /* POTENTIALSETUP_HPP */

