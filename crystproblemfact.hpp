/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   crystproblemfact.hpp
 * Author: posypkin
 *
 * Created on January 5, 2017, 9:48 PM
 */

#ifndef CRYSTPROBLEMFACT_HPP
#define CRYSTPROBLEMFACT_HPP

#include <string>
#include <common/sgerrcheck.hpp>
#include <latticedata.hpp>
#include <atoms.hpp>
#include "potentialnames.hpp"
#include "ljproblemfactory.hpp"
#include "tsofproblemfactory.hpp"

/**
 * Crystal problem factory
 */
class CrystallProblemFactory {
public:

    /**
     * Constructor 
     * @param nlayers number of layers
     * @param double pieceLength the horizontal length of the modeling piece
     * @param double radius the cutof radius
     */
    CrystallProblemFactory(std::string potname = LENNARD_JONES_POTENTIAL, int nlayers = 3, double pieceLength = 16, double radius = 3, int atom =  lattice::CARBON)
    : mPotName(potname) {
        setupLatticeData(nlayers, pieceLength, radius, atom);
    }

    COMPI::MPProblem<double>* get() {
        if (mPotName == LENNARD_JONES_POTENTIAL) {
            LJProblemFactory ljpf(mData);
            return ljpf.getLJProblem();
        } else if (mPotName == TERSOFF_POTENTIAL) {
            TsofProblemFactory tsofpf(mData);
            return tsofpf.getTsofProblem();            
        } else {
            SG_ERROR_REPORT("wrong potential specified\n");
        }
    }

private:

    void setupLatticeData(int nlayers, double pieceLength, double radius, int atom) {
        mData.mNumLayers = nlayers;
        mData.mLength = pieceLength;
        mData.mRadius = radius;
        mData.mLayersAtoms.assign(mData.mNumLayers, atom);        
    }

    lattice::LatticeData mData;
    std::string mPotName;
};

#endif /* CRYSTPROBLEMFACT_HPP */

