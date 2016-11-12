/* 
 * File:   energyfunc.hpp
 * Author: medved
 *
 * Created on February 29, 2016, 5:23 PM
 */

#ifndef ENERGYFUNC_HPP
#define ENERGYFUNC_HPP

#include <functor.hpp>
#include <energy.hpp>
namespace CRYSTAL {

    class EnergyFunc : public COMPI::Functor <double> {
    public:

        EnergyFunc(lattice::Energy& energy) : mEnergy(energy) {
        }

        double func(const double* x) {
            return mEnergy.energy(x);
        }

        /**
         * Retrieve energy 
         * @return reference to energy
         */
        lattice::Energy& getEnergy() {
            return mEnergy;
        }

    private:

        lattice::Energy& mEnergy;
    };


}

#endif /* ENERGYFUNC_HPP */

