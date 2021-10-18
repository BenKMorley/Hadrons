/*
 * SaveOneptTxt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MScalarSUN_SaveOneptTxt_hpp_
#define Hadrons_MScalarSUN_SaveOneptTxt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Save a list of operators to .txt files                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class SaveOneptTxtPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveOneptTxtPar,
                                    std::vector<std::string>, Ops);
};

template <typename SImpl>
class TSaveOneptTxt: public Module<SaveOneptTxtPar>
{
public:
    typedef typename SImpl::Field         Field;
    typedef typename SImpl::ComplexField  ComplexField;

public:
    // constructor
    TSaveOneptTxt(const std::string name);

    // destructor
    virtual ~TSaveOneptTxt(void) {};

    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);

    // execution
    virtual void execute(void);
};


MODULE_REGISTER_TMP(SaveOneptTxtSU2, TSaveOneptTxt<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(SaveOneptTxtSU3, TSaveOneptTxt<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(SaveOneptTxtSU4, TSaveOneptTxt<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(SaveOneptTxtSU5, TSaveOneptTxt<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(SaveOneptTxtSU6, TSaveOneptTxt<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                       TSaveOneptTxt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TSaveOneptTxt<SImpl>::TSaveOneptTxt(const std::string name)
: Module<SaveOneptTxtPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TSaveOneptTxt<SImpl>::getInput(void)
{   
    std::vector<std::string> in;

    for (auto &o: par().Ops)
    {
        in.push_back(o);
    }

    return in;
}

template <typename SImpl>
std::vector<std::string> TSaveOneptTxt<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TSaveOneptTxt<SImpl>::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TSaveOneptTxt<SImpl>::execute(void)
{
    const unsigned int                           nd      = env().getNd();
    const unsigned int                           traj    = vm().getTrajectory();
    std::string                                  op_string;
    FFT                                          fft(envGetGrid(Field));
    std::complex<double>                         peek_buf   (0.0, 0.0);
    std::vector<int>                             site       (nd, 0);

    for (auto &op_string: par().Ops)
    {
        LOG(Message) << "============= Begin Saving ===================" << std::endl;

        auto &op = envGet(ComplexField, op_string);

        LOG(Message) << "Saving Operator:" << op_string << std::endl;

        // Set up file for saving
        std::string string = "/mnt/drive2/Fourier-Laplace/data/g0.1/su2/L32/m2-0.031/config/";
        string += op_string;
        string += "_";
        string += std::to_string(traj);
        std::string string_real = string + "_real";
        std::string string_imag = string + "_imag";

        string_real += ".txt";
        string_imag += ".txt";

        std::ofstream myfile_real;
        std::ofstream myfile_imag;

        myfile_real.open(string_real);
        myfile_imag.open(string_imag);

        LOG(Message) << "About to save to " << string_real << std::endl;
        LOG(Message) << "About to save to " << string_imag << std::endl;

        // Iterate over all lattice sites
        for (int i = 0; i < env().getDim()[0]; i++){
            for (int j = 0; j < env().getDim()[1]; j++){
                for (int k = 0; k < env().getDim()[2]; k++){
                    site = {i, j, k};
                    peekSite(peek_buf, op, site);

                    myfile_real << std::setprecision(std::numeric_limits<double>::digits10 + 1) << peek_buf.real() << '\n';
                    myfile_imag << std::setprecision(std::numeric_limits<double>::digits10 + 1) << peek_buf.imag() << '\n';
                }
            }
        }

        myfile_real.close();
        myfile_imag.close();

        LOG(Message) << "============= End Saving ===================" << std::endl;
    }

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_SaveOneptTxt_hpp_
