/*
 * Hello_World.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MScalarSUN_Hello_World_hpp_
#define Hadrons_MScalarSUN_Hello_World_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Trace of powers of a scalar field                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class Hello_WorldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Hello_WorldPar,
                                    std::string, name);
};

template <typename SImpl>
class THello_World: public Module<Hello_WorldPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;

public:
    // constructor
    THello_World(const std::string name);

    // destructor
    virtual ~THello_World(void) {};

    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);

    // execution
    virtual void execute(void);
};

// Hello_WorldSU2 is instantiated by 59 with SImpl = ScalarNxN...., 
MODULE_REGISTER_TMP(Hello_WorldSU2, THello_World<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(Hello_WorldSU3, THello_World<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(Hello_WorldSU4, THello_World<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(Hello_WorldSU5, THello_World<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(Hello_WorldSU6, THello_World<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                          THello_World implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
THello_World<SImpl>::THello_World(const std::string name)
: Module<Hello_WorldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> THello_World<SImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> THello_World<SImpl>::getOutput(void)
{
    std::vector<std::string> out;

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void THello_World<SImpl>::setup(void)
{}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void THello_World<SImpl>::execute(void)
{
    LOG(Message) << "Hello World! " << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Hello_World_hpp_
