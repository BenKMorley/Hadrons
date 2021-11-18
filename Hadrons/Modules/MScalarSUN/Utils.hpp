/*
 * Utils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MScalarSUN_Utils_hpp_
#define Hadrons_MScalarSUN_Utils_hpp_

#define PI           3.14159265358979323846


#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MScalarSUN)

GRID_SERIALIZABLE_ENUM(DiffType, undef, forward, 1, backward, 2, central, 3);

template <typename Field>
inline void dmu(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRONS_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out = in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out = Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out = 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRONS_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

template <typename Field>
inline void dmuAcc(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRONS_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out += in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out += Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out += 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRONS_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

inline double WindowBeta(double a, double b, double r)
{
    if (r==a || r==b)
    {
        return 0.0;
    } else 
    {
        return exp(-1.0/(r-a)-1.0/(b-r));
    }
}


inline double trapezoidalIntegral(double a, double b, double r, int n, const std::function<double (double, double, double)> &f) {
    const double width = (r-a)/n;

    double trapezoidal_integral = 0;
    for(int step = 0; step < n; step++) {
        const double x1 = a + step*width;
        const double x2 = a + (step+1)*width;

        trapezoidal_integral += 0.5*(x2-x1)*(f(a, b, x1) + f(a, b, x2));
    }

    return trapezoidal_integral;
}


inline double windowFunction(double a, double b, double r, int n_step)
{
    if (r <= a)
    {
        return 0.0;
    }
    else if (r >= b)
    {
        return 1.0;
    } else
    {
        return trapezoidalIntegral(a, b, r, n_step, &WindowBeta)/trapezoidalIntegral(a, b, b, n_step, &WindowBeta);
    }
}


inline void MakeWindowField(LatticeComplex &field, double windowmin, double windowmax, int n_step)
{
    auto &env = Environment::getInstance();
    GridBase *grid = field.Grid();
    auto latt_size = grid->GlobalDimensions();

    std::vector<int> shift(env.getNd(), 0);
    int ri, rj, rk;
    double r;
    Complex buf; 
    for (int i = 0; i < latt_size[0]; i++)
    {
        for (int j = 0; j < latt_size[1]; j++)
        {
            for (int k = 0; k < latt_size[2]; k++)
            {
                ri = std::min(i, latt_size[0] - i);
                rj = std::min(j, latt_size[1] - j);
                rk = std::min(k, latt_size[2] - k);
                r = pow(ri*ri+rj*rj+rk*rk, 0.5);
                shift = {i, j, k};
                buf = windowFunction(windowmin, windowmax, r, n_step);
                pokeSite(buf, field, shift);
            }
        }
    }
}


inline double LaplaceTransform3D(LatticeComplex &field, LatticeComplex &temp, int offset, const int L, bool debug){
    std::vector<std::complex<double>>                        p_s(L, 0);
    std::vector<std::complex<double>>                        x_s(L, 0);
    int                                                      i;
    int                                                      j;
    int                                                      p;
    int                                                      x;
    int                                                      y;
    int                                                      z;
    int                                                      r;
    int                                                      d;
    auto&                                                    env = Environment::getInstance();
    int                                                      nd = env.getNd();
    std::vector<int>                                         qt(nd,0);
    std::vector<std::vector<std::complex<double>>>           exp_comb(L);
    std::vector<int>                                         fetch(nd,0);
    std::vector<int>                                         set(nd,0);
    std::complex<double>                                     sum;
    std::complex<double>                                     fetch_buf;

    for (p = 0; p < L; p ++){
        p_s[p] = (PI * 2 * p / L);
    }

    for (x = 0; x < L - offset; x++){
        x_s[x] = x;
    }

    for (x = L - offset; x < L; x++){
        x_s[x] = x - L;
    }

    // Set up the exponant pre-factors to the Laplace Transform
    for (i = 0; i < L; i++) {
        std::vector<std::complex<double>> e(L, 0);
        exp_comb[i] = e;
    }

    thread_for(r, L * L, {
        i = r / L;
        j = r % L;
        exp_comb[i][j] = exp(-p_s[i] * x_s[j]);
    })

    for (d = 0; d < 3; d++){
        /* Vectorize x and y directions of the for loop */
        // for(r = 0; r < L * L; r ++){
        thread_for(r, L * L, {
            // Each thread needs its own copy of these variables
            std::vector<int>                                         fetch_thread(nd,0);
            std::vector<int>                                         set_thread(nd,0);
            std::complex<double>                                     fetch_buf_thread;
            int                                                      i_thread = 0;
            int                                                      j_thread = 0;
            int                                                      x_thread = 0;
            int                                                      y_thread = 0;
            int                                                      z_thread = 0;
            std::complex<double>                                     sum_thread = 0;

            y_thread = r / L;
            z_thread = r % L;

            for (x_thread = 0; x_thread < L; x_thread++){
                /* Sum will count along the contracted index of the sum */
                sum_thread = 0;

                for (i_thread = 0; i_thread < L; i_thread++){
                    fetch_thread = {i_thread, x_thread, y_thread};
                    peekSite(fetch_buf_thread, field, fetch_thread);

                    // Note that x, y and z here are just being used as indices
                    sum_thread += exp_comb[z_thread][i_thread] * fetch_buf_thread;
                }

                /* Copy to the temp field */
                set_thread = {x_thread, y_thread, z_thread};
                pokeSite(sum_thread, temp, set_thread);
            }
        })
    
        if (debug == true){
            LOG(Message) << "DEBUG 1" << std::endl;

            for (int j = 0; j < L; j++){
                set = {0, 0, j};
                peekSite(fetch_buf, temp, set);
                LOG(Message) << fetch_buf << std::endl;
            }
        }
        
        /* Only after the loop has finished do we copy temp into field */
        thread_for(r, L * L, {
            // Each thread needs its own copy of these variables
            std::vector<int>                                         fetch_thread(nd,0);
            std::complex<double>                                     fetch_buf_thread;
            int                                                      x_thread = 0;
            int                                                      y_thread = 0;
            int                                                      z_thread = 0;

            x_thread = r / L;
            y_thread = r % L;

            for (z_thread = 0; z_thread < L; z_thread++){
                fetch_thread = {x_thread, y_thread, z_thread};
                peekSite(fetch_buf_thread, temp, fetch_thread);
                pokeSite(fetch_buf_thread, field, fetch_thread);
            }
        })

        if (debug == true){
            LOG(Message) << "DEBUG 2" << std::endl;

            for (int j = 0; j < L; j++){
                set = {0, 0, j};
                peekSite(fetch_buf, field, set);
                LOG(Message) << fetch_buf << std::endl;
            }
        }
    }

    return 0;
}


template <class SinkSite, class SourceSite>
std::vector<Complex> makeTwoPoint(const std::vector<SinkSite>   &sink,
                                  const std::vector<SourceSite> &source,
                                  const double factor = 1.)
{
    assert(sink.size() == source.size());
    
    unsigned int         nt = sink.size();
    std::vector<Complex> res(nt, 0.);
    
    for (unsigned int dt = 0; dt < nt; ++dt)
    {
        for (unsigned int t  = 0; t < nt; ++t)
        {
            res[dt] += trace(sink[(t+dt)%nt]*adj(source[t]));
        }
        res[dt] *= factor/static_cast<double>(nt);
    }
    
    return res;
}

inline std::string varName(const std::string name, const std::string suf)
{
    return name + "_" + suf;
}

inline std::string varName(const std::string name, const unsigned int mu)
{
    return varName(name, std::to_string(mu));
}

inline std::string varName(const std::string name, const unsigned int mu, 
                           const unsigned int nu)
{
    return varName(name, std::to_string(mu) + "_" + std::to_string(nu));
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Utils_hpp_
