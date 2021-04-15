
#include "minimizer.hpp"

#include <cmath>
#include <cstdio>

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>

// Lazy
using arr = std::vector<double>;
using arrref = std::vector<double>&;

void read_data(std::string fname, arrref t1, arrref y1, arrref dy1, arrref t2, arrref y2, arrref dy2);
int sin_f (const gsl_vector *x, void *data, gsl_vector *f);
int sin_df (const gsl_vector *x, void *data, gsl_matrix *J);
int sin_fvv (const gsl_vector *x, const gsl_vector *v, void *data, gsl_vector *fvv);

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: ./executable filename\n";
        return 1;
    }

    arr t1, t2, y1, y2, dy1, dy2;
    std::string input_fname(argv[1]);

    // The input is of the form A.BC_khz.dat
    // So if we parse it, we just take the first 4 chars, and parse as double
    // Convert filename frequency into angular frequency, in rad/sec
    double w = std::stod(input_fname.substr(0,4)) * 1000. * 2.*M_PI;

    // Hard coded circuit parameters, we can change these to file parameters later
    double R = 9. + 220. + 470.;
    double L = 0.01;
    double C = 0.22e-6;

    // Calculate the estimated circuit response parameters
    double w0 = 1. / std::sqrt(L * C);
    double gamma = R / L;
    double V0 = 0.12;

    std::cout << "Expected: w0 = " << w0 << "\tgamma = " << gamma << "\n";
    // Read chan1, chan2 data
    read_data(input_fname, t1, y1, dy1, t2, y2, dy2);

    std::vector<double> init_vals = {{V0, w, 0., 0.}};

    // Fit the first set of parameters, aka CH_A, the input signal
    minimizer::nonlinear_ls<std::vector<double>> fitter(t1, y1, dy1, sin_f, sin_df, sin_fvv, init_vals); 
    fitter.fit();
    fitter.print_results();

    // Use the estimated response parameters as our initial guess
    // Multiply input voltage by the ratio V_c/V_in (see the "Useful fitting quantities" module)
    init_vals[0] *= w0 * w0 / std::sqrt(std::pow(w0 * w0 - w * w, 2) + gamma * gamma * w * w);
    // Now, add the phase difference to our fitted phase
    init_vals[2] = M_PI / 2. + std::atan((w * w - w0 * w0)/(w * gamma));

    // And now run the fitter for V_c!
    auto fitter2 = minimizer::nonlinear_ls<std::vector<double>>(t2, y2, dy2, sin_f, sin_df, sin_fvv, init_vals); 
    fitter2.fit();
    fitter2.print_results();

    // We're going to loop over the files in bash, so just append for now
    std::ofstream amp_file("rel_amplitudes.dat", std::ios::app);
    std::ofstream phase_file("phase_diff.dat", std::ios::app);

    double rel_amp = fitter2.parameter(0) / fitter.parameter(0);
    // Propagate uncertainties...
    double drel_amp = rel_amp * std::sqrt( std::pow(fitter.uncertainty(0)/fitter.parameter(0), 2)
                                         + std::pow(fitter2.uncertainty(0)/fitter2.parameter(0), 2)
                                         );

    double phase_diff = fitter2.parameter(2) - fitter.parameter(2);
    double dphase = std::sqrt(std::pow(fitter.uncertainty(2),2) + std::pow(fitter2.uncertainty(2), 2));

    /*
    if (rel_amp < 0)
    {
        // sin(x) = -sin(x +/- pi), so the fitter can sometimes try to fit a shifted-inverted sine.
        // in that case, we have to correct this (and recalculate the phase modulo [-pi/2, 3pi/2].
        rel_amp *= -1;
        drel_amp *= -1;

        phase_diff += M_PI;
    }
    */

    /*
    // We can optimize this later, just get it running for now
    while (phase_diff < -M_PI/2.)
    {
        phase_diff += 2*M_PI;
    }
    while (phase_diff > 1.5 * M_PI)
    {
        phase_diff -= 2*M_PI;
    }*/

    amp_file << w << " " << rel_amp << " " << drel_amp << "\n";
    phase_file << w << " " << phase_diff << " " << dphase << "\n";

    return 0;
}

void read_data(std::string fname, arrref t1, arrref y1, arrref dy1, arrref t2, arrref y2, arrref dy2)
{
    std::size_t line_count = 0u;
    std::ifstream infile(fname);
    if (!infile.is_open())
    {
        throw std::runtime_error("Error opening input file!\n");
    }
    for (std::string temp; std::getline(infile, temp);)
    {
        line_count++;
    }

    // Clears EOF bit as of C++11
    infile.clear();
    infile.seekg(0);

    for (auto i = 0u; i < line_count; i++)
    {
        double temp;
        infile >> temp;
        t1.push_back(temp);
        infile >> temp;
        y1.push_back(temp);
        infile >> temp;
        dy1.push_back(1/(temp*temp));
        
        infile >> temp;
        t2.push_back(temp);
        infile >> temp;
        y2.push_back(temp);
        infile >> temp;
        dy2.push_back(1/(temp*temp));
    }
}

int sin_f (const gsl_vector *x, void *data, gsl_vector *f)
{
    // y = A * sin(b*t + c) + D

    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    std::size_t n = sdata->n;
    double *t = sdata->t;
    double *y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);
    double D = gsl_vector_get(x, 3);

    for (auto i = 0u; i < n; i++)
    {
        double yi = A * std::sin(B * t[i] + C) + D;
        gsl_vector_set(f, i, yi - y[i]);
    }

    return GSL_SUCCESS;
}

int sin_df (const gsl_vector *x, void *data, gsl_matrix *J)
{
    // y = A * sin(B*t + C) + D

    // dy/dA = sin(B*t+c)
    // dy/dB = A * t * cos(B*t + C)
    // dy/dC = A * cos(B*t + C)
    // dy/dD = 1

    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    std::size_t n = sdata->n;
    double *t = sdata->t;
    // double *y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);
    // double D = gsl_vector_get(x, 3);

    for (auto i = 0u; i < n; i++)
    {
        gsl_matrix_set (J, i, 0, std::sin(B * t[i] + C));
        gsl_matrix_set (J, i, 1, A * t[i] * std::cos(B * t[i] + C));
        gsl_matrix_set (J, i, 2, A * std::cos(B * t[i] + C));
        gsl_matrix_set (J, i, 3, 1);
    }

    return GSL_SUCCESS;
}

int sin_fvv (const gsl_vector *x, const gsl_vector *v, void *data, gsl_vector *fvv)
{
    // y = A * sin(B*t + C) + D

    // First partials
    // dy/dA = sin(B*t+c)
    // dy/dB = A * t * cos(B*t + C)
    // dy/dC = A * cos(B*t + C)
    // dy/dD = 1

    // Second non-mixed partials
    // d^2y/dA^2 = 0
    // d^2y/dB^2 = - A * t^2 * sin(B*t+C)
    // d^2y/dC^2 = - A * sin(B*t+C)
    // d^2y/dD^2 = 0
    
    // Mixed partials (equal by clairaut, so double them)
    // d^2y/dAdB = t * cos(B*t+C)
    // d^2y/dAdC = cos(B*t+C)
    // d^2y/dAdD = 0
    // d^2y/dBdC = - A * t * sin(B*t + C)
    // d^2y/dBdD = 0
    // d^2y/dCdD = 0
    
    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    auto n = sdata->n;
    auto t = sdata->t;
    // auto y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);
    // double D = gsl_vector_get(x, 3);

    for (auto i = 0u; i < n; i++)
    {
        auto va = gsl_vector_get(v, 0);
        auto vb = gsl_vector_get(v, 1);
        auto vc = gsl_vector_get(v, 2);

        auto fab = t[i] * std::cos(B*t[i]+C);
        auto fac = std::cos(B*t[i]+C);
        auto fbb = - A * t[i] * t[i] * std::sin(B*t[i] + C);
        auto fbc = - A * t[i] * std::sin(B*t[i]+C);
        auto fcc = - A * std::sin(B*t[i]+C); 

        auto fvvi = vb * vb * fbb + vc * vc * fcc
            + 2. * (va * vb * fab + va * vc * fac + vb * vc * fbc);

        gsl_vector_set(fvv, i, fvvi);
    }
    return GSL_SUCCESS;
}

