/*! @file minimizer.hpp
 *  @brief Weighted Nonlinear Least Squares (Levenberg-Marquardt with Geodesic Acceleration)
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @copyright GPLv3
 *  @date 2021-04-12
 */

#include <cmath>
#include <cstdio>

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

namespace minimizer
{
    enum trs
    {
        LM = 0,
        LMACCEL,
        DOGLEG,
        DDOGLEG,
        SUBSPACE2D
        // CGST // Only available for multilarge systems
    };

    enum scaling
    {
        MORE = 0,
        LEVENBERG = 1,
        MARQUARDT = 2
    };

    enum solver
    {
        QR = 0,
        CHOLESKY = 1,
        MCHOLESKY = 2,
        SVD = 3
    };

    enum fdtype
    {
        FWDIFF = 0,
        CTRDIFF = 1
    };

    struct user_data
    {
        std::size_t n;
        double *t;
        double *y;
    };

    template <class T>
    class nonlinear_ls
    {
        public:
            nonlinear_ls(T t, T y, T dy,
                    int (*uf)(const gsl_vector *, void *, gsl_vector *),
                    int (*udf)(const gsl_vector *, void *, gsl_matrix *),
                    int (*ufvv)(const gsl_vector *, const gsl_vector *, void *, gsl_vector *),
                    T init_vals)
                : n(t.size()), p(init_vals.size()), xtol(1e-8),
                gtol(std::pow(std::numeric_limits<double>::epsilon(), 1./3.)),
                ftol(0.), MAX_ITERS(1000)
            {
                type = gsl_multifit_nlinear_trust;
                fdf_params = gsl_multifit_nlinear_default_parameters();
                covar = gsl_matrix_alloc(p,p);

                data.n = n;
                data.t = t.data();
                data.y = y.data();
                x = gsl_vector_view_array (init_vals.data(), p);
                wts = gsl_vector_view_array(dy.data(), n);
                
                fdf.f = uf;
                fdf.df = udf;
                fdf.fvv = ufvv;
                fdf.n = n;
                fdf.p = p;
                fdf.params = &data;

                workspace = gsl_multifit_nlinear_alloc(type, &fdf_params, n, p);
                gsl_multifit_nlinear_winit(&x.vector, &wts.vector, &fdf, workspace);
                f = gsl_multifit_nlinear_residual(workspace);
                gsl_blas_ddot(f, f, &chisq);
            }

            int fit()
            {
                auto start = std::chrono::steady_clock::now();
                int result = gsl_multifit_nlinear_driver(MAX_ITERS, xtol, gtol, ftol, NULL, NULL, &info, workspace);
                auto end = std::chrono::steady_clock::now();
                last_time = end - start;


                J = gsl_multifit_nlinear_jac(workspace);
                gsl_multifit_nlinear_covar(J, 0., covar);

                gsl_blas_ddot(f, f, &chisq);


                return result;
            }

            void print_results()
            {
                std::cout << "Reduced chi-square = " << chisq/(n-p) << "\n";
                for (auto i = 0u; i < p; i++)
                {
                    std::cout << "Parameter " << i << " = " << gsl_vector_get(workspace->x, i) 
                        << " +/- " << sqrt(gsl_matrix_get(covar, i, i) * (chisq/(n-p))) << "\n";
                }
                std::cout << "Elapsed time = " << last_time.count() << " seconds.\n";
            }

            double reduced_chisq() const
            {
                return chisq/(n-p);
            }

            double parameter(std::size_t i) const
            {
                if (i >= p)
                {
                    throw std::runtime_error("Desired parameter index does not exist");
                }
                return gsl_vector_get(workspace->x, i);
            }

            double uncertainty(std::size_t i) const
            {
                if (i >= p)
                {
                    throw std::runtime_error("Desired uncertainty index does not exist");
                }
                return sqrt(gsl_matrix_get(covar, i, i) * reduced_chisq());
            }

            ~nonlinear_ls()
            {
                gsl_multifit_nlinear_free(workspace);
                gsl_matrix_free(covar);
            }

        private:

            int info;
            std::size_t n;
            std::size_t p;

            struct user_data data;
            const gsl_multifit_nlinear_type *type;
            gsl_multifit_nlinear_workspace *workspace;
            gsl_multifit_nlinear_fdf fdf;
            gsl_multifit_nlinear_parameters fdf_params;
            gsl_vector *f;
            gsl_matrix *J;
            gsl_matrix *covar;

            gsl_vector_view x;
            gsl_vector_view wts;

            double chisq;
            double xtol;
            double gtol;
            double ftol;
            std::size_t MAX_ITERS;

            std::chrono::duration<double> last_time;
    };
}
