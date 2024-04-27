#pragma once

#include <vector>
#include <cstdlib>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
///
/// Linear interpolation class (one-dimensional)
/// Header-only class
///
///////////////////////////////////////////////////////////////////////////////

class linearInterp {

    public:

        std::vector<double> *X;         ///< abscissas for Y(X)
        std::vector<double> *Y;         ///< ordinate values for Y(X)

        int nxy;                        ///< number of grid points

    private:

        int ilo;                        ///< lower of adjacent bounding grid points for desired interpolation point
        int ihi;                        ///< upper of adjacent bounding grid points for desired interpolation point

        //////////////////////////////////////////////////////////////////

    public:

        /// Interpolation function for arbitrary point x.
        /// @param x \input grid value we are interpolating to (abscissa)
        /// @return interpolated value y(x)
        double interp(double x){
            set_bounding_indicies(x);
            return (*Y)[ilo] + (x-(*X)[ilo])*((*Y)[ihi]-(*Y)[ilo])/((*X)[ihi]-(*X)[ilo]);
        }

    private: 

        /// find values of ilo and ihi using a search algorithm from stl
        /// @param x \input grid value we are interpolating to (abscissa)
        void set_bounding_indicies(double x){
            if(x <= (*X)[0])
                ihi = 1;
            else if(x >= (*X).back())
                ihi = nxy-1;
            else {
                std::vector<double>::iterator itHi = std::lower_bound((*X).begin(), (*X).end(), x); // lower_bound gives values >= x
                ihi = itHi - (*X).begin();
            }
            ilo = ihi-1;
        }

        //////////////////////////////////////////////////////////////////

    public:

        //------------- constructors

        /// Default (empty) constructor
        linearInterp(){}

        /// constructor function
        /// @param _X grid of abscissas forming the basis for interpolation
        /// @param _Y grid of ordinate values _Y(_X) forming the basis for interpolation
        linearInterp(std::vector<double> &_X, std::vector<double> &_Y){
            X = &_X;
            Y = &_Y;
            nxy = (*X).size();
            if((*X).size()!=(*Y).size()){
                std::cout << std::endl << "Error in interp_linear: X, Y need same size" << std::endl;
                exit(0);
            }

            for(int i=1; i<nxy-1; i++) 
                if( ((*X)[i] == (*X)[i-1]) || 
                    ((*X)[i] == (*X)[i+1]) || 
                    ((*X)[i]<(*X)[i-1] && (*X)[i]<(*X)[i+1]) ||
                    ((*X)[i]>(*X)[i-1] && (*X)[i]>(*X)[i+1]) ){
                    std::cout << std::endl << "Error in linearInterp: X should be monotonic" << std::endl;
                    exit(0);
                }
        }
};

