#pragma once

#include <vector>
#include <cstdlib>
#include <iostream>

//-----------------------------

class linearInterp {

    public:

        std::vector<double> *X;
        std::vector<double> *Y;

        int nxy;

    private:

        int ilo;
        int ihi;

        //////////////////////////////////////////////////////////////////

    public:

        double interp(double x){
            set_bounding_indicies(x);
            return (*Y)[ilo] + (x-(*X)[ilo])*((*Y)[ihi]-(*Y)[ilo])/((*X)[ihi]-(*X)[ilo]);
        }

    private: 

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

        linearInterp(){}

        linearInterp(std::vector<double> &X_p, std::vector<double> &Y_p){
            X = &X_p;
            Y = &Y_p;
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

