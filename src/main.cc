
#include <string>
#include <iostream>

using std::string;
using std::cout, std::endl;

////////////////////////////////////////////////////////////////////////////////

int driver_premixed();
int driver_diffusion();
int driver_diffusion_table();
int driver_flamelet();

////////////////////////////////////////////////////////////////////////////////

int main(int nargs, char *argv[]) {

    int ireturn;

    if(nargs < 2) {
        ireturn = driver_premixed();
    }
    else {
        if (string(argv[1]) == "premixed")
            ireturn = driver_premixed();
        else if (string(argv[1]) == "diffusion")
            ireturn = driver_diffusion();
        else if (string(argv[1]) == "diffusion_table")
            ireturn = driver_diffusion_table();
        else if (string(argv[1]) == "flamelet")
            ireturn = driver_flamelet();
        else {
            cout << "\n\nERROR: invalid case type argument\n" << endl;
            ireturn = 1;
        }
    }

    return ireturn;
}
