/** \name   interface_fortran.cpp
    \brief  Wrapper to some of AGAMA library's functions with added functionality
    \author Thor Tepper-Garcá
    \date   12 OCT 2025

This file works as a wrapper form some C++ Agama functions with C linkage, and provides additional functionality: it allows to delete a potential object created by invoking the AGAMA function agama_initfromfile_(), embedded in agama_initfromfile_wrapper_

This file is inteded to be used in conjunction with agama_interface.f90.
These two files demonstrate how to add functionality to an existing AGAMA installation without the need to modify its source and recompile.

**/

#include <cstring>
#include <iostream> /// only needed if debug message is used
#include "potential_factory.h"
#include "utils_config.h"
#include "utils.h"


// My own shared vector
std::vector<const potential::BasePotential*>& get_tracked_potentials() {
    static std::vector<const potential::BasePotential*> potentials;
    return potentials;
}

// Provided by ChatGPT (https://chatgpt.com/c/68d62367-ebf8-8320-8622-0f73bb0af83e)
// Wrapper to call original AGAMA init
extern "C" void agama_initfromfile_wrapper_(int verbose, void* c_obj, char* inifilename, long /*len(c_obj)=8*/, long len) {
    // Call the original AGAMA function
    extern void agama_initfromfile_(void* c_obj, char* inifilename, long /*len(c_obj)=8*/, long len);
    agama_initfromfile_(c_obj, inifilename, 8, len);

    // Track the pointer in your own vector
    const potential::BasePotential* ptr;
    memcpy(&ptr, c_obj, sizeof(void*));
    auto& potentials = get_tracked_potentials();
    potentials.push_back(ptr);

    if(verbose == 1){
        printf("[AGAMA DEBUG] Created potential at %p, tracked.size()=%zu\n", ptr, potentials.size());
        fflush(stdout);
    }
}

/// Provide by E.Vasiliev (2025)
/// delete a previously created potential stored in c_obj; obviously, it should not be used afterwards
extern "C" void agama_delete_(int verbose, void* c_obj, long)
{
    auto& potentials = get_tracked_potentials();
    const potential::BasePotential* pot;
    memcpy(&pot, c_obj, sizeof(void*));

    if(verbose == 1){
        printf("[AGAMA DEBUG] agama_delete_ called, potentials.size() = %zu\n", potentials.size());
        fflush(stdout);  // ensure immediate output
    }

//    for(size_t i=0; i<potentials.size(); i++)
//        if(potentials[i].get() == pot)
//            potentials[i].reset();

// ChatGPT (https://chatgpt.com/c/68d62367-ebf8-8320-8622-0f73bb0af83e) suggested this instead:
    for(size_t i = 0; i < potentials.size(); ++i) {
        if(potentials[i] == pot) {
            potentials.erase(potentials.begin() + i);  // removes from vector and deletes object
            if(verbose == 1){
                printf("[AGAMA DEBUG] Deleted potential at %p\n", pot);
                fflush(stdout);  // ensure immediate output
            }
            break;  // stop searching after deleting
        }
    }
    // zero out Fortran handle
    memset(c_obj, 0, sizeof(void*));
}
