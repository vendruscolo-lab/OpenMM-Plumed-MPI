%module openmmplumed


%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"
%include "mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);

%{
#include "PlumedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include <mpi.h>
%}

%pythoncode %{
import simtk.openmm as mm
%}

namespace PlumedPlugin {

class PlumedForce : public OpenMM::Force {
public:
    PlumedForce(const std::string& script, const MPI_Comm intra_comm, const MPI_Comm inter_comm);
    const std::string& getScript() const;
};

}
