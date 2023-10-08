# OpenMM-Plumed with MPI Support

## Installation

```
conda create -n pl python=3.9
conda activate pl

conda install mpich
conda install mpi4py
conda install openmm -c conda-forge
conda install plumed=2.8.2=mpi_mpich_h7ded119_0 -c conda-forge
conda install py-plumed -c conda-forge
conda install cmake swig
conda install pandas
conda install mdtraj -c conda-forge

git clone https://github.com/husseinmur/OpenMM-Plumed-MPI

cd OpenMM-Plumed-MPI
mkdir build install openmm -p plumed/include -p plumed/lib
unzip openmm.zip -d openmm
unzip plumed_include.zip -d plumed/include
unzip plumed_lib.zip -d plumed/lib

cd build
ccmake ..
```

click c to configure and set (replace paths as required):
```
CMAKE_INSTALL_PREFIX             {git clone path}/OpenMM-Plumed-MPI/install
MPI4PY_DIR                       {user path}/miniconda3/envs/pl/lib/python3.9/site-packages/mpi4py/include
OPENMM_DIR                       {git clone path}/OpenMM-Plumed-MPI/openmm
PLUMED_INCLUDE_DIR               {git clone path}/OpenMM-Plumed-MPI/plumed/include
PLUMED_LIBRARY_DIR               {git clone path}/OpenMM-Plumed-MPI/plumed/lib
```
Keep everything else as is and click c then g then:
```
make
make install
make PythonInstall

cd ../install/lib
cp -r * ~/miniconda3/envs/pl/lib
```

## Running the simulation
In the script folder, run `simulate.py` as follows:

```
mpirun -np {number of replicas} python simulate.py
```

Note: you might need to replace the checkpoint, the forcefield xml, and other files if you need to run the simulation on a different system. The provided files are prepared for TDP-43.

To generate the forcefield xml and the exclusions pickle for a different system, you can use the `gen_xml_and_constraints.py` script in the scripts folder. It takes a fasta file as a parameter, and the file should only include a fasta sequence.
