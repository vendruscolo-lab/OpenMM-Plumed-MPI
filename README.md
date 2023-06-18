# How to install OpenMM-Plumed-MPI

#### NOTE: if using anaconda replace miniconda3 with your anaconda folder.

```
conda create -n pl python=3.9
conda activate pl

conda install mpich
conda install mpi4py
conda install openmm -c conda-forge
conda install plumed=2.8.2=mpi_mpich_h7ded119_0 -c conda-forge
conda install py-plumed -c conda-forge
conda install cmake swig

cd /pool/work/faidon/
git clone https://github.com/husseinmur/OpenMM-Plumed-MPI

cd OpenMM-Plumed-MPI
mkdir build install openmm -p plumed/include -p plumed/lib
unzip openmm.zip -d openmm
unzip plumed_include.zip -d plumed/include
unzip plumed_lib.zip -d plumed/lib

cd build
ccmake ..
```

click c to configure and set:
```
CMAKE_INSTALL_PREFIX             /pool/work/faidon/OpenMM-Plumed-MPI/install
CUDA_USE_STATIC_CUDA_RUNTIME     ON
MPI4PY_DIR                       /home/faidon/miniconda3/envs/pl/lib/python3.9/site-packages/mpi4py/include
OPENMM_DIR                       /pool/work/faidon/OpenMM-Plumed-MPI/openmm
PLUMED_BUILD_CUDA_LIB            ON
PLUMED_BUILD_OPENCL_LIB          ON
PLUMED_INCLUDE_DIR               /pool/work/faidon/OpenMM-Plumed-MPI/plumed/include
PLUMED_LIBRARY_DIR               /pool/work/faidon/OpenMM-Plumed-MPI/plumed/lib
```
Keep everything else as is and click c then g then:
```
make
make install
make PythonInstall

cd ../install/lib
cp -r * ~/miniconda3/envs/pl/lib
```
DONE
