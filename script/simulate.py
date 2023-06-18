from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import pandas as pd
import numpy as np
import pickle
import mdtraj as md
from mpi4py import MPI
from openmmplumed import PlumedForce
comm1 = MPI.COMM_SELF
comm2 = MPI.COMM_WORLD

platform = Platform.getPlatformByName('CUDA')

residues = pd.read_csv("residues.csv").set_index('one')

fasta = """MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGWGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGWGSASNAGSGSGFNGGFGSSMDSKSSGWGM""".replace('\n', '')

def adjust_terminals_HIS(r,pH,fasta):
    r.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
    r.loc['X','MW'] = r.loc[fasta[0],'MW'] + 2.
    r.loc['X','three'] = 'X'
    
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
    r.loc['Z','MW'] = r.loc[fasta[-1],'MW'] + 16.
    r.loc['Z','three'] = 'Z'
    
    return r.set_index('three')

residues = adjust_terminals_HIS(residues,7.4,fasta)

atomic_number = 117
for i,r in residues.iterrows():
    name = 'carbon-{}'.format(i)
    symbol = 'C-{}'.format(i)
    new_element = Element(atomic_number,name,symbol,r.MW*dalton)
    atomic_number += 1

pdb = PDBFile('input.pdb')

top = pdb.topology

for chain in top._chains:
    for residue in chain._residues:
        for atom in residue._atoms:
            atom.element = Element._elements_by_symbol['C-{}'.format(residue.name)]
            
atoms = list(top.atoms())

for i in range(len(fasta)-1):
    top.addBond(atoms[i],atoms[i+1])
    
forcefield = ForceField("calvados.xml")

system = forcefield.createSystem(top,nonbondedMethod=CutoffNonPeriodic)

system.removeForce(3)


with open('r1_excl.pkl', 'rb') as fp:
    r1_exclusions = pickle.load(fp)
    
forces = system.getForces()

forces[1].setCutoffDistance(4*nanometers)
for bond in r1_exclusions:
    forces[1].addExclusion(bond[0],bond[1])
    
forces[2].setForceGroup(1)
forces[2].setCutoffDistance(2*nanometers)
for bond in r1_exclusions:
    forces[2].addExclusion(bond[0],bond[1])

script = """MOLINFO MOLTYPE=protein STRUCTURE=input.pdb
WHOLEMOLECULES ENTITY0=1-414

distance_map:  CONTACTMAP ...
ATOMS1=1,401
ATOMS2=1,402
ATOMS3=1,403
ATOMS4=1,404
ATOMS5=1,405
ATOMS6=1,406
ATOMS7=1,407
ATOMS8=1,408
ATOMS9=1,409
ATOMS10=1,410
ATOMS11=1,411
ATOMS12=1,412
ATOMS13=1,413
SWITCH={CUSTOM FUNC=x R_0=1}
...

af_dist: CONSTANT VALUES=2.1632848754525185,2.160191326960921,2.153502482734621,2.1511839145794514,2.15038164164871,2.1478565398603684,2.142725614644587,2.140135044977069,2.133207448013127,2.1274217223748564,2.1308686502277854,2.116384238377213,2.115404123440385 NODERIV

rg: GYRATION TYPE=RADIUS ATOMS=1-414

PRINT FILE=COLVAR ARG=rg STRIDE=200
PRINT FILE=DISTANCE_MAP ARG=distance_map.* STRIDE=200

ENSEMBLE ...
    ARG=(distance_map.*) 
    LABEL=ens
... ENSEMBLE

STATS ...
    ARG=(ens\.distance_map.*) PARARG=(af_dist.*)
    LABEL=af_fw_comp
... STATS

PRINT ARG=af_fw_comp.*,(ens\.distance_map.*)      STRIDE=200 FILE=ST.AF_FW

METAINFERENCE ...
    ARG=(distance_map.*)
    PARARG=(af_dist.*)
    NOISETYPE=MGAUSS
    #OPTSIGMAMEAN=SEM
    AVERAGING=200
    SIGMA_MEAN0=1
    SIGMA0=50.0 SIGMA_MIN=0.0001 SIGMA_MAX=150.0 DSIGMA=0.1
    WRITE_STRIDE=10000
    LABEL=af_mi
    TEMP=298
... METAINFERENCE

FLUSH STRIDE=200
PRINT ARG=af_mi.*   STRIDE=200 FILE=BAYES.af
ENDPLUMED"""

system.addForce(PlumedForce(script, comm1, comm2))
                        
N_res = len(fasta)
N_save = 3000 if N_res < 100 else int(np.ceil(3e-4*N_res**2)*1000)
N_steps = 1100*N_save
DCD_output_file = "output.dcd"
stats_output_file = "stats.csv"

integrator = LangevinIntegrator(298*kelvin, 0.01/picosecond, 0.005*picoseconds)

simulation = Simulation(top, system, integrator, platform)
    
simulation.context.setPositions(pdb.positions)

#simulation.minimizeEnergy()

simulation.loadCheckpoint("checkpoint_CUDA")

simulation.reporters.append(DCDReporter(DCD_output_file, N_save))

simulation.reporters.append(StateDataReporter(stats_output_file, N_save, step=True, potentialEnergy=True, temperature=True))

simulation.step(100000000000)

