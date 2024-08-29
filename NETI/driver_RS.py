import ase, ase.build, ase.io
import numpy
from lammps_python_interface import ThermoIntParams
from free_energy_getter import get_free_energy, get_RS_path
import os, sys, yaml
import numpy as np



#lmp is the statically compiled version of lammps. We use a ACE potential in this example.
#LAMMPS_RUN_COMMAND = "mpirun -np 16 ./lmp_g++_openmpi"
LAMMPS_RUN_COMMAND = f"srun ./lmp"


yamlfile=str(sys.argv[1])
with open(yamlfile, 'r') as f:
    yamlinput=yaml.safe_load(f)





#hcp=ase.io.read('data.Ti_hcp_ortho', format='lammps-data', style='atomic') #N=4
#bcc=ase.io.read('data.Ti_bcc_ortho', format='lammps-data', style='atomic') #N=2
#hex=ase.io.read('data.Ti_hex_ortho', format='lammps-data', style='atomic') #N=6



atoms=ase.io.read(yamlinput['primitive_cell'], format='lammps-data', style='atomic')

#atoms=ase.build.bulk('Ti', orthorhombic=True, crystalstructure='bcc', a=3.2)
sc=yamlinput['supercell']
atoms=atoms*sc
atoms_unrattled=atoms.copy()
atoms.rattle(0.05)


inittemp=yamlinput['start_temp']
finaltemp=yamlinput['final_temp']

Trange=np.linspace(inittemp,finaltemp,51)
#pair_style = "quip"
#pair_coeff = "* * Ti_turbo_awfr_0.02_min_0.025_cwvr_0.05_min_0.10_SP3000_curpoints.xml \"Potential xml_label=GAP_2023_8_10_60_16_29_2_871\" 22"
pair_style = 'hybrid/overlay pace table spline 6000'
pair_coeff = ' * *  pace Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights.yace Ti \n pair_coeff * *  table Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights_pairpot.table Ti_Ti'


#pair_style = "eam/fs"
#pair_coeff = "* * Ti1.eam.fs Ti"
name=yamlinput['name']


aniso=True
if yamlinput['iso']:
    anisostring='iso'
    aniso = False
else:
    anisostring='aniso'

folder=f"RS-{name}-{sc[0]}x{sc[1]}x{sc[2]}-{Trange[0]}-{Trange[-1]}-{anisostring}"
if not(os.path.exists(folder)):
    os.mkdir(folder)

TIP=ThermoIntParams(
    LAMMPS_RUN_COMMAND=LAMMPS_RUN_COMMAND, 
    atoms=atoms,
    atoms_unrattled=atoms_unrattled, 
    toploc=folder, 
    pair_style=pair_style,
    pair_coeff=pair_coeff,
    mass=48,
    pressure=float(yamlinput['pressure']), 
    temperature=Trange[0], 
    timestep=0.002, 
    nstep=25000,
    nstep_setup=25000,
    nstep_eq=25000,  
    thermostat=0.95, 
    barostat=2.00,
    averaging_setup=500,
    averaging=1000,
    thermoprint=20,
    aniso=aniso
)


get_free_energy(TIP, fl_file=f'thermoint_{folder}-P{TIP.pressure:.1f}.fl')
get_RS_path(TIP, finaltemp=finaltemp, fl_file=f'thermoint_{folder}-P{TIP.pressure:.1f}.fl', rs_file=f'thermoint_{folder}-P{TIP.pressure:.1f}.rs')

#for T in Trange:
#    TIP.update_PT(new_pressure=0, new_temperature=T)
#    get_free_energy(TIP, logfilename=f'thermoint_{folder}.log')
