from lammps_input_writer import *
from calculate_mu import *
import ase, ase.io, ase.build
import os, sys, yaml
from warnings import simplefilter
simplefilter("ignore")

yamlfile=str(sys.argv[1])
with open(yamlfile, 'r') as f:
    yamlinput=yaml.safe_load(f)


### lmp should be the statically compiled version of lammps. In this example we attach a ACE potential.
# e.g. "mpirun -np 4 ./lmp_mpi" or "srun ./lmp_mpiicpc"
LAMMPS_RUN_COMMAND = f"srun ./lmp"


#run folder
folder=f'GAP_6x6x32_3node_mdstep_25k_dt_0.002_colderstart'
if not(os.path.exists(folder)):
    os.mkdir(folder)

#As defined in a LAMMPS input file: make sure to end each line with \n 
potential_str=("pair_style      hybrid/overlay pace table spline 6000 \n"
               "pair_coeff       * *  pace Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights.yace Ti \n pair_coeff * *  table Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights_pairpot.table Ti_Ti \n")


#pair_style = 'hybrid/overlay pace table spline 6000'
#pair_coeff = ' * *  pace Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights.yace Ti \n pair_coeff * *  table Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_customweights_pairpot.table Ti_Ti'


#building things
topdir=os.getcwd()
atoms=ase.io.read('bcc_Ti_geomopti.castep')
a = atoms.cell.cellpar()[0]*np.sqrt(4/3)
atoms = ase.build.bulk('Ti', crystalstructure='bcc', a=a, orthorhombic='True')
N_x=6
N_y=6
N_z=32
sc=[N_x,N_y,N_z]
atoms = atoms*sc
atoms.rattle(0.05)

pressures=yamlinput['pressures']
start_temp=yamlinput['start_temp']

##starting temp init: 1800K

#initialise interface pinning parameters
IP=IPparams(toploc=folder,
            potential=potential_str,
            pressure=pressures[0],
            temperature=start_temp,
            step=25000,
            step_setup=20000,
            hkl=[N_x,N_y,0],
            thermostat=0.35,
            barostat=0.75,
            dt=0.002,
	    Nz=N_z,
            thinned=50,
            atoms=atoms.copy(),
            LAMMPS_RUN_COMMAND=LAMMPS_RUN_COMMAND,
            auto=True)
IP.write_IP_params()


#run one IP routine.

#default is melt_steps=100, tol=10, samples=10
run_till_converged(IP, IP.traj_name, melt_steps=100, tol=15, samples=10)

#phase diagram via isobars
for i, p in enumerate(pressures[1:]):
    crystal_cc = ReadLammps(f'{IP.location}/crystal_auto.out')
    liquid_cc = ReadLammps(f'{IP.location}/liquid_auto.out')
    next_temp=classius_clapeyron_next_temp(current_T=IP.temperature, dP=p-pressures[i], crystal=crystal_cc, liquid=liquid_cc, thinned=IP.thinned)
    IP.next_isobar_start(pressure=p, temperature=next_temp)
    run_till_converged(IP, IP.traj_name, melt_steps=100, tol=15, samples=10)


