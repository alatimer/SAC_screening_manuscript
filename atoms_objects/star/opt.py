#!/home/vossj/suncat/bin/python

#SBATCH -p iric,owners
#SBATCH --job-name=opt.py
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=48:00:00									#default is 20 hours
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL							#get emailed about job BEGIN, END, or FAIL
#SBATCH  --mail-user=allegralatimer@gmail.com
#SBATCH --ntasks-per-node=16 							#task to run per node; each node has 16 cores

from ase.constraints import *
from ase.utils.geometry import *
from ase.lattice.spacegroup import *
from math import *
from ase.lattice.surface import *
from ase import *
from ase.io import read,write,Trajectory
from ase.optimize import QuasiNewton
from espresso import espresso
from ase.dft.bee import BEEF_Ensemble
import cPickle as pickle
import os.path
import os
import pdos,bader

xc = 'BEEF'
kpts=(6,6,1)
pw = 550
dw = 5500
psppath = '/home/alatimer/bin/psp/esp_psp_w_gbrvRu/'
dipole = {'status':True}
spinpol = True
output = {'removesave':True}
outdir = 'calcdir'

save_pdos_pkl = True
save_cube = False
save_cd = False

if os.path.exists('qn.traj')==True  and os.stat("qn.traj").st_size > 0:
    atoms =read('qn.traj')
else:
    atoms = read('init.traj')
    ##atoms = atoms.repeat([1,1,2])
    ##atoms.cell[2][2] += 15

if spinpol == True:
    for atom in atoms:
        if atom.symbol != 'H':
            atom.magmom = 3
        else:
            atom.magmom = 0
else:
    for atom in atoms:
        atom.magmom = 0



######################### Calculator ##########################

calc = espresso(pw=pw,					#plane-wave cutoff
                dw=dw,			 #density cutoff
                xc=xc,		#exchange-correlation functional
                kpts=kpts,       #k-point sampling;
                nbands=-100,     #extra bands besides the bands needed to hold valence
                sigma=0.1,		#default value, fd smearing
		          spinpol=spinpol,
                psppath=psppath,
                convergence= {'energy':1e-5,
							  'mixing':0.05,
						     'nmix':10,
							  'mix':4,
							  'maxsteps':300,  #Change to 500 for spin polarized
							  'diag':'david'	  #Can change to 'cg' for spin polarized
										 },	#convergence parameters
                output = output,
                dipole=dipole, #dipole correction to account for periodicity in z
                outdir=outdir #output directory for Quantum Espresso files
                )	

if xc=='BEEF':
    calc.beefensemble = True
    calc.printensemble = True

####################   Constraints   ####################

# highest_surf_index = 47 
# z_list = [atom.z for atom in atoms if atom.index <= highest_surf_index]
# thresh = sum(z_list)/len(z_list)
# mask = [atom.z<thresh for atom in atoms]
#highest_index_to_fix = 23
#mask = [atom.index<=highest_index_to_fix for atom in atoms]  #use if manually constraining 

# fixatoms = FixAtoms(mask=mask)
# constraints = [fixatoms]
#atoms.set_constraint(constraints)

#Uncomment if need to fix bond lengths during opt
#Oa = 48
#Ob = 49
#H = 54
#CH3 = 50
#fixbondlengths = FixBondLengths([[Oa,H],[Oa,CH3],[CH3,H]])
#constraints.append(fixbondlengths)

# atoms.set_constraint(constraints)
write('pre.traj',atoms)				#Can check all constraints after calc begins to run in pre.traj

####################### Optimization ###################
atoms.set_calculator(calc)
qn = QuasiNewton(atoms, trajectory='qn.traj', logfile='qn.log')
qn.run(fmax=0.05)

#################### BEEF Error #########################
if xc == 'BEEF':
    ens = BEEF_Ensemble(calc)
    ens_e = ens.get_ensemble_energies()
    ens.write('ensemble.bee')
    pickle.dump(ens_e,open('ensemble.pkl','w'))

####################### PDOS ###################

pdos_traj = Trajectory('pdos_qn.traj','w',atoms)
pdos.pdos(atoms,outdir=outdir,save_pkl=save_pdos_pkl,spinpol=spinpol)
pdos_traj.write() #final image will have charges and duplicate geometry

bader_traj = Trajectory('bader_qn.traj','w',atoms)
bader.bader(atoms,outdir=outdir,save_cube=save_cube,save_cd=save_cd,spinpol=spinpol)
bader_traj.write()
