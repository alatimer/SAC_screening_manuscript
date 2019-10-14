#!/home/vossj/suncat/bin/python

#SBATCH -p iric
#SBATCH --job-name=adsFreq.py
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=48:00:00                                 #default is 20 hours
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL                            #get emailed about job BEGIN, END, or FAIL
#SBATCH  --mail-user=allegralatimer@gmail.com
#SBATCH --ntasks-per-node=16                            #task to run per node; each node has 16 cores

from ase.io import read
from ase.constraints import FixAtoms
from ase.vibrations import Vibrations
from espresso import espresso
from espresso.vibespresso import vibespresso
from ase.thermochemistry import HarmonicThermo
import os

#############
## more information here: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
#############

xc = 'BEEF'
kpts=(6,6,1)
pw = 550
dw = 5500
psppath = '/home/alatimer/bin/psp/esp_psp_w_gbrvRu/'
dipole = {'status':True}
spinpol = True
output = {'removesave':True}
outdir = 'calcdir'
sigma = 0.1
smearing = 'fd' #default

vibrateatoms=[48,49,50,51,52,53]	      	      	      	       # calculate the vibration modes of atoms #12 and #13

if os.path.exists('qn.traj')==True  and os.stat("qn.traj").st_size > 0:
    atoms =read('qn.traj')
else:
    atoms = read('init.traj')

if spinpol == True:
    for atom in atoms:
        if atom.symbol != 'H':
            atom.magmom = 3

# special calculator for the vibration calculations
calcvib = vibespresso(pw = pw,
                dw = dw,
                kpts = kpts, 
                nbands = -100,
                xc = xc, 
                psppath=psppath,
                convergence = {'energy':1e-5,
                               'mixing':0.1,
                               'nmix':10,
                               'maxsteps':200,
                               'diag':'david'
                                },
                spinpol = spinpol,
                smearing = smearing,
                sigma = sigma,
                outdirprefix = 'vibdir',
                )  	      	      	      	   # log file                                         


#energy = atoms.get_potential_energy()                      # caclulate the energy, to be used to determine G

atoms.set_calculator(calcvib)    	      	      	       # attach vibrations calculator to the atoms                   

# Calculate vibrations                                                                                        
vib = Vibrations(atoms,indices=vibrateatoms,delta=0.03)    # define a vibration calculation                   
vib.run()     	      	      	      	      	      	   # run the vibration calculation                    
vib.summary(method='standard')	      	      	      	   # summarize the calculated results                 

for mode in range(len(vibrateatoms)*3):                    # Make trajectory files to visualize the modes.    
    vib.write_mode(mode)

# Calculate free energy
vibenergies=vib.get_energies()
#gibbs = HarmonicThermo(vib_energies = vibenergies, potentialenergy = energy)

# At 300K and 101325 Pa
# change for your operating conditions 
#freeenergy = gibbs.get_free_energy(300,101325)

#f=open('out.energy','w')
#f.write('Potential energy: '+str(energy)+'\n'+'Free energy: '+str(freeenergy)+'\n')
#f.close
