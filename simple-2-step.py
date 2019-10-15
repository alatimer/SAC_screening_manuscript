import sys
import math
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cmx
from mpmath import mpf,mp

sys.path.append('./src/')
sys.path.append('./src/ase')
from ase.units import kB
from delta import Reactant,Reaction

def reader(data_file,ftype='separate',from_nmat=False):
    """ Reads in all data points from data file provided
    and returns a dictionary to plot on volcano
    """
    new_dict = {}
    for i,line in enumerate(open('./data/%s'%(data_file),'r').readlines()):
        if i == 0:
            labels = line.split()
        else:
            if line.startswith('#'):
                continue
            entries = line.split()
            if ftype == 'together':
                surftype = entries[0]
                surf = entries[1]
                if from_nmat==True:
                    surftype='nmat_pub'
                if surftype not in new_dict.keys():
                    new_dict[surftype] = {}
                #key = entries[0]+'_'+entries[1]
                new_dict[surftype][surf]={}
                energies = [float(eng) for eng in entries[2:]]
                for label,energy in zip(labels,energies):
                    if from_nmat==True and label=='EH':
                        energy-=0.5*H2O_form.get_dE() #change to O2/H2 ref
                    new_dict[surftype][surf][label]=energy
                new_dict[surftype][surf]['from_nmat']=from_nmat
            else:
                surf = entries[0]
                new_dict[surf] = {}
                energies = [float(eng) for eng in entries[1:]]
                for label,energy in zip(labels,energies):
                    new_dict[surf][label]=energy
    return new_dict

#setting reaction conditions. Temperature can be changed, but pressures are not functional.
annot_sites=True
kB=0.0000862
h=4.1356e-15
T=550
PCH4 = 1
PN2O = 1
PN2 = 1

#maximum and minimum logRates to plot
rmin = -30
rmax= 0

EO_vec = np.arange(-7,5,12./25.)
EH_vec = np.arange(-4,3,7/25.)

#Reading in atoms objects and using their free energy corrections
gasDB = pickle.load(open('gases.pkl','rb'))
gasDB = gasDB.filter(lambda x: x.calc_params['xc']=='BEEF')
gasDB = gasDB.filter(lambda x: x.calc_params['pw']=='650')
N2O = gasDB.filter(lambda x: x.surf_name=='N2O').data[0]
CH3 = gasDB.filter(lambda x: x.surf_name=='CH3').data[0]
H2 = gasDB.filter(lambda x: x.surf_name=='H2').data[0]
H2O = gasDB.filter(lambda x: x.surf_name=='H2O').data[0]
O2 = gasDB.filter(lambda x: x.surf_name=='O2').data[0]
O3 = gasDB.filter(lambda x: x.surf_name=='O3').data[0]
N2 = gasDB.filter(lambda x: x.surf_name=='N2').data[0]
CH3OH = gasDB.filter(lambda x: x.surf_name=='CH3OH').data[0]
CH4 = gasDB.filter(lambda x: x.surf_name=='CH4').data[0]
N2 = gasDB.filter(lambda x: x.surf_name=='N2').data[0]

#N2O activation from Pd on Gr
N2_O = Reactant.Reactant(
        species_type='adsorbate',
        traj_loc='./atoms_objects/N2-O/init.traj',
        vib_loc='./atoms_objects/N2-O/' )

#Others on TiO2
CH3_H_O = Reactant.Reactant(
        species_type='adsorbate',
        traj_loc='./atoms_objects/CH3-H-O/init.traj',
        vib_loc='./atoms_objects/CH3-H-O/' )

O_star = Reactant.Reactant(
        species_type='adsorbate',
        traj_loc='./atoms_objects/O_star/init.traj',
        vib_loc='./atoms_objects/O_star/' )

CH3OH_star = Reactant.Reactant(
        species_type='adsorbate',
        traj_loc='./atoms_objects/CH3OHc_TiO2/init.traj',
        vib_loc='./atoms_objects/CH3OHc_TiO2/' )

OH_star = Reactant.Reactant(
        species_type='adsorbate',
        traj_loc='./atoms_objects/OH_star/init.traj',
        vib_loc='./atoms_objects/OH_star/' )

star = Reactant.Reactant(
        species_type='slab',
        traj_loc='./atoms_objects/star/init.traj' )


#Define methanol partial pressure based on collector dE, in this case -1.6 for ~Al2O3
ECH3OH = -1.6
GCH3OH = ECH3OH + CH3OH_star.get_Gcorr(T=T) - CH3OH.get_Gcorr(T=T)
print ECH3OH,GCH3OH
PCH3OH = np.exp(GCH3OH/kB/T)
print PCH3OH

#ECH3g = 2.33 #DFT
#Calculating gas formation energies
N2O_form = Reaction.Reaction(ISs=[O2,N2],FSs=[N2O])
H2O_form = Reaction.Reaction(ISs=[],FSs=[H2O])
CH3OH_form = Reaction.Reaction(ISs=[],FSs=[CH3OH])
CH4_form = Reaction.Reaction(ISs=[],FSs=[CH4])
#CH3_form = Reaction.Reaction(ISs=[],FSs=[CH3])

#ECH3g = CH3_form.get_dG(T=T,P=[101325*PCH3])
ECH3OHg = CH3OH_form.get_dG(T=T,P=[101325*PCH3OH])

######### Reaction network considered #############
## N2O_g + * <-> O* + N2_g
## O* + CH4_g <-> CH3OH_g + *

#@Arvin, barrier scaling should be double checked
#First step scales as EO
m1 = mpf(0.4661); b1 = 1.236
#Second step scale as EH
m3 = mpf(0.75); b3 = 1.96

#Set decimal points to keep track of
mp.dps = 50

count=0

#Calculate rate of methanol production
def fun_RCH3OH(EO,EH,T):
    
    global count
    
    EO = mpf(EO)
    EH = mpf(EH)
    #Calc GO (still relative to O2g)
    GO = EO + O_star.get_Gcorr(T) - 0.5*O2.get_Gcorr(T,101325)
    
    dG1 = GO - N2O_form.get_dG(T,P=101325*PN2O) 
    dG2 = CH3OH_form.get_dG(T,P=101325*PCH3OH)  - GO - CH4_form.get_dG(T,P=101325*PCH4)

    #1
    Ea1 = EO*m1+b1
    Ga1 = Ea1+N2_O.get_Gcorr(T)-N2O.get_Gcorr(T,PN2O) 
    Ga1 = max(Ga1,dG1)
    Ga1 = max(Ga1,mpf(0))
    k1 = kB*T/h*mp.exp(-Ga1/kB/T)
    
    #2
    Ea2 = EH*m3+b3
    Ga2 = Ea2+CH3_H_O.get_Gcorr(T)-O_star.get_Gcorr(T)-CH4.get_Gcorr(T,PCH4)
    Ga2 = max(Ga2,dG2)
    Ga2 = max(Ga2,mpf(0))
    k2 = kB*T/h*mp.exp(-Ga2/kB/T)

    #calculate reverse rates using equilib constants, no bigger than kbT/h
    k_1=min(k1/mp.exp(-dG1/kB/T),kB*T/h)
    k_2=min(k2/mp.exp(-dG2/kB/T),kB*T/h)
    
    a = PN2O*k1 + PCH3OH*k_2
    b = PN2O*k1 + PN2*k_1 + PCH4*k2 + PCH3OH*k_2
    t_O = a/b
    t_s = 1 - t_O 
    
    r1 = t_s*PN2O*k1 - t_O*PN2*k_1
    r2 = t_O*PCH4*k2 - t_s*PCH3OH*k_2

    thetas = [t_O,t_s]
    rates = [r1,r2]

    #print 'EO, EH, Ga1, dG1, k1, Ga2, dG2, k2, t_O, t_s'
    #print '%4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f'%(EO,EH,Ga1,dG1,k1,Ga2,dG2,k2, t_O, t_s)

    return thetas,rates


####################### MAKE THETA AND ELEM RATE GRIDS ###############################################

tO_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))
ts_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))
R1_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))
R2_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))
rR1_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))
rR2_grid = np.empty((EO_vec.shape[0],EH_vec.shape[0]))

grids = {}
grids['R1'] = R1_grid
grids['R2'] = R2_grid
grids['rR1'] = rR1_grid
grids['rR2'] = rR2_grid
grids['t_s'] = ts_grid
grids['t_O'] = tO_grid

for i,EO in enumerate(EO_vec):
    for j, EH in enumerate(EH_vec):
        thetas,rates = fun_RCH3OH(EO,EH,T)
        
        tO_grid[j][i] = thetas[0]
        ts_grid[j][i] = thetas[1]
        
        #Before taking log, take negative of rates to visualize any reverse reactions. 
        #Each elem. step should have valid entries for either the forward or reverse rate (but not both)
        R1_grid[j][i] =  min(rmax,mp.log(max(10**rmin,rates[0]),b=10))
        R2_grid[j][i] =  min(rmax,mp.log(max(10**rmin,rates[1]),b=10))
        rR1_grid[j][i] =  min(rmax,mp.log(max(10**rmin,-rates[0]),b=10))
        rR2_grid[j][i] =  min(rmax,mp.log(max(10**rmin,-rates[1]),b=10))

####################### PLOT GRIDS (VOLCANOS) ###############################################

for grid in grids:
    plt.close()
    fig = plt.figure(1,figsize=(9,7))
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$E_H$ (eV)',fontsize=20)
    ax.set_xlabel(r'$E_O$ (eV)',fontsize=20)
    ax.set_ylim((-3.5,2.))
    ax.set_xlim((-6,4))
    if 't_' in grid:
        ax.set_title('%s, T=%i'%(grid,T))
        temp_plot =  plt.contourf(EO_vec,EH_vec,grids[grid],25,cmap='jet',vmin=0,vmax=1)
    else:
        ax.set_title('log(%s), T=%i'%(grid,T))
        temp_plot =  plt.contourf(EO_vec,EH_vec,grids[grid],25,cmap='jet',vmin=rmin,vmax=rmax)
    plt.colorbar(temp_plot)
    plt.tight_layout()
    plt.savefig('%s.pdf'%(grid))


##################### PLOT CH3OH PRODUCTION VOLCANO WITH DATA ###############################
plt.close()
fig = plt.figure(1,figsize=(9,7))
ax = fig.add_subplot(111)
ax.set_ylabel(r'$E_H$ (eV)',fontsize=20)
ax.set_xlabel(r'$E_O$ (eV)',fontsize=20)
ax.set_title('T=%i'%(T))
ax.set_ylim((-3.5,2.))
ax.set_xlim((-6,4))
temp_plot =  plt.contourf(EO_vec,EH_vec,R2_grid,25,cmap='jet')
clb = plt.colorbar(temp_plot)
clb.ax.set_title('log(Rate)')

# read in arvin's data, refd to H2, O2
surf_dicts = reader('arvin-all-new.dat',ftype='together')
nmat_dicts = reader('univ_rad.full.dat',ftype='together',from_nmat=True)
surf_dicts.update(nmat_dicts)

NUM_COLORS = len(surf_dicts)
cm = plt.cm.Paired
cNorm = colors.Normalize(vmin=0,vmax=NUM_COLORS)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

markers = ['o','^','v','s','*','>','<']
color_idx = 0
allEO=[]
allEH=[]
for surf_type in surf_dicts:
    color_idx+=1
    color = scalarMap.to_rgba(color_idx)
    fillstyle = 'full'
    if surf_type=='nmat_pub':
        color='grey'
    for i,surf in enumerate(surf_dicts[surf_type]):
        EO = surf_dicts[surf_type][surf]['EO']
        y = surf_dicts[surf_type][surf]['EH']
       # y-=0.5*h2o_form_dft #change to h2o H ref if desired
        if i!=0:
            leg_lab='_nolegend_'
        else:
            leg_lab=surf_type
        ax.plot(EO,y,color=color,marker='o',fillstyle=fillstyle,label=leg_lab,markersize=10,markeredgecolor='k')
        allEO.append(EO)
        allEH.append(y)
        if annot_sites == True and surf_type!='nmat_pub':
            ax.annotate(surf,[EO,y],ha='center',va='center',color = 'k',size=6)

p,C_p = np.polyfit(allEO,allEH,1,cov=True)
sigma_m, sigma_b = np.sqrt(np.diag(C_p))
m,b = p

plt.plot(EO_vec,m*EO_vec+b,'b',label='Best fit: EH=%4.2f*EO+%4.2f'%(m,b))

plt.legend(loc='upper right',fontsize=10)
plt.tight_layout()
plt.savefig('Rate-w-data.pdf')
