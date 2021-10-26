#- Simulating the Rn220 decay chain approximately using discrete time steps --#
#- Developed by Chris Nedlik ------------------------------------------------------------------#

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
#                                             (64%) Po-212                      #
#                                            /            \                     #
#    Rn-220 -> Po-216 -> Pb-212 -> Bi-212 ->                -> Pb-208 (stable)  #
#                                            \            /                     #
#                                             (36%) Tl-208                      #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import matplotlib.pyplot as plt
import math
import time as stopwatch
t0 = stopwatch.time()

#------------------------------Settings--------------------------------#

A_0   = 6.0              # Rn220 (Pb212) Atoms Passing the Getter (Entering the Detector)
time  = 1.0              # Total time to simulate (hours)
dt    = 0.2              # Time step - seconds. Small for better approximation
inj_L = 1.0              # Length of injection in hours (only used when Continuous is True)

ER_Background = 19.06     # Events/keV in 60d between 1.5-6.5 keV
pur_tau = 2.36            # Purification by re-circulation 'half-life' (10,000 kg LXe at 500 SLPM) (used if Purification == True)
#pur_tau = 0.0424         # Purification by re-circulation 'half-life' (    90 kg GXe at 250 SLPM)
pur_dt  = 1.0             # Purification time step in seconds (if Purification = True)

Pb212_only_mode  = False   # If True, Rn220 and Po216 are not shown (don't make it into the TPC)
Purification     = False   # If True, isotopes are removed according to LZ purification time, t_1/2 = 2.35 days 
Continuous       = True   # If True, the Rn220 rate will stay fixed at the injection value
Cryopumping      = False   # If True, performs additional injection 24 hours after the first one, to account for a cryopumping period
Charge_Fractions = False  # If True, charged Pb212 atoms will be removed from decay counting and daughter production
Gas_only_mode    = False  # If True, only alpha decays are shown in the plots

#-----------------------------Constants--------------------------------#

#-- Decay Constants --#
Rn220_C = math.log(2)/55.6                # Halflife =  55.6 seconds
Po216_C = math.log(2)/0.145               # Halflife = 145.0 milliseconds
Pb212_C = math.log(2)/(10.64*60.0*60.0)   # Halflife = 10.64 hours
Bi212_C = math.log(2)/(60.55*60.0)        # Halflife = 60.55 minutes
Po212_C = math.log(2)/(0.000000299)       # Halflife = 0.299 microseconds
Tl208_C = math.log(2)/(3.053*60.0)        # Halflife = 3.053 minutes

alpha_ion_frac = 0.503  # Charged fraction of alpha daughters (from EXO-200)
beta_ion_frac  = 0.764  # Charged fraction of beta daughters (from EXO-200)

Bi_Po = 0.6406  # Bi212->Po212 Branching Fraction
Bi_Tl = 0.3594  # Bi212->Tl208 Branching Fraction


#-------------------Defining Convienient Variables----------------------#

alpha_neutral_frac = 1.0 - alpha_ion_frac # Neutral fraction of alpha daughters (from EXO-200)
beta_neutral_frac  = 1.0 - beta_ion_frac  # Neutral fraction of beta daughters (from EXO-200)

secs  = int(60.0*60.0*time)                # Run time in seconds
steps = int(secs/dt)                       # Number of dt steps
pur_frac = pur_dt/(pur_tau*24.0*60.0*60.0) # fraction of isotopes removed each pur_dt 


#-----------------------Initializing Arrays-----------------------------#

t = np.zeros(steps)  # time array

#-- Number of Atoms --#  #--- Decayed Atoms ---#   #-- Cumulative Events --#
Rn220 = np.zeros(steps); dRn220 = np.zeros(steps); Rn220C = np.zeros(steps);
Po216 = np.zeros(steps); dPo216 = np.zeros(steps); Po216C = np.zeros(steps);
Pb212 = np.zeros(steps); dPb212 = np.zeros(steps); Pb212C = np.zeros(steps);
Bi212 = np.zeros(steps); dBi212 = np.zeros(steps); Bi212C = np.zeros(steps);
Po212 = np.zeros(steps); dPo212 = np.zeros(steps); Po212C = np.zeros(steps);
Tl208 = np.zeros(steps); dTl208 = np.zeros(steps); Tl208C = np.zeros(steps);
Pb208 = np.zeros(steps); 

Rn220[0]  = A_0           # Number of atoms 
dRn220[0] = A_0*Rn220_C   # Initial dN/dt


#--------------------------- Main Loop ---------------------------------#

for i in range(1,steps):   
    t[i]   = dt*i
      
    #---------------- Rn220 -----------------#
    dRn220_tmp = (Rn220_C*Rn220[i-1]*dt)  # Number of decaying atoms in dt
    if(dRn220_tmp > Rn220[i-1]):          # If the # decaying > # remaining
        Rn220[i]  = 0                       # Then the # of atoms remaining is 0 (not negative)
        dRn220[i] = Rn220[i-1]              # The # of decayed atoms = the # that were left
    else:
        Rn220[i]  = Rn220[i-1]-dRn220_tmp
        dRn220[i] = dRn220_tmp
    Rn220C[i] = Rn220C[i-1]+dRn220[i]

    if Continuous==True:
        if Cryopumping==True:
            full   = ((t[i]<3600.0*inj_L)) | \
                     ((t[i]>3600.0*24.0)  & (t[i]<3600.0*(24.0+inj_L)))  
            if(full == True): 
                Rn220[i] = Rn220[i]+(A_0*dt)
                
        elif(t[i]<3600.0*inj_L):
            Rn220[i] = Rn220[i]+(A_0*dt)
            
    #---------------- Po216 ----------------#
    dPo216_tmp = (Po216_C*Po216[i-1]*dt)
    if(dPo216_tmp > Po216[i-1]):
        Po216[i]  = 0
        dPo216[i] = Po216[i-1]
    else:
        Po216[i]  = Po216[i-1]+dRn220[i]-dPo216_tmp
        dPo216[i] = dPo216_tmp
    Po216C[i] = Po216C[i-1]+dPo216[i]
    
    #---------------- Pb212 ----------------#
    dPb212_tmp = (Pb212_C*Pb212[i-1]*dt)
    if(dPb212_tmp > Pb212[i-1]):
        Pb212[i]  = 0
        dPb212[i] = Pb212[i-1]
    else:
        if(Charge_Fractions == True):
            Pb212[i] = Pb212[i-1]+(dPo216[i]*alpha_neutral_frac)-dPb212_tmp
        else:
            Pb212[i] = Pb212[i-1]+dPo216[i]-dPb212_tmp
        dPb212[i] = dPb212_tmp
    Pb212C[i] = Pb212C[i-1]+dPb212[i]
        
    #---------------- Bi212 ----------------#
    dBi212_tmp = (Bi212_C*Bi212[i-1]*dt)
    if(dBi212_tmp > Bi212[i-1]):
        Bi212[i]  = 0
        dBi212[i] = Bi212[i-1]
    else:
        Bi212[i] = Bi212[i-1]+dPb212[i]-dBi212_tmp
        dBi212[i] = dBi212_tmp
    Bi212C[i] = Bi212C[i-1]+dBi212[i]

    #---------------- Po212 ----------------#
    # These decay instantly in the time-steps of this sim, so never accumulate
    dPo212[i] = dBi212[i]*Bi_Po
    Po212[i]  = 0
    Po212C[i] = Po212C[i-1]+dPo212[i]

    #---------------- Tl208 ----------------#
    dTl208_tmp = (Tl208_C*Tl208[i-1]*dt)
    if(dTl208_tmp > Tl208[i-1]):
        Tl208[i]  = 0
        dTl208[i] = Tl208[i-1]
    else:
        Tl208[i]  = Tl208[i-1]+(dBi212[i]*Bi_Tl)-dTl208_tmp
        dTl208[i] = dTl208_tmp
    Tl208C[i] = Tl208C[i-1]+dTl208[i]

    #---------------- Pb208 ----------------#
    Pb208[i] = Pb208[i-1]+dPo212[i]+dTl208[i]

    #------------- Purification ------------#
    if Purification == True:
        if(t[i]%pur_dt == 0):
            Rn220[i] = Rn220[i]-(Rn220[i]*pur_frac)
            Po216[i] = Po216[i]-(Po216[i]*pur_frac)
            Pb212[i] = Pb212[i]-(Pb212[i]*pur_frac)
            Bi212[i] = Bi212[i]-(Bi212[i]*pur_frac)
            Po212[i] = Po212[i]-(Po212[i]*pur_frac)
            Pb208[i] = Pb208[i]-(Pb208[i]*pur_frac)


#------------------------- Final Analysis ------------------------------#
            
# Total Activities
Total = (dRn220[0::2]/dt)+(dPo216[0::2]/dt)+(dPb212[0::2]/dt)+(dBi212[0::2]/dt)+(dPo212[1::2]/dt)+(dTl208[0::2]/dt) 
Pb212Total = (dPb212[0::2]/dt)+(dBi212[0::2]/dt)+(dPo212[1::2]/dt)+(dTl208[0::2]/dt)   # Total Activity w/o Rn220, Po216
Gas_Total = (dRn220[0::2]/dt)+(dPo216[0::2]/dt)+((dBi212[0::2]/dt)*Bi_Tl)+(dPo212[1::2]/dt)


# Total Cumulative Events
TotalC = Rn220C[0::2]+Po216C[0::2]+Pb212C[0::2]+Bi212C[0::2]+Po212C[1::2]+Tl208C[0::2] 
Pb212TotalC = Pb212C[0::2]+Bi212C[0::2]+Po212C[1::2]+Tl208C[0::2] # Total cumulative events w/o Rn220,Po216

# Total Atoms
TotalA = Rn220[0::2]+Po216[0::2]+Pb212[0::2]+Bi212[0::2]+Po212[1::2]+Tl208[0::2]+Pb208[0::2]
Pb212TotalA = Pb212[0::2]+Bi212[0::2]+Po212[1::2]+Tl208[0::2]+Pb208[0::2] # Total Atoms w/o Rn220,Po216

print("Computation Time = "+str(round(stopwatch.time()-t0,1))+" sec")

X = 20  # must be even

T_X = int(X/2.0)
t=t/(60.0*60.0*24.0)

BG = (ER_Background/(60.0*24.0*60.0*60.0*0.002))*1.25 # Pb212 Equivalent BG Rate
BG_Counts = (ER_Background*5*1.25)/0.002

#---------------------------- Plotting ---------------------------------#

#---------------- Activity Plot -----------------#
plt.figure(1)
if(Pb212_only_mode==True):
    plt.plot(t[0::X], Pb212Total[0::T_X], 'grey', label='Total (No Rn220, Po216)')
    plt.plot(t[0::X], dPb212[0::X]/dt,'k',label='Pb212 (beta Q=569.1 keV)', linewidth=4)
    plt.plot(t[0::X], dBi212[0::X]/dt,'c',label='Bi212 (beta[64%] Q=2.3 MeV; alpha[36%] 6.2 MeV)')
    plt.plot(t[0::X], (dPo212[0::X]/dt),'b',label='Po212 (alpha 9.0 MeV)')
    plt.plot(t[0::X], dTl208[0::X]/dt,'y',label='Tl208 (beta Q=5.0 MeV)')

elif(Gas_only_mode==True):
    plt.plot(t[0::X], Gas_Total[0::T_X], 'grey', label='Total (Alphas Only')
    plt.plot(t[0::X], dRn220[0::X]/dt,'k',label='Rn220 (alpha 6.4 MeV)', linewidth="3")
    plt.plot(t[0::X], (dPo216[0::X]/dt),'g',label='Po216 (alpha 6.9 MeV)')
    plt.plot(t[0::X], (dBi212[0::X]/dt)*Bi_Tl,'c',label='Bi212 (alpha[36%] 6.2 MeV)')
    plt.plot(t[0::X], (dPo212[0::X]/dt),'b',label='Po212 (alpha 9.0 MeV)')

else:   
    plt.plot(t[0::X], Total[0::T_X], 'grey', label='Total')
    plt.plot(t[0::X], dRn220[0::X]/dt,'r',label='Rn220 (alpha 6.4 MeV)', linewidth="3")
    plt.plot(t[0::X], (dPo216[0::X]/dt),'g',label='Po216 (alpha 6.9 MeV)')
    plt.plot(t[0::X], dPb212[0::X]/dt,'k',label='Pb212 (beta Q=569.1 keV)', linewidth=4)
    plt.plot(t[0::X], dBi212[0::X]/dt,'c',label='Bi212 (beta[64%] Q=2.3 MeV; alpha[36%] 6.2 MeV)')
    plt.plot(t[0::X], (dPo212[0::X]/dt),'b',label='Po212 (alpha 9.0 MeV)')
    plt.plot(t[0::X], dTl208[0::X]/dt,'y',label='Tl208 (beta Q=5.0 MeV)')

#plt.plot([0, 0], [0, 100], color='g', label="Start of Injection")
plt.plot([inj_L/24.0, inj_L/24.0], [0, 100], color='r', label="End of Injection")
plt.plot([0, time/24.0], [BG, BG], color = "firebrick", linestyle='dashed', label="(Pb212 Equivalent) ER Background Rate") #  1.5-6.5 keV for Pb212
    
if Cryopumping==True:
    plt.plot([1, 1], [0, 100], color='g')
    plt.plot([(24.0+inj_L)/24.0, (24.0+inj_L)/24.0], [0, 100], color='r')

plt.legend(loc="upper right")
plt.ylabel("Activities (Bq)")
plt.xlabel("Time (days)")
plt.yscale("log")
plt.ylim(0.001, int(max(Total)*1.5))
plt.grid()
plt.xlim(0,math.ceil(max(t)))
if Continuous==True:
    plt.title("Rn220 Chain Activities ("+str(A_0)+" Hz Pb212 Accumulation)")
else:
    plt.title("Rn220 Chain Activities ("+str(A_0)+" Bq Rn220 Injection)")

#----------- Cumulative Events Plot -------------#
plt.figure(2)
if(Pb212_only_mode==False):
    plt.plot(t[0::X], Rn220C[0::X],'r',label='Rn220 (alpha 6.4 MeV)', linewidth="3")
    plt.plot(t[0::X], Po216C[0::X],'g',label='Po216 (alpha 6.9 MeV)', linewidth="2")
plt.plot(t[0::X], Pb212C[0::X],'k',label='Pb212 (beta Q=569.1 keV)', linewidth=4)
plt.plot(t[0::X], Bi212C[0::X],'c',label='Bi212 (beta[64%] Q=2.3 MeV; alpha[36%] 6.2 MeV)')
plt.plot(t[0::X], Po212C[0::X],'b',label='Po212 (alpha 9.0 MeV)')
plt.plot(t[0::X], Tl208C[0::X],'y',label='Tl208 (beta Q=5.0 MeV)')
plt.plot([0, time/24.0], [BG_Counts, BG_Counts], color = "firebrick", linestyle='dashed')
plt.plot([0, time/24.0], [10*BG_Counts, 10*BG_Counts], color = "firebrick", linestyle='dashed')
plt.plot([0, 0], [0, 12000000], color='g', label="Start of Injection")
plt.plot([inj_L/24.0, inj_L/24.0], [0, 12000000], color='r', label="End of Injection")
if Cryopumping==True:
    plt.plot([1, 1], [0, 12000000], color='g')
    plt.plot([(24.0+inj_L)/24.0, (24.0+inj_L)/24.0], [0, 12000000], color='r')
plt.legend(loc="lower right")
plt.ylabel("Number of Events")
plt.xlabel("Time (days)")
plt.ylim(0, 1000000)
plt.grid()
plt.xlim(0,math.ceil(max(t)))
if Continuous==True:
    plt.title("Rn220 Chain Cumulative Events ("+str(A_0)+" Hz Pb212 Accumulation)")
else:
    plt.title("Rn220 Chain Cumulative Events ("+str(A_0)+" Bq Rn220 Injection)")

#------------ Number of Atoms Plot --------------#
plt.figure(3)
if(Pb212_only_mode==False):
    plt.plot(t[0::X],Rn220[0::X],'r',label='Rn220')
    plt.plot(t[1::X],Po216[1::X],'g',label='Po216')
plt.plot(t[0::X],Pb212TotalA[0::T_X],'grey',label='Total')
plt.plot(t[0::X],Pb212[0::X],'k',label='Pb212')
plt.plot(t[0::X],Bi212[0::X],'c',label='Bi212')
plt.plot(t[0::X],Po212[0::X],'b',label='Po212')
plt.plot(t[0::X],Tl208[0::X],'y',label='Tl208')
plt.plot(t[0::X],Pb208[0::X],'m',label='Pb208')
plt.plot([0, 0], [0, 10000000000], color='g', label="Start of Injection")
plt.plot([inj_L/24.0, inj_L/24.0], [0, 10000000000], color='r', label="End of Injection")
if Cryopumping==True:
    plt.plot([(2*inj_L)/24.0, (2*inj_L)/24.0], [0, 100000000000], color='g')
    plt.plot([(3*inj_L)/24.0, (3*inj_L)/24.0], [0, 100000000000], color='r')
plt.legend(loc='lower right')
plt.ylabel('Number of Atoms')
plt.xlabel('Time (days)')
plt.yscale('log')
plt.grid()
plt.ylim(0.1, max([max(Rn220),max(Po216),max(Bi212),max(Po212),max(Pb212),max(Tl208),max(Pb208)])*1.5)
plt.xlim(0,math.ceil(max(t)))
if Continuous==True:
    plt.title("# of Rn220 chain atoms ("+str(A_0)+" Hz Pb212 Accumulation)")
else:
    plt.title("# of Rn220 chain atoms ("+str(A_0)+" Bq Rn220 Injection)")

plt.show()