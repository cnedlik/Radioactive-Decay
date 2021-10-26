#- Simulating the Rn222 decay chain approximately using discrete time steps -#
#- Developed by Chris Nedlik -----------------------------------------------------------------#

#                                      (99.979%) Po-214 
#                                     /               \
# Rn-222 -> Po-218 -> Pb-214 -> Bi-214                 Pb-210 -> Bi-210 -> Po-210 -> Pb-206
#                                     \               / 
#                                      (0.021%) Tl-210 

import numpy as np
import matplotlib.pyplot as plt
import math
import time as stopwatch
t0 = stopwatch.time()
#------------------------------Settings--------------------------------#

A_0   = 0.000025   # Rn222 Emanation rate (Bq)
time  = 568        # Time to simulate - hours - 
dt    = 0.2        # Time step - seconds. Small for better approximation

Purification    = False  # If true, isotopes are removed according to LZs purification time, tau = 2.35 days 
Cryopumping     = True  # The bottle will be emptied in about 16 hours. It takes 8 to refill. 

pur_dt = 10.0   # Purification time step in seconds (if Purification = True)

#-----------------------------Constants--------------------------------#

#-- Decay Constants --#
Rn222_C = math.log(2)/330177.6            # Halflife = 3.823 days
Po218_C = math.log(2)/182.4               # Halflife = 3.04 minutes
Pb214_C = math.log(2)/1584.0              # Halflife = 26.4 minutes
Bi214_C = math.log(2)/1194.0              # Halflife = 19.9 minutes
Po214_C = math.log(2)/0.1643              # Halflife = 0.1643 seconds
Tl210_C = math.log(2)/78.0                # Halflife = 1.30 minutes
Pb210_C = math.log(2)/(22.3*365.24*24.0*60.0*60.0) # Halflife = 22.3 years
Bi210_C = math.log(2)/433123.2            # Halflife = 5.013 days

alpha_ion_frac = 0.503
beta_ion_frac = 0.764

Bi_Po   = 0.99979 # Bi212->Po212 Branching Fraction
Bi_Tl   = 0.00021 # Bi212->Tl208 Branching Fraction

#----------------------------------------------------------------------#

secs  = int(60.0*60.0*time)             # Run time in seconds
steps = int(secs/dt)                    # Number of dt steps
pur_frac = pur_dt/(2.35*24.0*60.0*60.0) # fraction of isotopes removed each pur_dt (tau = 2.35 days)

t     = np.zeros(steps)      # time array

#-- Number of Atoms --#  #--- Decayed Atoms ---#   #-- Cumulative Events --#
Rn222 = np.zeros(steps); dRn222 = np.zeros(steps); Rn222C = np.zeros(steps);
Po218 = np.zeros(steps); dPo218 = np.zeros(steps); Po218C = np.zeros(steps);
Pb214 = np.zeros(steps); dPb214 = np.zeros(steps); Pb214C = np.zeros(steps);
Bi214 = np.zeros(steps); dBi214 = np.zeros(steps); Bi214C = np.zeros(steps);
Po214 = np.zeros(steps); dPo214 = np.zeros(steps); Po214C = np.zeros(steps);
Tl210 = np.zeros(steps); dTl210 = np.zeros(steps); Tl210C = np.zeros(steps);
Pb210 = np.zeros(steps); dPb210 = np.zeros(steps); Pb210C = np.zeros(steps);
Bi210 = np.zeros(steps); dBi210 = np.zeros(steps); Bi210C = np.zeros(steps);
Po210 = np.zeros(steps);

Rn222[0]  = A_0/Rn222_C  # Number of atoms 
dRn222[0] = A_0*dt       # Initial dN/dt

for i in range(1,int(secs/dt)):   
    t[i]   = dt*i
    
    #---------------- Rn222 -----------------#
    dRn222_tmp = (Rn222_C*Rn222[i-1]*dt)  # Number of decaying atoms
    if(dRn222_tmp > Rn222[i-1]):          # If the # decaying > # remaining
        Rn222[i]  = 0                     # Then the # of atoms remaining is 0 (not negative)
        dRn222[i] = Rn222[i-1]            # The # of decayed atoms = the # that were left
    else:
        Rn222[i]  = Rn222[i-1]-dRn222_tmp
        dRn222[i] = dRn222_tmp
    Rn222C[i] = Rn222C[i-1]+dRn222[i]

    if Cryopumping==True:

        full   = ((t[i]<3600.0*16.0)) | \
        ((t[i]>3600.0*24.0)  & (t[i]<3600.0*40.0))  | \
        ((t[i]>3600.0*48.0)  & (t[i]<3600.0*64.0))  | \
        ((t[i]>3600.0*72.0)  & (t[i]<3600.0*88.0))  | \
        ((t[i]>3600.0*96.0)  & (t[i]<3600.0*112.0)) | \
        ((t[i]>3600.0*120.0) & (t[i]<3600.0*136.0)) | \
        ((t[i]>3600.0*144.0) & (t[i]<3600.0*160.0)) | \
        ((t[i]>3600.0*168.0) & (t[i]<3600.0*184.0)) | \
        ((t[i]>3600.0*192.0) & (t[i]<3600.0*208.0)) | \
        ((t[i]>3600.0*216.0) & (t[i]<3600.0*232.0))
        
        if(full == True): 
            Rn222[i] = A_0/Rn222_C #Rn222[i-1]-dRn222_tmp+(A_0*dt)            
            
    #---------------- Po218 ----------------#
    dPo218_tmp = (Po218_C*Po218[i-1]*dt)
    if(dPo218_tmp > Po218[i-1]):
        Po218[i]  = 0
        dPo218[i] = Po218[i-1]
    else:
        Po218[i]  = Po218[i-1]+dRn222[i]-dPo218_tmp
        dPo218[i] = dPo218_tmp
    Po218C[i] = Po218C[i-1]+dPo218[i]
    
    #---------------- Pb214 ----------------#
    dPb214_tmp = (Pb214_C*Pb214[i-1]*dt)
    if(dPb214_tmp > Pb214[i-1]):
        Pb214[i]  = 0
        dPb214[i] = Pb214[i-1]
    else:
        Pb214[i]  = Pb214[i-1]+dPo218[i]-dPb214_tmp
        dPb214[i] = dPb214_tmp
    Pb214C[i] = Pb214C[i-1]+dPb214[i]
        
    #---------------- Bi214 ----------------#
    dBi214_tmp = (Bi214_C*Bi214[i-1]*dt)
    if(dBi214_tmp > Bi214[i-1]):
        Bi214[i]  = 0
        dBi214[i] = Bi214[i-1]
    else:
        Bi214[i] = Bi214[i-1]+dPb214[i]-dBi214_tmp
        dBi214[i] = dBi214_tmp
    Bi214C[i] = Bi214C[i-1]+dBi214[i]

    #---------------- Po214 ----------------#
    dPo214_tmp = (Po214_C*Po214[i-1]*dt)
    if(dPo214_tmp > Po214[i-1]):
        Po214[i]  = 0
        dPo214[i] = Po214[i-1]
    else:
        Po214[i]  = Po214[i-1]+(dBi214[i]*Bi_Po)-dPo214_tmp
        dPo214[i] = dPo214_tmp
    Po214C[i] = Po214C[i-1]+dPo214[i]

    #---------------- Tl210 ----------------#
    dTl210_tmp = (Tl210_C*Tl210[i-1]*dt)
    if(dTl210_tmp > Tl210[i-1]):
        Tl210[i]  = 0
        dTl210[i] = Tl210[i-1]
    else:
        Tl210[i]  = Tl210[i-1]+(dBi214[i]*Bi_Tl)-dTl210_tmp
        dTl210[i] = dTl210_tmp
    Tl210C[i] = Tl210C[i-1]+dTl210[i]

    #---------------- Pb210 ----------------#
    dPb210_tmp = (Pb210_C*Pb210[i-1]*dt)
    if(dPb210_tmp > Pb210[i-1]):
        Pb210[i]  = 0
        dPb210[i] = Pb210[i-1]
    else:
        Pb210[i]  = Pb210[i-1]+dPo214[i]+dTl210[i]-dPb210_tmp
        dPb210[i] = dPb210_tmp
    Pb210C[i] = Pb210C[i-1]+dPb210[i]

    #---------------- Bi210 ----------------#
    dBi210_tmp = (Bi210_C*Bi210[i-1]*dt)
    if(dBi210_tmp > Bi210[i-1]):
        Bi210[i]  = 0
        dBi210[i] = Bi210[i-1]
    else:
        Bi210[i]  = Bi210[i-1]+dPb210[i]-dBi210_tmp
        dBi210[i] = dBi210_tmp
    Bi210C[i] = Bi210C[i-1]+dBi210[i]
    
    #---------------- Po210 ----------------#
    Po210[i] = Pb210[i-1]+dBi210[i]
    
    #------------- Purification ------------#

    if Purification == True:
        if(t[i]%pur_dt == 0):
            Po218[i] = Po218[i]-(Po218[i]*pur_frac)
            Pb214[i] = Pb214[i]-(Pb214[i]*pur_frac)
            Bi214[i] = Bi214[i]-(Bi214[i]*pur_frac)
            Po214[i] = Po214[i]-(Po214[i]*pur_frac)
            Tl210[i] = Tl210[i]-(Tl210[i]*pur_frac)
            Pb210[i] = Pb210[i]-(Pb210[i]*pur_frac)
            Bi210[i] = Bi210[i]-(Bi210[i]*pur_frac)
            Po210[i] = Po210[i]-(Po210[i]*pur_frac)            

# Total Activities
Total = (dRn222[0::2]/dt)+(dPo218[0::2]/dt)+(dPb214[0::2]/dt)+(dBi214[0::2]/dt)+(dPo214[1::2]/dt)+(dTl210[0::2]/dt)+(dPb210[0::2]/dt)+(dBi210[0::2]/dt) 
Pb214Total = (dPb214[0::2]/dt)+(dBi214[0::2]/dt)+(dPo214[1::2]/dt)+(dTl210[0::2]/dt)+(dPb210[0::2]/dt)+(dBi210[0::2]/dt)   # Total Activity w/o Rn222, Po218
BetaActivities = (dPb214[0::2]/dt)+(dBi214[0::2]/dt)+(dTl210[0::2]/dt)+(dPb210[0::2]/dt)+(dBi210[0::2]/dt) # Total activity of beta decayers

# Total Cumulative Events
TotalC = Rn222C[0::2]+Po218C[0::2]+Pb214C[0::2]+Bi214C[0::2]+Po214C[1::2]+Tl210C[0::2]+Pb210C[0::2]+Bi210[0::2]

BA = (dPb214[-1]/dt)+(dBi214[-1]/dt)+(dTl210[-1]/dt)+(dPb210[-1]/dt)+(dBi210[-1]/dt)

print("Computation Time = "+str(round(stopwatch.time()-t0,1))+" sec")
print("Final Total Activity = "+str(Total[-1]/dt)+" Bq")
print("Final 222Rn Activity = "+str(dRn222[-1]/dt)+" Bq")
print("Final Beta Activity = "+str(BA)+" Bq")


#---------------- Activity Plot -----------------#

X = 20  # must be even

T_X = int(X/2.0)
t=t/(60.0*60.0*24.0)

plt.figure(1)
plt.plot(t[0::X], Total[0::T_X], 'grey', label='Total')
plt.plot(t[0::X], dRn222[0::X]/dt,'k',label='Rn222')
plt.plot(t[0::X], (dPo218[0::X]/dt),'b',label='Po218')
plt.plot(t[0::X], dPb214[0::X]/dt,'g',label='Pb214')
plt.plot(t[0::X], dBi214[0::X]/dt,'r',label='Bi214')
plt.plot(t[1::X], (dPo214[1::X]/dt),'c',label='Po214')
plt.plot(t[0::X], dTl210[0::X]/dt,'y',label='Tl210')
plt.plot(t[0::X], dPb210[0::X]/dt,'olive',label='Pb210')
plt.plot(t[0::X], dBi210[0::X]/dt,'dodgerblue',label='Bi210')
plt.legend(loc="lower right")
plt.ylabel("Activities (Bq)")
plt.xlabel("Time (days)")
plt.yscale("log")
plt.ylim(0.0000000001, int(max(Total)*1.5)+1)
plt.grid()
plt.xlim(0,math.ceil(max(t)))
plt.title("Rn222 Chain Activities with Continuous Injection ("+str(A_0)+" Bq Rn222)")

#----------- Cumulative Events Plot -------------#
plt.figure(2)

plt.plot(t[0::X], TotalC[0::T_X], 'grey', label='Total')
plt.plot(t[0::X], Rn222C[0::X],'k',label='Rn222')
plt.plot(t[0::X], Po218C[0::X],'b',label='Po218')
plt.plot(t[0::X], Pb214C[0::X],'g',label='Pb214')
plt.plot(t[0::X], Bi214C[0::X],'r',label='Bi214')
plt.plot(t[1::X], Po214C[1::X],'c',label='Po214')
plt.plot(t[0::X], Tl210C[0::X],'y',label='Tl210')
plt.plot(t[0::X], Pb210C[0::X],'olive',label='Pb210')
plt.plot(t[0::X], Bi210C[0::X],'dodgerblue',label='Bi210')
plt.legend(loc="lower right")
plt.ylabel("Number of Events")
plt.xlabel("Time (days)")
plt.yscale("log")
plt.ylim(0.0000001, int(max(TotalC)*1.5))
plt.grid()
plt.xlim(0,math.ceil(max(t)))
plt.title("Rn222 Chain Cumulative Events w/ Continuous Injection ("+str(A_0)+" Bq Rn222)")

#------------ Number of Atoms Plot --------------#

plt.figure(3)
plt.plot(t[0::X],Rn222[0::X],'k',label='Rn222')
plt.plot(t[1::X],Po218[1::X],'b',label='Po218')
plt.plot(t[0::X],Pb214[0::X],'g',label='Pb214')
plt.plot(t[0::X],Bi214[0::X],'r',label='Bi214')
plt.plot(t[0::X],Po214[0::X],'c',label='Po214')
plt.plot(t[0::X],Tl210[0::X],'y',label='Tl210')
plt.plot(t[0::X],Pb210[0::X],'olive',label='Pb210')
plt.plot(t[0::X],Bi210[0::X],'dodgerblue',label='Bi210')
plt.legend(loc='lower right')
plt.ylabel('Number of Atoms')
plt.xlabel('Time (days)')
plt.yscale('log')
plt.grid()
plt.ylim(0.1, max([max(Rn222),max(Po218),max(Bi214),max(Po214),max(Pb214),max(Tl210),max(Pb210),max(Bi210)])*1.5)
plt.xlim(0,math.ceil(max(t)))
plt.title("# of 222Rn chain atoms w/ Continuous Injection ("+str(A_0)+" Bq Rn222)")

plt.show()
