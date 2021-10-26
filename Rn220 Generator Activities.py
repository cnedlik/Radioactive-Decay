import math
import datetime
import matplotlib.pyplot as plt
import numpy as np

# INPUTS ########################################################################################

A     = 20110                      # starting Th228 activity (Bq) on day generator was produced
d0    = datetime.date(2020, 5, 1)  # date that generator was dosed with Th228
d1    = datetime.date(2021, 10, 8) # injection date

#A     = 260                         # Injected Rn220 activity at UMass
#d0    = datetime.date(2020, 7, 11)  # Injection Date
#d1    = datetime.date(2021, 10, 21) # injection date

time = 7 # time since generator pumpout to plot [days]

#################################################################################################

delta = d1 - d0

#-- Decay Constants --#
Th228_C = math.log(2)/(698.21*24.0*60.0*60.0)# Halflife = 1.911 years
Ra224_C = math.log(2)/(3.632*24.0*60.0*60.0) # Halflife = 3.632 days
Rn220_C = math.log(2)/55.6                   # Halflife =  55.6 seconds

A_f = round(A*math.exp(-Th228_C*delta.days*24.0*60.0*60.0), 2)

print("Initial Th228 Activity = "+str(round(A,2))+" Bq")
print("Th288 Rate "+str(delta.days)+" days later = "+str(A_f)+" Bq")

print("Number of Rn220 atoms in "+str(52)+" Bq = "+str(round(52.0/Rn220_C)))



t_convert  = int(time*24.0*60.0*60.0)
t          = np.zeros(t_convert)
Th_Rate    = np.zeros(t_convert)
Ra_Rate    = np.zeros(t_convert)
Rn_Rate    = np.zeros(t_convert)

  #-- Number of Atoms --#       #--- Decayed Atoms ---#     #-- Cumulative Events --#
Th228 = np.zeros(t_convert); dTh228 = np.zeros(t_convert); Th228C = np.zeros(t_convert);
Ra224 = np.zeros(t_convert); dRa224 = np.zeros(t_convert); Ra224C = np.zeros(t_convert);
Rn220 = np.zeros(t_convert); dRn220 = np.zeros(t_convert); Rn220C = np.zeros(t_convert);

t[0]       = 0
Th228[0]  = A_f/Th228_C # Initial number of atoms 
dTh228[0] = A_f         # Initial dN/dt
Ra_Rate[0] = 0
Rn_Rate[0] = 0

for i in range(1,t_convert):
    t[i]       = i
    Th_Rate[i] = A_f*math.exp(-Th228_C*t[i])
    Ra_Rate[i] = Th_Rate[i]-Ra_Rate[i-1]
    Rn_Rate[i] = Ra_Rate[i]-Rn_Rate[i-1]
    
 #---------------- Th228 ----------------#
    dTh228_tmp = (Th228_C*Th228[i-1])
    if(dTh228_tmp > Th228[i-1]):
        Th228[i]  = 0
        dTh228[i] = Th228[i-1]
    else:
        Th228[i]  = Th228[i-1]-dTh228_tmp
        dTh228[i] = dTh228_tmp
    Th228C[i] = Th228C[i-1]+dTh228[i]

 #---------------- Ra224 ----------------#
    dRa224_tmp = (Ra224_C*Ra224[i-1])
    if(dRa224_tmp > Ra224[i-1]):
        Ra224[i]  = 0
        dRa224[i] = Ra224[i-1]
    else:
        Ra224[i]  = Ra224[i-1]+dTh228[i]-dRa224_tmp
        dRa224[i] = dRa224_tmp
    Ra224C[i] = Ra224C[i-1]+dRa224[i]
    
    #---------------- Rn220 ----------------#
    dRn220_tmp = (Rn220_C*Rn220[i-1])
    if(dRn220_tmp > Rn220[i-1]):
        Rn220[i]  = 0
        dRn220[i] = Rn220[i-1]
    else:
        Rn220[i]  = Rn220[i-1]+dRa224[i]-dRn220_tmp
        dRn220[i] = dRn220_tmp
    Rn220C[i] = Rn220C[i-1]+dRn220[i]

plt.figure(1)
plt.plot(t/60.0, dTh228, 'k', label='Th228')
plt.plot(t/60.0, dRa224, 'b', label='Ra224', linewidth="3")
plt.plot(t/60.0, dRn220, 'r', label='Rn220')
plt.yscale('log')
plt.legend(loc='lower right')
plt.ylabel('Rate (Bq)')
plt.xlabel('Minutes Since Generator Pumpout')
plt.title("Rn220 Generator Activities After Pumpout")
plt.ylim(0.1,100000)
plt.grid()