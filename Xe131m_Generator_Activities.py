##############################################
### Calculating I131 and Xe131m Activities ###
##############################################

### INFO ###########################################

#I131 halflife  = 8.02 days = 692,928 seconds
#Xe131m halflife =  11.93 days = 1,030,752 seconds

####################################################

import math
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np

# INPUTS ########################################################################################

#A    = 200540                               # Calibrated I131 activity (Bq) of first I131 capsule
#d0   = dt.datetime(2020, 1, 15, 0, 0, 0)    # Calibration Date of first I131 capsule
#A    = 185000                               # Calibrated I131 activity (Bq) of second I131 capsule
#d0   = dt.datetime(2020, 3, 5, 0, 0, 0)     # Calibration Date of second I131 capsule
#A    = 393680                               # Calibrated I131 activity (Bq) of third I131 capsule
#d0   = dt.datetime(2020, 11, 2, 0, 0, 0)    # Calibration Date of third I131 capsule
#A    = 2006300                              # Calibrated I131 activity (Bq) of fourth I131 capsule
#d0   = dt.datetime(2021, 1, 6, 15, 0, 0)    # Calibration date of the fourth I131 capsule
#A    = 189070                               # Calibrated I131 activity (Bq) of fifth I131 capsule
#d0   = dt.datetime(2021, 8, 3, 14, 0, 0)    # Calibration date of the fifth I131 capsule

A    = 301920                               # Calibrated I131 activity (Bq) of fifth I131 capsule
d0   = dt.datetime(2021, 9, 2, 13, 0, 0)    # Calibration date of the fifth I131 capsule

d1   = dt.datetime(2021, 9, 21, 12, 0, 0)  # Date on which to start Wait Time

stopwatch = 60*24*60      # Minutes since last pump-out of generator stopped

#################################################################################################

BF    = 0.0083      # Branching fraction for Xe131m from I131 (0.0083 is from a best fit to previous injections)

I_dc  = 0.000001    # Decay Constant (ln(2)/halflife) [sec^-1]
Xe_dc = 0.00000067  # Decay Constant (ln(2)/halflife) [sec^-1]

d0_str = d0.strftime("%b-%d-%Y"); d1_str = d1.strftime("%b-%d-%Y")

#Rate Outputs
delta = d1 - d0;
A_I = A*math.exp(-I_dc*((delta.days*86400.0)+delta.seconds)) # Activity at start of wait time
A_f = A*math.exp(-I_dc*((delta.days*86400.0)+delta.seconds+(stopwatch*60.0)))                  # Activity at Injection time
print("\nI131 Rate at calibration           = "+str(round(A))+" Bq = "+str(round(A*.00002703,4))+" uCi - ("+d0_str+")")
print("I131 Rate in Generator Start of Wait Time = "+str(round(A_I))+" Bq = "+str(round(A_I*.00002703,4))+" uCi - ("+d1_str+") ("+str(delta.days)+" days later)")

#Calculating Activities 
time = int(delta.days*24.0*60.0*1.2) # Defining the amount of time to calculate I131 activities for in minutes
wait = int((stopwatch*60.0)+96)      # Defining the amount of time to calculate Xe131m activities after end of pumpout, in seconds

t = np.zeros(time); I_Rate = np.zeros(time); w_t = np.zeros(wait*2); Xe_Rate = np.zeros(wait*2)
t[0] = 0.0;         I_Rate[0] = A;           w_t[0] = 0.0;           Xe_Rate[0] = 0.0

for i in range(1,time):
    t[i]       = i
    I_Rate[i]  = A*math.exp(-I_dc*t[i]*60.0)
 
for i in range(1,wait*2):
    w_t[i]     = i    
    Xe_Rate[i] = (A_I/(Xe_dc-I_dc))*(math.exp(-I_dc*w_t[i])-math.exp(-Xe_dc*w_t[i]))*Xe_dc*BF

print("\nExpected Injected Activity from "+str(stopwatch)+" min Wait Time = "+str(round(Xe_Rate[wait-1], 3))+" Bq")

### Figures ################################################################################################
############################################################################################################

plt.figure(1)
plt.plot(t/(60.0*24.0),I_Rate,'k')
plt.plot([delta.days, delta.days], [1, A], 'r')
plt.xlabel('Days Since Calibration Date, '+d0_str); plt.ylabel('Rate (Bq)')
plt.xlim([0, delta.days*1.2]);
plt.title("I-131 Activity"); plt.grid(); plt.yscale('log')  

plt.figure(2)
plt.plot(w_t/60.0,Xe_Rate,'b', label="Xe131m")
plt.plot([stopwatch+1.6, stopwatch+1.6],[0, Xe_Rate[wait-1]*1.2], color='r', label="Stopwatch Time")
plt.xlim([0.1,(stopwatch+1.6)*1.1]); plt.ylim([0, Xe_Rate[wait-1]*1.2])
plt.xlabel('Time Since Generator Pump Out ('+d1_str+') (minutes)'); plt.ylabel('Xe131m Activity (Bq)')
plt.title("Generator Xe-131m Activity After Pump Out of "+str(round(A_I))+" Bq I-131")
plt.legend(); plt.grid()
plt.show()