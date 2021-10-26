#Input starting Rb83 activity, starting date, injection date, and minutes since generator pumpout.
#Script outputs injection-date-Rb-activity, and Kr83m activity in generator.

#Rb83 halflife  = 86.2 days  = 7,447,680 seconds
#Kr83m halflife = 1.83 hours = 6588 seconds

import math
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np

# INPUTS ########################################################################################

#A  = 17172                       # starting Rb83 activity (Bq) on day generator was dosed - Generator that we dosed @ Yale
#d0 = datetime.date(2018, 8, 16)  # date that generator was dosed with Rb83 or calibration date

#A  = 29700                       # starting Rb83 activity (Bq) on day generator was dosed - Generator that Sid Shipped to us
#d0 = dt.date(2020, 11, 23)       # date that generator was dosed with Rb83 or calibration date

#A  = 37000000                    # Starting activity of UMass Rb83 Bottle on February 21, 2021 
#A  = 74000                       # Starting specific activity of UMass Rb83 Bottle on February 21, 2021 (Bq/uL)
#d0 = dt.date(2021, 2, 21)        # Calibration date of Rb83 Bottle

A  = 45500                       # Starting Rb83 activity (Bq) on day generator was dosed (First OVC charcoal generator)
d0 = dt.date(2021, 7, 14)        # Date of generator dosing 

#A  = 34200                       # Starting Rb83 activity (Bq) on day generator was dosed (Minimal charcoal [OVC] generator)
#d0 = dt.date(2021, 8, 10)  # Date of generator dosing

#---------------------------------
A  = 100000
d1 = dt.date(2021, 7, 21)  # injection date or date that we want to know the activity

delta = d1 - d0 # = time between the dates
time  = 720     # time since generator pumpout to plot [MINUTES] (Only for Figures 1 and 2)

#################################################################################################

Rb_dc  = 0.000000093       # Rb83 decay constant (ln(2)/halflife) [seconds]
Kr_dc  = 0.000105214       # Kr83m decay constant (ln(2)/halflife) [seconds]
Rb_dc_days = 0.00804115    # Rb83 decay constant [days]
BF     = 0.779             # Branching fraction for Kr83m
A_Ci   = A*.00000000002703 # convert activity in Bq to Ci
A_f    = round(A*math.exp(-Rb_dc*((delta.days*86400.0)+delta.seconds)), 12) # Activity on end date
A_f_Ci = round(A_Ci*math.exp(-Rb_dc*((delta.days*86400.0)+delta.seconds)), 12) #Activity on end date

d1_str = d1.strftime("%b-%d-%Y")
d0_str = d0.strftime("%b-%d-%Y")

print("Initial Rb Activity = "+str(round(A,2))+" Bq = "+str(round(A_Ci*1000000.0,4))+" uCi")
print ("Rb Rate "+str(delta.days)+" days later = "+str(A_f)+" Bq = "+str(round(A_f_Ci*1000000.0,4))+" uCi")
print ("Number of Kr83m atoms flushed out of generator = "+str(round(A_f*(1/Kr_dc))))

#Activity = 1675 # Bq
#print("Number of Kr83m atoms in "+str(Activity)+" Bq = "+str(round(Activity*(1/Kr_dc))))
#print("Mass of Kr83m atoms in "+str(Activity)+" Bq = "+str(((Activity*(1/Kr_dc))/(6.0221409e+23))*83.8)+" g.")

t_convert  = int(time*60.0)      # convert the input time (minutes) to seconds
t          = np.zeros(t_convert)
Rb_Rate    = np.zeros(t_convert)
Kr_Rate    = np.zeros(t_convert)
t[0]       = 0
Rb_Rate[0] = A_f
Kr_Rate[0] = 0

for i in range(1,t_convert):
    t[i]       = i
    Rb_Rate[i] = (A_f*math.exp(-Rb_dc*t[i]))
    Kr_Rate[i] = (((A_f)/(Kr_dc-Rb_dc))*(math.exp(-Rb_dc*t[i])-math.exp(-Kr_dc*t[i])))*Kr_dc*BF

t_days       = np.zeros(delta.days)
Rb_rate_days = np.zeros(delta.days)
t_days[0]       = 0
Rb_rate_days[0] = A

for i in range(1,delta.days):
    t_days[i]       = i
    Rb_rate_days[i] = A*math.exp(-Rb_dc_days*t[i])


# Plot of Rb83 and Kr83m activities in the generator after pump-out has stopped
plt.figure(1)
plt.plot(t/(60.0),Rb_Rate,'k',label='Rb83')
plt.plot(t/(60.0),Kr_Rate,'b',label='Kr83m')
plt.xlabel('Minutes Since Generator Pumpout'); plt.ylabel('Rate (Bq)')
plt.title("Kr83m Generator Activities After Pumpout"); plt.legend(loc='lower right');
#plt.yscale('log'); plt.ylim(1000,10000)
plt.grid()

# Plot of Kr83m activity as a % of equillibrium activity after pump-out has stopped
plt.figure(2)
plt.plot(t/(60.0),(Kr_Rate/(Rb_Rate*BF))*100.0,'r')
plt.xlabel('Minutes Since Generator Pumpout'); plt.ylabel('% Equillibrium Activity')
plt.title("Kr83m Generator Activity After Pumpout"); plt.ylim(0,100); plt.grid()
 
y1 = Rb_rate_days*0.0001
y2 = Rb_rate_days*0.1

# Plot of activity of the Rb83 solution 
plt.figure(3)
plt.plot(t_days,Rb_rate_days,'k',label='Rb83 Source Activity')
plt.plot(t_days,Rb_rate_days*0.1,'r',label='Max Inj Kr83m Activity')
plt.plot(t_days,Rb_rate_days*0.0001,'b',label='Min Inj Kr83m Activity')
plt.fill_between(t_days, y1, y2, where=y2>=y1, facecolor='g', interpolate=True, alpha=0.5)
#plt.plot([delta.days, delta.days],[100, A_f], color='red', linestyle='dashed')
#plt.plot([0, delta.days],[A_f, A_f], color='red', linestyle='dashed')
plt.yscale('log'); plt.legend(loc='upper right'); plt.xlabel('Days Since '+d0_str); plt.ylabel('Activity (Bq)')
plt.xlim(0,t_days[-1]+(0.05*t_days[-1])); #plt.ylim(min(Rb_Rate)-(0.4*min(Rb_Rate)), 2*A)
plt.title("Rb83 & Inj Kr83m Activities from "+d0_str+" to "+d1_str); plt.grid();

