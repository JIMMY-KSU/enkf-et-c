
import numpy as np
from datetime import date
import csv
import os

savepath=os.path.dirname(os.path.abspath(__file__)) 

# flux tower ET data
fluxfile=savepath+'\\AmeriFlux_US-Goo_reduced_JJA02.csv'

#conversion factor from latent heat flux (W/m2) to ET in mm/day
rho=1000 #kg/m3 density of water
lv=2257000 #J/kg latent heat of 
secsperday=86400
letoet=rho**(-1.)*lv**(-1)*secsperday*1000

with open(fluxfile, 'rb') as csvfile:
    reading = csv.reader(csvfile, delimiter='\t')
    a=0
    data=[]
    for row in reading:
        if a==2:
            headings0=row
        elif a > 2 : 
            b=row[0].split(',')
            data.append(np.array(b).astype(np.float))
        a=a+1
    
    headings=headings0[0].split(',')    
    dm=np.array(data)
    tdm=dm.transpose()
    labels={}
    for i2 in range(len(headings)):
        foo=headings[i2]
        labels[foo]=i2
        assert labels[headings[i2]]==i2

    if 'TIMESTAMP_START' in labels:
        timestart0=tdm[labels['TIMESTAMP_START']]
    if 'TIMESTAMP_END' in labels:
        timeend0=tdm[labels['TIMESTAMP_END']]
    if 'LE' in labels:
        le0=tdm[labels['LE']]
    if 'TA' in labels:
        TA0=tdm[labels['TA']]
    if 'NETRAD' in labels:
        NETRAD0=tdm[labels['NETRAD']]
    if 'PA' in labels:
        PA0=tdm[labels['PA']]
    
wherele=np.where(le0 != -9999)[0] #only keep points with ET data. "-9999" is NoData value.
le=le0[np.array(wherele)]
TA=TA0[np.array(wherele)]
NETRAD=NETRAD0[np.array(wherele)]
PA=PA0[np.array(wherele)]
timestart=timestart0[np.array(wherele)]
timeend=timeend0[np.array(wherele)]

#split up yr-month-day-hour-minute format
yrstart=[]
mostart=[]
daystart=[]
hrstart=[]
minstart=[]
hrminstart=[]
for j in range(len(timestart)):
    yrstart.append(str(timestart[j])[:4])
    mostart.append(str(timestart[j])[4:6])
    daystart.append(str(timestart[j])[6:8])
    hrstart.append(str(timestart[j])[8:10])
    minstart.append(str(timestart[j])[10:12])
    hrminstart.append(str(timestart[j])[8:12])

yrstart=np.array(yrstart)
mostart=np.array(mostart)
daystart=np.array(daystart)
hrstart=np.array(hrstart)
minstart=np.array(minstart)
hrminstart=np.array(hrminstart)

# make an nd.array of the measurements at each timestep:
#rows: timestamp in julian day with fraction, ET (mm/day), TA (air temp in C), PA (atmospheric pressure in kPa), NETRAD (net radiation in W/m^2)
d0 = date(int(yrstart[0]), 1, 1) #start time for calculating julian day

fluxdata=np.ndarray(shape=(5,len(le)))
for i in range(len(le)):
    d1 = date(int(yrstart[0]), int(mostart[i]), int(daystart[i]))
    delta = d1 - d0
    jtime=delta.days+1+float(hrstart[i])/24+(float(minstart[i])+15)/(60*24)
    fluxdata[0,i]=jtime
    fluxdata[1,i]=le[i]*letoet
    fluxdata[2,i]=TA[i]
    fluxdata[3,i]=PA[i]
    fluxdata[4,i]=NETRAD[i]

np.savez(savepath+'\\fluxtower_et.npz', fluxdata=fluxdata) 

