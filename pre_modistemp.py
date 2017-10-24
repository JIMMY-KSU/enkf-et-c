

#Import modules
import numpy as np
import arcpy
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

savepaths='E:\\for_KF\\mmts\\' 

thepoint=savepaths+'goodwin.shp'

# MODIS temperature: 8-day averages
#   Metadata: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod11a2_v006
#   Units for the below: T_day and T_night are land surface temp in kelvin. Need to multiply raw data by scaling factor of 0.02. Average over the following 8 days.
#    Time_day and night are times of observation, ranging from 0-240 for 0-24 hr
#    QC_T_day and night are quality control, must convert integer to binary '{0:08b}'.format(thenumber) , and then follow chart:
#      bit # [-2:] 00: LST error <=1K, 01: <=2K, 10: <=3K, 11: >3K
datafolder='E:\\for_KF\\MODIS_temp\\' 
datanames=["T_day", "QC_T_day", "Time_day","T_night", "QC_T_night","Time_night"]
dataindices=["0", "1","2", "4","5","6"]


# MODIS PET and ET data: 8 day totals
#Metadata: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod11a2_v006
#Units for the below: ET and PET data are total over the 8 days following date. Need to multiply raw data by scaling factor of 0.1 to get mm/8day.
#  QC is quality control 
datafolder2='E:\\for_KF\\MODIS_PET\\'
datanames2=["ET", "PET", "QC"]
dataindices2=["0","2", "4"]

tile_name=['h10v05'] 

lat=34.2547
lon=-89.8735
#below - location of the above coordinates in the map of this tile. 
locsr="-8260048.662  3808960.264"

days=np.arange(145,249,8)
days=np.delete(days,np.where(days==177)[0]) # a scene was missing from the MODIS data

#extract subdatasets
arcpy.CreateFileGDB_management(datafolder, "extracts.gdb")
arcpy.CreateFileGDB_management(datafolder2, "extracts.gdb")

#Temperature
arcpy.env.workspace=datafolder
alltiles=arcpy.ListRasters('*.hdf') 
         
for j in range(len(alltiles)):
    for i in range(len(dataindices)):
        theday=alltiles[j][13:16]
        out_name=datafolder+'extracts.gdb\\'+datanames[i]+"_"+theday
        arcpy.ExtractSubDataset_management(datafolder+'\\'+alltiles[j], out_name, dataindices[i])

#ET and PET
arcpy.env.workspace=datafolder2
alltiles=arcpy.ListRasters('*.hdf')
         
for j in range(len(alltiles)):
    for i in range(len(dataindices2)):
        theday=alltiles[j][13:16]
        out_name=datafolder2+'extracts.gdb\\'+datanames2[i]+"_"+theday
        arcpy.ExtractSubDataset_management(datafolder2+'\\'+alltiles[j], out_name, dataindices2[i])

tempdata=np.ndarray(shape=(len(dataindices),len(days)))
for i in range(len(dataindices)):
    for j in range(len(days)):
        thisras=datafolder+'extracts.gdb\\'+datanames[i]+"_"+str(days[j])
        result = arcpy.GetCellValue_management(thisras, locsr)[0] #, "2;3")
        if result == 'NoData':
            tempdata[i,j]=np.nan
        else:
            tempdata[i,j]=np.float(result)
            if i in [0,3]:
                tempdata[i,j]=tempdata[i,j]*0.02-273.15 #to degrees C
        
etdata=np.ndarray(shape=(len(dataindices2),len(days)))
for i in range(len(dataindices2)):
    for j in range(len(days)):
        thisras=datafolder2+'extracts.gdb\\'+datanames2[i]+"_"+str(days[j])
        result = arcpy.GetCellValue_management(thisras, locsr)[0] #, "2;3")
        if result == 'NoData':
            etdata[i,j]=np.nan
        else:
            etdata[i,j]=np.float(result)
            if i in [0,1]:
                etdata[i,j]=etdata[i,j]*0.1

alltemp=[]
allttimes=[]
allqc=[]   
errorcodes=['00','01','10','11']
errorcodes=np.array(errorcodes)
errorbounds=[1,2,3,6]   #in average Kelvin error bound. '6' stands in for >=3
for i in range(len(days)):
    allttimes.append(days[i]+4+tempdata[2][i]/240.)
    allttimes.append(days[i]+4+tempdata[5][i]/240.)
    alltemp.append(tempdata[0][i])
    alltemp.append(tempdata[3][i])
    if not np.isnan(tempdata[1][i]):
        qc1='{0:08b}'.format(int(tempdata[1][i]))
        allqc.append(errorbounds[np.where(errorcodes==qc1[-2:])[0]])
    else:
        allqc.append('NaN')
    if not np.isnan(tempdata[4][i]):       
        qc2='{0:08b}'.format(int(tempdata[4][i]))
        allqc.append(errorbounds[np.where(errorcodes==qc2[-2:])[0]])        
    else:
        allqc.append('NaN')

alltemp=np.array(alltemp)
allttimes=np.array(allttimes)
allqc=np.array(allqc)

np.savez(savepaths+'MODIS_temp_et2.npz',etdata=etdata,tempdata=tempdata,tempdatanames=datanames,etdatanames=datanames2,days=days, alltemp=alltemp,allttimes=allttimes,allqc=allqc)
