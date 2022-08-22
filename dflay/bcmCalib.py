#!/usr/bin/python
# A script to run BCM analysis

import os 
import sys
import json
import csv
import datetime

import util 

# read in the configuration file 
NARG = len(sys.argv) 
if(NARG<2): 
  print("[bcmCalib]: Invalid input!") 
  print("[bcmCalib]: Usage: ./bcmCalib.py config-name") 
  sys.exit(0) 

tag     = sys.argv[1]
homeDir = os.getcwd()  

startTime = util.utcNow() 

# determine JSON config file path  
jConfPath = os.getcwd()+ "/input/{0}/conf.json".format(tag)
# read JSON config file 
# confData = json.loads( open(jConfPath).read() )

# create output directories 
outDir  = os.getcwd() + "/output/{0}".format(tag); 
logDir  = os.getcwd() + "/output/{0}/log".format(tag); 
plotDir = os.getcwd() + "/output/{0}/plots/".format(tag); 

print( "[bcmCalib]: Creating directory '{0}'".format(outDir) ) 
util.createDir(outDir)
print( "[bcmCalib]: Creating directory '{0}'".format(logDir) ) 
util.createDir(logDir)
print( "[bcmCalib]: Creating directory '{0}'".format(plotDir) ) 
util.createDir(plotDir)

#===============================================================================
# process all cuts (on all BCM variables)  
varList    = ["unser","u1","unew","d1","d3","d10","dnew"]
scriptName = "bcmCalibProcessCuts.cxx"
for var in varList: 
  cmd = "analyzer -q -b -l '{0}(\"{1}\",\"{2}\")'".format(scriptName,jConfPath,var)
  print(cmd)
  os.system(cmd) 
#===============================================================================
# compute BCM calibration coefficients
scriptName = "bcmCalibrate.cxx"
cmd        = "analyzer -q -b -l '{0}(\"{1}\")'".format(scriptName,jConfPath)
print(cmd)
os.system(cmd) 
#===============================================================================
# compute BCM calibration coefficients for EPICS 
scriptName = "bcmCalibrateEPICS.cxx"
cmd        = "analyzer -q -b -l '{0}(\"{1}\")'".format(scriptName,jConfPath)
print(cmd)
os.system(cmd) 
#===============================================================================
# create DB files 
scriptName = "bcmCalibDBFiles.cxx"
cmd        = "analyzer -q -b -l '{0}(\"{1}\")'".format(scriptName,jConfPath)
print(cmd)
os.system(cmd) 
#===============================================================================
# cleanup output directory
os.chdir(outDir)
# make csv dir
if not os.path.exists("./csv"):  
   os.mkdir("csv",0755) 
   print("[bcmCalib]: Directory csv created")
else:  
   print("[bcmCalib]: Directory csv already exists")
# move files to csv directory  
os.system("mv *.csv csv")
os.system("mv csv/*result* .")  # retrieve the final result file 
# move back to home directory 
os.chdir(homeDir)
print("[bcmCalib]: Done.")   
#===============================================================================
# compute job stats 
endTime         = util.utcNow() 
elapsedTime     = float(endTime-startTime) 
elapsedTime_min = elapsedTime/60. 
print("[bcmCalib]: Elapsed time = {0:.3f} sec ({1:0.3f} min)".format(elapsedTime,elapsedTime_min))
#===============================================================================
