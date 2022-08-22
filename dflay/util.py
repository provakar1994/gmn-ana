# helpful utility functions
import os 
import sys
import json
import csv
import shutil
import datetime

#_______________________________________________________________________________
def utcNow(): 
   theTime = datetime.datetime.utcnow() 
   return int(theTime.strftime("%s")) 
#_______________________________________________________________________________
def deleteDir(path): 
   shutil.rmtree(path) 
#_______________________________________________________________________________
def createDir(path): 
   if not os.path.exists(path): 
      os.makedirs(path)
