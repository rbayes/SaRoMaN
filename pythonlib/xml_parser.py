#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 12/11-2015
#Created on 12/11-2015
#######################################################################################################################
#General python import
#######################################################################################################################
from xml.dom import minidom
import sys
#######################################################################################################################
#Importing own python files
#######################################################################################################################

#######################################################################################################################
#Class generation
#######################################################################################################################
class xml_parser:
    
    def __init__(self,parsingFile,outFile):
        self.parsingFile = parsingFile
        self.outFile = outFile
        self.xmldoc = minidom.parse(self.parsingFile)


    def Parse_file(self):
        #Variables used
        volumeNodeList = self.xmldoc.getElementsByTagName('volume')
        solidRefRefList = []
        volumeRefRefList = []

        boxNodeList = self.xmldoc.getElementsByTagName('box')
        boxXList =[]
        boxYList =[]
        boxZList =[]
        boxXListEval =[]
        boxYListEval =[]
        boxZListEval =[]
        boxDict = {}
        posDict = {}

        
        constantNodeList = self.xmldoc.getElementsByTagName('constant')
        constantDict= {}

        #Get constants
        for constantNode in constantNodeList:
            try:
                constantDict[str(constantNode.attributes['name'].value)] = float(constantNode.attributes['value'].value)
            except ValueError: #If read element is not well defined. Could intead eval the full dict.
                constantDict[str(constantNode.attributes['name'].value)] = eval(str(constantNode.attributes['value'].value),constantDict)

                continue

        
        #Get volume ref & position name
        for volumeNode in volumeNodeList:
            if(volumeNode.attributes['name'].value == 'Detector'):
                volumeRefList= volumeNode.getElementsByTagName("volumeref")
                for volumeRef in volumeRefList:
                    volumeRefRefList.append(volumeRef.attributes['ref'].value)
                posNameRefList= volumeNode.getElementsByTagName("position")
                for posNameRef in posNameRefList:
                    name_val = str(posNameRef.attributes['name'].value)
                    x_val = eval(str(posNameRef.attributes['x'].value),constantDict)
                    y_val = eval(str(posNameRef.attributes['y'].value),constantDict)
                    z_val = eval(str(posNameRef.attributes['z'].value),constantDict)
                    posDict[name_val] = ((x_val,y_val,z_val))
                    
        #Get solid ref
        for volumeRefRef in volumeRefRefList: 
            for volumeNode in volumeNodeList:
                if(volumeNode.attributes['name'].value == volumeRefRef):
                    solidRefList= volumeNode.getElementsByTagName("solidref")
                    for solidRef in solidRefList:
                        solidRefRefList.append(solidRef.attributes['ref'].value)

        #Get solid ref xyz
        for solidRefRef in solidRefRefList: 
            for boxNode in boxNodeList:
                if(boxNode.attributes['name'].value == solidRefRef):
                    x_val = eval(str(boxNode.attributes['x'].value),constantDict)
                    y_val = eval(str(boxNode.attributes['y'].value),constantDict)
                    z_val = eval(str(boxNode.attributes['z'].value),constantDict)
                    boxDict[str(solidRefRef)] = ((x_val,y_val,z_val))
                
        outfile = open(self.outFile,'w+')
        for key,value in boxDict.iteritems():
            outfile.write("SOLID\t%s\t%s\t%s\t%s\n"%(key[:-6],value[0],value[1],value[2]))
        for key,value in posDict.iteritems():
            outfile.write("POS\t%s\t%s\t%s\t%s\n"%(key[:-4],value[0],value[1],value[2]))
        outfile.close()

#######################################################################################################################
