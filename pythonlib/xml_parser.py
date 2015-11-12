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
        boxDictEval = {}
        
        constantNodeList = self.xmldoc.getElementsByTagName('constant')
        constantDict= {}
        
        #Get volume ref
        for volumeNode in volumeNodeList:
            if(volumeNode.attributes['name'].value == 'Detector'):
                volumeRefList= volumeNode.getElementsByTagName("volumeref")
                for volumeRef in volumeRefList:
                    volumeRefRefList.append(volumeRef.attributes['ref'].value)
                    
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
                    x_val = str(boxNode.attributes['x'].value)
                    y_val = str(boxNode.attributes['y'].value)
                    z_val = str(boxNode.attributes['z'].value)
                    boxXList.append(str(x_val))
                    boxYList.append(str(y_val))
                    boxZList.append(str(z_val))
                    boxDict[str(solidRefRef)] = ((x_val,y_val,z_val))

        #Get constants
        for constantNode in constantNodeList:
            try:
                constantDict[str(constantNode.attributes['name'].value)] = float(constantNode.attributes['value'].value)
            except ValueError: #If read element is not well defined. Could intead eval the full dict.
                continue
            
        #Eval the values
        for key,value in boxDict.iteritems():
            evalList = []
            for i in value:
                evalList.append(eval(i,constantDict))
                boxDictEval[key] = evalList
                
        outfile = open(self.outFile,'w+')
        for key,value in boxDictEval.iteritems():
            outfile.write("%s\t%s\t%s\t%s\n"%(key,value[0],value[1],value[2]))
        outfile.close()

#######################################################################################################################
