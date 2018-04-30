######################################################################
# motifSearch.py                                                     #
# Author:  Dario Ghersi                                              #
# Version: 20180418                                                  #
# Goal:    A GUI program that calculates gene family distributions   #
#          and clonotypes for CDR regions that match a given motif   #
######################################################################

import re
import sys
from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtGui import *
from gui import * 

######################################################################
# GLOBAL FUNCTIONS                                                   #
######################################################################

def runMotifSearch(data, params):

  ## output variables
  numClonot = 0
  numReads = 0
  numTotClonot = 0
  numTotReads = 0
  numMatchClonot = 0
  numMatchReads = 0
  clonotypes = ""
  output = ""
  
  ## process the header
  fields = data[0].split("\t")

  fieldPos = {}
  fields = data[0][:-1].split("\t")
  count = 0
  for field in fields:
    fieldPos[field] = count
    count += 1

  ## process the data
  famFreqV = {}
  famFreqJ = {}
  alleleFreqV = {}
  alleleFreqJ = {}
  for i in range(1, len(data)):
    fields = data[i][:-1].split("\t")

    if fields[fieldPos["sampleID"]].find(params["sampleID"]) != -1 and\
       (fields[fieldPos["type"]].upper() == params["type"] or\
        params["type"] == "ALL") and\
        (fields[fieldPos["visit"]].upper() == params["visit"] or\
         params["visit"] == "ALL") and\
        fields[fieldPos["epitope"]].upper() == params["epitope"] :

      # clonotype matching everything but not length necessarily
      numTotClonot += 1
      numTotReads += int(fields[fieldPos["reads"]])

      if fields[fieldPos["cdr3_length"]] == params["cdrLen"]:
        numClonot += 1
        numReads += int(fields[fieldPos["reads"]])

        if params["motif"].search(fields[fieldPos["amino_acid"]]):

          numMatchClonot += 1
          numMatchReads += int(fields[fieldPos["reads"]])
        
          # store clonotype data
          clonotypes += " ".join((fields[fieldPos["sampleID"]],
                                  fields[fieldPos["rearrangement"]],
                                  fields[fieldPos["amino_acid"]],
                                  fields[fieldPos["v_family"]],
                                  fields[fieldPos["v_gene"]],
                                  fields[fieldPos["j_family"]],
                                  fields[fieldPos["j_gene"]],
                                  "\n"))

          # calculate V family freq
          if famFreqV.has_key(fields[fieldPos["v_family"]]):
            famFreqV[fields[fieldPos["v_family"]]] += 1
          else:
            famFreqV[fields[fieldPos["v_family"]]] = 1

          # calculate V allele freq
          if alleleFreqV.has_key(fields[fieldPos["v_gene"]]):
            alleleFreqV[fields[fieldPos["v_gene"]]] += 1
          else:
            alleleFreqV[fields[fieldPos["v_gene"]]] = 1          
        
          # calculate J family freq
          if famFreqJ.has_key(fields[fieldPos["j_family"]]):
            famFreqJ[fields[fieldPos["j_family"]]] += 1
          else:
            famFreqJ[fields[fieldPos["j_family"]]] = 1

          if alleleFreqJ.has_key(fields[fieldPos["j_gene"]]):
            alleleFreqJ[fields[fieldPos["j_gene"]]] += 1
          else:
            alleleFreqJ[fields[fieldPos["j_gene"]]] = 1
        
  ## calculate percentages and sort families and genes
  tot = sum(map(lambda x: famFreqV[x], famFreqV.keys()))

  for fam in famFreqV:
    famFreqV[fam] = float(famFreqV[fam]) / tot

  for gene in alleleFreqV:
    alleleFreqV[gene] = float(alleleFreqV[gene]) / tot
    
  for fam in famFreqJ:
    famFreqJ[fam] = float(famFreqJ[fam]) / tot

  for gene in alleleFreqJ:
    alleleFreqJ[gene] = float(alleleFreqJ[gene]) / tot

  for fam, value in sorted(famFreqV.iteritems(), key=lambda (k, v): (v, k),
                           reverse=True):
    if fam == "":
      fam = "N/A"
    output += "%s: %.3f%%" % (fam, 100.0 * value)
    output += "\n"

  output += "\n-------------------\n"

  for gene, value in sorted(alleleFreqV.iteritems(), key=lambda (k, v): (v, k),
                            reverse=True):
    if gene == "":
      gene = "N/A"
    output += "%s: %.3f%%" % (gene, 100.0 * value)
    output += "\n"

  output += "\n-------------------\n"
    
  for fam, value in sorted(famFreqJ.iteritems(), key=lambda (k, v): (v, k),
                           reverse=True):
    if fam == "":
      fam = "N/A"
    output += "%s: %.3f%%" % (fam, 100.0 * value)
    output += "\n"

  output += "\n-------------------\n"

  for gene, value in sorted(alleleFreqJ.iteritems(), key=lambda (k, v): (v, k),
                            reverse=True):
    if gene == "":
      gene = "N/A"
    output += "%s: %.3f%%" % (gene, 100.0 * value)
    output += "\n"
    
  output += "\n********************\n\n"

  ## calculate the number of unique clonotypes and reads matching the
  ## search criteria (by CDR3 length)
  output += "Unique clonotypes (by length):\n%d / %d = %f%%" %\
            (numMatchClonot, numClonot,\
             100.0 * numMatchClonot / numClonot)
  output += "\n\nNumber of reads (by length):\n%d / %d = %f%%" %\
            (numMatchReads, numReads,\
             100.0 * numMatchReads / numReads)
  output += "\n\n********************\n\n"

  ## calculate the number of unique clonotypes and reads matching the
  ## search criteria (total)
  output += "Unique clonotypes (over total):\n%d / %d = %f%%" %\
            (numMatchClonot, numTotClonot,\
             100.0 * numMatchClonot / numTotClonot)
  output += "\n\nNumber of reads (over total):\n%d / %d = %f%%" %\
            (numMatchReads, numTotReads,\
             100.0 * numMatchReads / numTotReads)
  
  return clonotypes, output

######################################################################
# CLASSES                                                            #
######################################################################

class MyApp(QtGui.QMainWindow, Ui_MainWindow):
    
  def __init__(self):
    QtGui.QMainWindow.__init__(self)
    Ui_MainWindow.__init__(self)

    ## set up variables
    self.setupUi(self)
    self.typeBox.addItem("Persistent")
    self.typeBox.addItem("Non-Persistent")
    self.typeBox.addItem("De Novo")
    self.typeBox.addItem("All")
    self.visitBox.addItem("Acute")
    self.visitBox.addItem("Conv")
    self.visitBox.addItem("All")
    self.cdrLength = ""
    self.data = []

    ## go button event
    self.goButton.clicked.connect(self.RunSearch)

    ## file menu
    self.loadData.triggered.connect(self.ReadFile)
    
  ####################################################################
      
  def CheckParameters(self):
    """
    make sure all parameters are properly set
    """
    
    errorMsg = ""
    if self.motif == "":
      errorMsg = "Please enter a motif"

    elif not isinstance(self.cdrLength, int):
      errorMsg = "Please enter the CDR length"

    elif len(self.data) == 0:
      errorMsg = "Please load some data"

    elif self.sample == "":
      errorMsg = "Please enter a sample ID"
      
    if errorMsg != "":
      msg = QMessageBox()
      msg.setText(errorMsg)
      msg.setWindowTitle("Error")
      msg.setStandardButtons(QMessageBox.Ok)
      msg.exec_()
      return False

    return True

  ####################################################################

  def ReadFile(self):
    """
    Store the data into a variable
    """

    name = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
    file = open(name, "r")
    for line in file:
      self.data.append(line)
    
    file.close()

    ## display a pop-up window to confirm that data has been loaded
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setText("File successfully loaded:")
    msg.setInformativeText(name)
    msg.setWindowTitle("MotifSearch")
    msg.setStandardButtons(QMessageBox.Ok)
	
    msg.exec_()

  ####################################################################
    
  def RunSearch(self):

    self.outputClonotypes.setText("")
    self.genesOutput.setText("")
    
    ## get the parameters
    self.motif = self.motifEntry.text()
    self.cdrLength = self.cdrLengthEntry.text()
    try:
      self.cdrLength = int(self.cdrLength)
    except:
      self.cdrLength = ""
    self.epitope = self.epitopeEdit.text()
    self.phase = self.typeBox.currentText()
    self.sample = self.sampleEdit.text()
    self.visit = self.visitBox.currentText()

    ## check that the parameters are properly set
    if self.CheckParameters():
      
      ## run the motif search
      params = {}
      params["motif"] = str(self.motif)
      params["motif"] = params["motif"].replace("X", "[A-Z]")
      params["motif"] = re.compile(str(params["motif"]))
      params["cdrLen"] = str(self.cdrLength)
      params["epitope"] = str(self.epitope).upper()
      params["sampleID"] = str(self.sample)
      params["type"] = str(self.phase).upper()
      params["visit"] = str(self.visit).upper()
      clonotypes, families = runMotifSearch(self.data, params)

      if clonotypes != "":
        self.outputClonotypes.setText(clonotypes)
        self.genesOutput.setText(families)
      
######################################################################
# MAIN PROGRAM                                                       #
######################################################################
    
if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  window = MyApp()
  window.show()
  sys.exit(app.exec_())
