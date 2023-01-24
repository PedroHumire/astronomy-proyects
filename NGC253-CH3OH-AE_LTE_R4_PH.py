#======================================================================================#
# This script is an example showing how to use the LineAnalysisScripting module with : #
# - one molecule (CN)                                                               #
# - chi2 computation with the Monte-Carlo Markov Chain algorithm                       #
# - seven components with the radex model                                                  #
# Cassis scripting documentation : http://cassis.irap.omp.eu/docs/script/README.html   #
#======================================================================================#

import time
import subprocess
import ScriptEnvironment
from Range import Range
from Component import Component
from FileReader import FileReader
from LineAnalysisScripting import UserInputs
from eu.omp.irap.cassis.properties import Software
from cassisStats import writeStats
from java.io import File
from Plot import Plot

#======================================================================================#
# Inputs :                                                                             #
#======================================================================================#
# Define the velocity ranges for the lines
Range.unit          = "km/s"   # Possible units are GHz and km/s, MHz, cm-1 and micrometer to be implemented

r1 = Range(203,400)   # E-CH3OH
r5 = Range(180,400)   # A-CH3OH
r10 = Range(44,400)   # A-CH3OH
r11 = Range(180,470)  # A-CH3OH
r12 = Range(187,330)  # E-CH3OH
r14 = Range(120,400)  # A-CH3OH
r15 = Range(216,400)  # E-CH3OH
r16 = Range(120,460)  # E-CH3OH
r21 = Range(20,400)  # A-CH3OH
r30 = Range(30,470)  # A-CH3OH
r33 = Range(100,400)  # E-CH3OH (33,35,37) A-CH3OH (34,36)
r38 = Range(120,330)  # E-CH3OH
r40 = Range(20,450)   # E-CH3OH
#r46 = Range(40,400)   # E-CH3OH
r48 = Range(160,427)  # E-CH3OH
r49 = Range(60,280)  # E-CH3OH
r50 = Range(136,410)  # E-CH3OH (50-51)
r53 = Range(90,418)  # E-CH3OH
r63 = Range(39,334)  # A-CH3OH + E-CH3OH
r67 = Range(90,480)   # A-CH3OH + E-CH3OH
r76 = Range(120,460)  # A-CH3OH (76-77)
r79 = Range(67,460)  # A-CH3OH
r80 = Range(16,270)  # E-CH3OH
r81 = Range(180,450)  # E-CH3OH
r83 = Range(180,335)  # E-CH3OH
r84 = Range(130,435)  # E-CH3OH
r85 = Range(120,400)  # E-CH3OH
r86 = Range(176,395)  # A-CH3OH + E-CH3OH
r88 = Range(188,474)  # A-CH3OH
r106 = Range(150,477) # A-CH3OH
r107 = Range(134,262) # A-CH3OH
#r108 = Range(170,320) # A-CH3OH
r111 = Range(116,300)  # A-CH3OH
r121 = Range(20,370)  # A-CH3OH
r124 = Range(116,370) # E-CH3OH (124-125)
r127 = Range(80,330)  # E-CH3OH 
r130 = Range(130,296) # A-CH3OH 
r131 = Range(34,447) # A-CH3OH 
r134 = Range(165,450) # A-CH3OH + E-CH3OH 
r140 = Range(70,430)  # A-CH3OH + E-CH3OH 
r150 = Range(210,400) # A-CH3OH 
r153 = Range(120,400) # E-CH3OH (153) + A-CH3OH (154)
r155 = Range(120,300) # A-CH3OH 
r156 = Range(120,400) # A-CH3OH 
r158 = Range(120,400) # A-CH3OH 
r159 = Range(120,320) # E-CH3OH 
r160 = Range(120,400) # A-CH3OH 
r161 = Range(225,400) # E-CH3OH
r162 = Range(120,400) # E-CH3OH (162, 164) + A-CH3OH (163)
r167 = Range(224,380) # A-CH3OH 
r168 = Range(120,400) # A-CH3OH 
r171 = Range(137,500) # E-CH3OH + A-CH3OH 
r182 = Range(50,400)  # E-CH3OH + A-CH3OH 
r186 = Range(176,400) # A-CH3OH 
r190 = Range(192,338) # A-CH3OH 
r191 = Range(196,350) # E-CH3OH 
r192 = Range(140,400) # E-CH3OH 
r193 = Range(225,400) # E-CH3OH 

# Define the step of the walkers values in the MCMC
rpp                 = 200

# Set all what is needed here
sourceName          = "NGC253"
speciesName         = "CH3OH_AE"
myModel             = "2C_LTE_R4_PH"
myName              = sourceName+"_"+speciesName+"_"+myModel
cassisDataPath      = Software.getCassisPath()+"/delivery/data/"
cassisScriptsPath   = Software.getCassisPath()+"/delivery/script/examples/"
myDirInput          = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/contsubs/"
myDirOutput         = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/scripts_results/R4/LTE/"
inputFile           = myDirInput+"contsub_R4.fus"
outputFile          = myDirOutput+myName+".dat"

# =============================================================================
# USER INPUTS
# =============================================================================

userInputs          = UserInputs(
inputFile           = inputFile,
telescope           = {"1,5,10,11,12,14-16,21,30,34,35,37,38,40,48-50,53,63,76-77,79-81,83,84,86,88,106-107,111,124-125,130-131,134,140,150,153-156,158-164,167,171,182,186,190-193":"alma_170m"},
#telescope           = {"6-8,39-154": "alma_209m"},
tuningRange         = [84.1573788, 373.078890], #GHz
tuningBand          = 500.0,  #km/s
aijMin              = 1.0e-06,
eup                 = [0.0, 150.0],
kup                 = ["*","*"],
template            = "Full VASTEL",
moltags             = [32083,32093],
tmb2ta              = False,
tbg                 = 2.73,
continuum           = Software.getCassisPath()+"/delivery/continuum/continuum-0",
isoUnique           = False,
outputFile          = outputFile,
plotTitle           = myName,
warning             = True,

# Enter here the lines and the corresponding ranges to be taken into account in the computation
# The lines should be sorted by increasing frequency.
selectedLines    = {"1":r1,"5":r5,"10":r10,"11":r11,"12":r12,"14":r14,"15":r15,"16":r16,"21":r21,"30":r30,"34,35,37":r33,"38":r38,"40":r40,"48":r48,"49":r49,"50":r50,"53":r53,"63":r63,"76-77":r76,"79":r79,"80":r80,"81":r81,"83":r83,"84":r84,"86":r86,"88":r88,"106":r106,"107":r107,"111":r111,"124-125":r124,"130":r130,"131":r131,"134":r134,"140":r140,"150":r150,"153-154":r153,"155":r155,"156":r156,"158":r158,"159":r159,"160":r160,"161":r161,"162-164":r162,"167":r167,"171":r171,"182":r182,"186":r186,"190":r190,"191":r191,"192":r192,"193":r193},
#selectedLines    = {"6-8,39-154":r11},
# rms of the data around each selected line (in K)
rmsLines = {"1,5,10,11,12,14-16,21,30,34,35,37,38,40,48-50,53,63,76-77,79-81,83,84,86,88,106-107,111,124-125,130-131,134,140,150,153-156,158-164,167,171,182,186,190-193":0.005},
#rmsLines = {"6-8,39-154":0.007},
# calibration accuracy of the data around each selected line  
calLines         = {"1,5,10,11,12,14-16,21,30,34,35,37,38,40,48-50,53,63,76-77,79-81,83,84,86,88,106-107,111,124-125,130-131,134,140,150,153-156,158-164,167,171,182,186,190-193": 0.15},
#calLines         = {"6-8,39-154":0.2},
)

#======================================================================================#
# MODEL INPUTS                                                                         #
#======================================================================================#

# Type of models :
# LTE  : nmol, tex, fwhm, size, vlsr and iso if there are different moltags
# RADEX: nmol, collisionFile, n_collisioners, tkin, fwhm, size, vlsr and iso

#======================================================================================#
# Parameters to be defined for all components                                          #
#======================================================================================#

comp_1              = Component(
# Needed for any model
nmol                = {'min':1E14, 'max':1e16, 'nstep':1, 'log_mode':False},
temp                = {'min':10,   'max':100.0,'nstep':1, 'log_mode':False},
fwhm                = {'min':30.,  'max':90.0, 'nstep':1, 'log_mode':False},
size                = {'min':1.6,  'max':1.6,  'nstep':1, 'log_mode':False},
vlsr                = {'min':245,  'max':255,  'nstep':1, 'log_mode':False},
iso                 = {'min':0.1,  'max':10.0, 'nstep':1, 'log_mode':False},
interacting         = True,
model               = "lte",
# only needed for RADEX
#collisioner         = ["p_H2","o_H2"],
#n_p_H2              = {'min':2.499999e8, 'max':2.500001e8, 'nstep':1, 'log_mode':True},
#n_o_H2              = {'min':7.499999e8, 'max':7.500001e8,   'nstep':1, 'log_mode':True},
                                #collisionFile       = ["hnco_op-h2_cassis.dat"],
geometry            = "slab",
reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp},
)

comp_2              = Component(
# Needed for any model
nmol                = {'min':1E14, 'max':1e16, 'nstep':1, 'log_mode':False},
temp                = {'min':10,   'max':100.0,'nstep':1, 'log_mode':False},
fwhm                = {'min':10.,  'max':70.0, 'nstep':1, 'log_mode':False},
size                = {'min':1.6,  'max':1.6,  'nstep':1, 'log_mode':False},
vlsr                = {'min':245,  'max':250,  'nstep':1, 'log_mode':False},
iso                 = {'min':0.1,  'max':10.0, 'nstep':1, 'log_mode':False},
interacting         = True,
model               = "lte",
# only needed for RADEX
#collisioner         = ["p_H2","o_H2"],
#n_p_H2              = {'min':2.499999e8, 'max':2.500001e8, 'nstep':1, 'log_mode':True},
#n_o_H2              = {'min':7.499999e8, 'max':7.500001e8,   'nstep':1, 'log_mode':True},
                                #collisionFile       = ["hnco_op-h2_cassis.dat"],
geometry            = "slab",
reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp},
)

# Initialisation of the starting point
params_1             = {"nmol": 8E14, "fwhm": 83.5, "vlsr": 252, "temp": 26, "size" : 1.6, "iso" : 1.6}
params_2             = {"nmol": 4E14, "fwhm": 32, "vlsr": 248, "temp": 35, "size" : 1.6, "iso" : 1.5}


# set the walker and burning
drawNumber         = 600
cutOff             = drawNumber/2.
ratioAtCutOff      = 0.5

# Execution time beginning
timeStart = time.time()

# Computation of the MCMC
#userInputs.initComponentsForMCMC([comp_1, params_1],[comp_2, params_2],[comp_continuum, params_continuum],[comp_3, params_3],[comp_4, params_4])
userInputs.initComponentsForMCMC([comp_1, params_1],[comp_2, params_2])

# Computation of the minimum chi2
userInputs.computeChi2MinUsingMCMC(drawNumber, cutOff, ratioAtCutOff)

# =============================================================================
# ANALYSIS OF THE RESULTS
# =============================================================================

# A. Plot the best model and save the corresponding spectra and config files
#lineModel = userInputs.plotBestModel(moltag = [43511], overSampling=5.0, tuningBand = 20, telescope = "apex")
#lineModel.saveConfig(File(myDirOutput+myName+".lam"))
#bestPhysicalModels = userInputs.getBestPhysicalModels()
#userInputs.saveBestPhysicalModels(myDirOutput+myName+"_bestModel.lis")

# A. Plot the best model and save the corresponding spectra and config files
lineModel           = userInputs.plotBestModel(moltag = [32083,32093], overSampling=10, tuningBand = 500, telescope = "alma_170m",lines="1,5,10,11,12,14-16,21,30,34,35,37,38,40,48-50,53,63,76-77,79-81,83,84,86,88,106-107,111,124-125,130-131,134,140,150,153-156,158-164,167,171,182,186,190-193")

lineModel.saveConfig(File(myDirOutput+myName+".lam"))
bestPhysicalModels  = userInputs.getBestPhysicalModels()
userInputs.saveBestPhysicalModels(myDirOutput+myName+"_bestModel.lis")

# B. Compute and write the statistics of the parameters
# (for more details, print the documentation with 'print writeStats.__doc__')
writeStats(userInputs, sigmaClip = 3)

# C. Plot the acceptance rates
# (print documentation for "Plot" with "print Plot.__doc__", documentation for "FileReader" coming soon)
reader = FileReader(userInputs.outputFile)
reader.read_file(lines_to_skip=1)
columns = reader.columns
iterations = columns['#']
rates = columns['rate']
sp = Plot(sourceName)
sp.plot(x=iterations, y=rates, xlabel="Iteration number", ylabel=" ", legend="Acceptance rate", lineStyle="line", lineWidth=2, plotType="line", lineColor = "blue")

# CASSIS Execution time ending
timeEnd = time.time()
print "CASSIS execution time =", timeEnd - timeStart

# D. Launch the triangle plot in python
# Set the correct path for your python here !
#myPython            = "/Users/antonio/anaconda/bin/python3.6"
#myPython            = "/Users/antonio/anaconda/bin/python"
myPython            = "/usr/bin/python3.8"
myPythonScript      = cassisScriptsPath+"Plots_MCMC.py"
# Set the fraction of rejected walkers here !
fracOfRejWalkers    = "0.1"
trianglePlot        = [myPython+" "+myPythonScript+" "+myDirOutput+" "+myName+" "+fracOfRejWalkers]
subprocess.Popen(trianglePlot, shell=True)
# ==============================================================================
