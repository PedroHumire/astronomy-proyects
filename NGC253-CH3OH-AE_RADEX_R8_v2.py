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
#r4   = Range(150,350)
r5   = Range(50,350)
r10  = Range(50,350)
#r11  = Range(50,350)
r27  = Range(50,550)
#r28  = Range(50,350)
r29  = Range(61,355)    #new
r37  = Range(-100,450)
r38  = Range(115,280)   #new
r44  = Range(20,400)    #new 
r54  = Range(0,400)
r77  = Range(-50,450)
#r87  = Range(50,350)
r90  = Range(-200,600)
r97  = Range(50,400) #new 
r98  = Range(50,400)  #new
r100 = Range(50,350)
#r107 = Range(-150,550)
r118 = Range(0,400)   #new
r119 = Range(50,350)
r120 = Range(0,400)  #new
r121 = Range(0,400)  #new
r122 = Range(50,350)
r123 = Range(50,400)   #new
r124 = Range(0,400)   #new
#r135 = Range(-200,600)
r140 = Range(122,295)
r142 = Range(100,400)   #new


# Define the step of the walkers values in the MCMC
rpp                 = 50

# Set all what is needed here
sourceName          = "NGC253"
speciesName         = "CH3OH_AE"
myModel             = "2C_RADEX_R8"
myName              = sourceName+"_"+speciesName+"_"+myModel
cassisDataPath      = Software.getCassisPath()+"/delivery/data/"
cassisScriptsPath   = Software.getCassisPath()+"/delivery/script/examples/"
myDirInput          = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/contsubs/"
myDirOutput         = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/scripts_results/R8/RADEX/"
inputFile           = myDirInput+"contsub_R8.fus"
outputFile          = myDirOutput+myName+".dat"

# =============================================================================
# USER INPUTS
# =============================================================================

userInputs          = UserInputs(
inputFile           = inputFile,
telescope           = {"5,10,27,29,37,38,44,54,77,90,97,98,100,118,119,120,121,122,123,124,140,142":"alma_170m"}, # 42
tuningRange         = [84.1573788, 373.078890], #GHz
tuningBand          = 500.0,  #km/s
aijMin              = 1.0e-6,
eup                 = [0.0, 100.0],
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
selectedLines    = {"5":r5,"10":r10,"27":r27,"29":r29,"37":r37,"38":r38,"44":r44,"54":r54,"77":r77,"90":r90,"97":r97,"98":r98,"100":r100,"118":r118,"119":r119,"120":r120,"121":r121,"122":r122,"123":r123,"124":r124,"140":r140,"142":r142},
# rms of the data around each selected line (in K)
rmsLines         = {"5,10,27,29,37,38,44,54,77,90,97,98,100,118,119,120,121,122,123,124,140,142":0.005},
# calibration accuracy of the data around each selected line  
calLines         = {"5,10,27,29,37,38,44,54,77,90,97,98,100,118,119,120,121,122,123,124,140,142": 0.15},
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
nmol                = {'min':1e13, 'max':3e15, 'nstep':1, 'log_mode':False},
temp                = {'min':10.0, 'max':40.0, 'nstep':1, 'log_mode':False},
fwhm                = {'min':70.,  'max':90.0, 'nstep':1, 'log_mode':False},
size                = {'min':1.6,  'max':1.6,  'nstep':1, 'log_mode':False},
vlsr                = {'min':190,  'max':210,  'nstep':1, 'log_mode':False},
iso                 = {'min':0.1,  'max':10.0, 'nstep':1, 'log_mode':False},
interacting         = True,
model               = "radex",
# only needed for RADEX
collisioner         = ["p_H2"],
n_p_H2              = {'min':1.0e6, 'max':1.0e9, 'nstep':1, 'log_mode':False},
#n_o_H2              = {'min':7.499999e8, 'max':7.500001e8,   'nstep':1, 'log_mode':True},
collisionFile       = ["e-CH3OH-pH2.dat","a-CH3OH-pH2.dat"],
#collisionFile       = ["e-ch3oh-ph2.dat","a-ch3oh-ph2.dat"],
geometry            = "slab",
reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp, "n_p_H2": rpp},
)

comp_2              = Component(
# Needed for any model
nmol                = {'min':1e14, 'max':1e16, 'nstep':1, 'log_mode':False},
temp                = {'min':10.0, 'max':20.0, 'nstep':1, 'log_mode':False},
fwhm                = {'min':40,   'max':80.0, 'nstep':1, 'log_mode':False},
size                = {'min':1.6,  'max':1.6,  'nstep':1, 'log_mode':False},
vlsr                = {'min':190,  'max':210,  'nstep':1, 'log_mode':False},
iso                 = {'min':0.1,  'max':10.0, 'nstep':1, 'log_mode':False},
interacting         = True,
model               = "radex",
# only needed for RADEX
collisioner         = ["p_H2"],
n_p_H2              = {'min':1.0e4, 'max':1.0e6, 'nstep':1, 'log_mode':False},
#n_o_H2              = {'min':7.499999e8, 'max':7.500001e8,   'nstep':1, 'log_mode':True},
#collisionFile       = ["e-ch3oh-ph2.dat","a-ch3oh-ph2.dat"],
collisionFile       = ["e-CH3OH-pH2.dat","a-CH3OH-pH2.dat"],
geometry            = "slab",
reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp, "n_p_H2": rpp},
)

# Initialisation of the starting point
#params_1             = {"nmol": 9.4e13, "fwhm": 76, "vlsr": 204, "temp": 23.0, "size" : 1.6, "iso" : 2.0, "n_p_H2": 2.5e8}
params_1             = {"nmol": 1.1082E14, "fwhm": 75.809, "vlsr": 204.61, "temp": 24.802, "size" : 1.6, "iso" : 1.9841, "n_p_H2": 2.7772e08}
params_2             = {"nmol": 1.2905E15, "fwhm": 62.286, "vlsr": 199.68, "temp": 14.358, "size" : 1.6, "iso" : 3.141, "n_p_H2": 5.949e+05}

# set the walker and burning
drawNumber         = 1000
cutOff             = drawNumber/2.
ratioAtCutOff      = 1.0

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
#lineModel           = userInputs.plotBestModel(moltag = [32083,32093], overSampling=10, tuningBand = 500, telescope = "alma_170m",lines="4,5,10,11,27-28,37,54,77,87,90,100,107,119,122,135")

lineModel           = userInputs.plotBestModel(moltag = [32083,32093], overSampling=10, tuningBand = 500, telescope = "alma_170m",lines="5,10,27,29,37,38,44,54,77,90,97,98,100,118,119,120,121,122,123,124,140,142")

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
myPython            = "/usr/bin/python3.8"
myPythonScript      = cassisScriptsPath+"Plots_MCMC.py"
# Set the fraction of rejected walkers here !
fracOfRejWalkers    = "0.25"
trianglePlot        = [myPython+" "+myPythonScript+" "+myDirOutput+" "+myName+" "+fracOfRejWalkers]
subprocess.Popen(trianglePlot, shell=True)
# ==============================================================================

