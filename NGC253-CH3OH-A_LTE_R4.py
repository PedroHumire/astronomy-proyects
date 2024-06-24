#======================================================================================#
# This script is an example showing how to use the LineAnalysisScripting module with : #
# - one molecule (CN)                                                                  #
# - chi2 computation with the Monte-Carlo Markov Chain algorithm                       #
# - seven components with the radex model                                              #
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

r2                 = Range(191,410)
r4                 = Range(0,500)
r20                = Range(0,400)
r26                = Range(0,500)
r29                = Range(0,500)
r30                = Range(181,391)
r31                = Range(191,356)
r47                = Range(48,464)
r48                = Range(78,388)
r49                = Range(79,392)
r50                = Range(84,401)
r51                = Range(91,406)
r52                = Range(102,417)
r53                = Range(110,430)
r54                = Range(187,500)
r56                = Range(71,400)
r68                = Range(0,500)
r69                = Range(0,303)
r70                = Range(0,500)
r71                = Range(0,500)
r72                = Range(81,400)
r73                = Range(0,500)
r74                = Range(0,500)
r76                = Range(0,371)
r77                = Range(0,500)
r86                = Range(112,500)
r89                = Range(189,329)

# Define the step of the walkers values in the MCMC
rpp                 = 200

# Set all what is needed here
sourceName          = "NGC253"
speciesName         = "CH3OH_A"
myModel             = "1C_LTE_R4"
myName              = sourceName+"_"+speciesName+"_"+myModel
cassisDataPath      = Software.getCassisPath()+"/delivery/data/"
cassisScriptsPath   = Software.getCassisPath()+"/delivery/script/examples/"
myDirInput          = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis"
myDirOutput         = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/scripts_results/R4/LTE/A/"
inputFile           = myDirInput+"/contsubs/contsub_R4.fus"
outputFile          = myDirOutput+myName+".dat"

# =============================================================================
# USER INPUTS
# =============================================================================

userInputs          = UserInputs(
inputFile           = inputFile,
telescope           = {"2,4,20,26,29,30,31,47,48,49,50,51,52,53,54,56,68,69,70,71,72,73,74,76,77,86,89":"alma_170m"},
tuningRange         = [84.1573788, 373.078890], #GHz
tuningBand          = 500.0,  #km/s
aijMin              = 1e-6,
eup                 = [0.0, 150.0],
kup                 = ["*","*"],
template            = "Full VASTEL",
moltags             = [32093],
tmb2ta              = False,
tbg                 = 2.73,
continuum           = Software.getCassisPath()+"/delivery/continuum/continuum-0",
isoUnique           = False,
outputFile          = outputFile,
plotTitle           = myName,
warning             = True,

# Enter here the lines and the corresponding ranges to be taken into account in the computation
# The lines should be sorted by increasing frequency.
selectedLines    = {"2":r2,"4":r4,"20":r20,"26":r26,"29":r29,"30":r30,"31":r31,"47":r47,"48":r48,"49":r49,"50":r50,"51":r51,"52":r52,"53":r53,"54":r54,"56":r56,"68":r68,"69":r69,"70":r70,"71":r71,"72":r72,"73":r73,"74":r74,"76":r76,"77":r77,"86":r86,"89":r89},

# rms of the data around each selected line (in K)
rmsLines = {"2":0.0012,"4":0.0009,"20":0.0053,"26":0.0035,"29":0.0046,"30":0.0026,"31":0.004,"47":0.0053,"48":0.0034,"49":0.0037,"50":0.0038,"51":0.0038,"52":0.0031,"53":0.0033,"54":0.0034,"56":0.0031,"68":0.0049,"69":0.0046,"70":0.0042,"71":0.0054,"72":0.0055,"73":0.0048,"74":0.0034,"76":0.0055,"77":0.0064,"86":0.0049,"89":0.0064},

# calibration accuracy of the data around each selected line  
calLines         = {"2,4,20,26,29,30,31,47,48,49,50,51,52,53,54,56,68,69,70,71,72,73,74,76,77,86,89":0.15},
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
nmol                = {'min':5E14, 'max':4e15, 'nstep':1, 'log_mode':False},
temp                = {'min':18.0,  'max':50.0,  'nstep':1, 'log_mode':False},
fwhm                = {'min':40,    'max':90.0,    'nstep':1, 'log_mode':False},
size                = {'min':1.6,    'max':1.6,    'nstep':1, 'log_mode':False},
vlsr                = {'min':245,    'max':260,    'nstep':1, 'log_mode':False},
iso                 = {'min':1.0,    'max':1.0,    'nstep':1, 'log_mode':False},
interacting         = False,
model               = "lte",
# only needed for RADEX
#collisioner         = ["p_H2","o_H2"],
#n_p_H2              = {'min':2.499999e8, 'max':2.500001e8, 'nstep':1, 'log_mode':True},
#n_o_H2              = {'min':7.499999e8, 'max':7.500001e8,   'nstep':1, 'log_mode':True},
                                #collisionFile       = ["hnco_op-h2_cassis.dat"],
geometry            = "slab",
#reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp,  "n_p_H2": rpp, "n_o_H2": rpp},
)
reducePhysicalParam = {"nmol": rpp, "temp": rpp, "fwhm": rpp, "size": rpp, "vlsr": rpp, "iso": rpp},



# Initialisation of the starting point
params_1             = {"nmol": 9.73E14, "fwhm": 64.641, "vlsr": 252.76, "temp": 34.73, "size" : 1.6}


# set the walker and burning
drawNumber         = 1000
cutOff             = drawNumber/2.
ratioAtCutOff      = 0.5

# Execution time beginning
timeStart = time.time()


# Execution time beginning
timeStart = time.time()

# Computation of the MCMC
#userInputs.initComponentsForMCMC([comp_1, params_1],[comp_2, params_2],[comp_continuum, params_continuum],[comp_3, params_3],[comp_4, params_4])
userInputs.initComponentsForMCMC([comp_1, params_1])

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
lineModel           = userInputs.plotBestModel(moltag = 32093, overSampling=10, tuningBand = 500, telescope = "alma_170m",lines="2,4,20,26,29,30,31,47,48,49,50,51,52,53,54,56,68,69,70,71,72,73,74,76,77,86,89")

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
#myPython            = '/soft/astro/casa/casa-release-5.7.0-134.el7/lib/casa/bin/'
myPython            = "/usr/bin/python3.8"
#myPython            = "/usr/bin/python3"
myPythonScript      = cassisScriptsPath+"Plots_MCMC.py"
# Set the fraction of rejected walkers here !
fracOfRejWalkers    = "0.25"
trianglePlot        = [myPython+" "+myPythonScript+" "+myDirOutput+" "+myName+" "+fracOfRejWalkers]
subprocess.Popen(trianglePlot, shell=True)
# ==============================================================================
