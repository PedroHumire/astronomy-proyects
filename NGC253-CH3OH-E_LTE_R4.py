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

#from StringIO import StringIO
#string ="2,6,7,8,20,26,29,30,31,32,33,34,37,38,39,40,41,42,43,44,45,46,53,95,97,101,106,114,122,127,135,136,143,151,153,154,165,172,177,178,179,190,196"
#file = StringIO(string)
##print(file.getvalue())

#======================================================================================#
# Inputs :                                                                             #
#======================================================================================#
# Define the velocity ranges for the lines
Range.unit          = "km/s"   # Possible units are GHz and km/s, MHz, cm-1 and micrometer to be implemented

r2 = Range(83,324)
r9 = Range(171,500)
r19 = Range(206,450)
r20 = Range(163,340)
r21 = Range(172,331)
r22 = Range(215,436)
r23 = Range(224,298)
r24 = Range(43,274)
r25 = Range(48,279)
r26 = Range(54,283)
r31 = Range(170,330)
r52 = Range(171,317)
r54 = Range(177,297)
r56 = Range(188,299)
r59 = Range(217,320)
r66 = Range(0,310)
r70 = Range(120,376)
r71 = Range(50,282)
r77 = Range(191,493)
r85 = Range(165,275)
r86 = Range(0,467)
r87 = Range(43,318)
r90 = Range(0,467)
r93 = Range(204,417)
r98 = Range(191,283)
r99 = Range(0,286)
r100 = Range(0,284)
r101 = Range(145,485)
r104 = Range(211,466)


#my_lines = '"2,6,7,8,20,26,29,30,31,32,33,34,37,38,39,40,41,42,43,44,45,46,53,95,97,101,106,114,122,127,135,136,143,151,153,154,165,172,177,178,179,190,196"'

# Define the step of the walkers values in the MCMC
rpp                 = 400

# Set all what is needed here
sourceName          = "NGC253"
speciesName         = "CH3OH_E"
myModel             = "1C_LTE_R4"
myName              = sourceName+"_"+speciesName+"_"+myModel
cassisDataPath      = Software.getCassisPath()+"/delivery/data/"
cassisScriptsPath   = Software.getCassisPath()+"/delivery/script/examples/"
myDirInput          = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis"
myDirOutput         = "/home/pedro/Documentos/Doctorado/CH3OH_paper/get_spec/CASSIS_analysis/scripts_results/R4/LTE/E/"
inputFile           = myDirInput+"/contsubs/contsub_R4.fus"
outputFile          = myDirOutput+myName+".dat"

# =============================================================================
# USER INPUTS
# =============================================================================

userInputs          = UserInputs(
inputFile           = inputFile,
telescope           = {"2,9,19,20,21,22,23,24,25,26,31,52,54,56,59,66,70,71,77,85,86,87,90,93,98,99,100,101,104": "alma_170m"},
tuningRange         = [84.1573788, 373.078890], #GHz
tuningBand          = 500.0,  #km/s
aijMin              = 1e-6,
eup                 = [0.0, 150.0],
kup                 = ["*","*"],
template            = "Full VASTEL",
moltags             = [32083],
tmb2ta              = False,
tbg                 = 2.73,
continuum           = Software.getCassisPath()+"/delivery/continuum/continuum-0",
isoUnique           = False,
outputFile          = outputFile,
plotTitle           = myName,
warning             = True,

# Enter here the lines and the corresponding ranges to be taken into account in the computation
# The lines should be sorted by increasing frequency.
selectedLines    = {"2":r2,"9":r9,"19":r19,"20":r20,"21":r21,"22":r22,"23":r23,"24":r24,"25":r25,"26":r26,"31":r31,"52":r52,"54":r54,"56":r56,"59":r59,"66":r66,"70":r70,"71":r71,"77":r77,"85":r85,"86":r86,"87":r87,"90":r90,"93":r93,"98":r98,"99":r99,"100":r100,"101":r101,"104":r104},
# rms of the data around each selected line (in K)
#rmsLines = {"2":0.0009,"9":0.0009,"19":0.0011,"20":0.0036,"21":0.0046,"22":0.0046,"23":0.0044,"24":0.0045,"25":0.0044,"26":0.0044,"31":0.0042,"52":0.027,"54":0.004,"56":0.0027,"59":0.0039,"66":0.0038,"70":0.0052,"71":0.0039,"77":0.0047,"85":0.0036,"86":0.0039,"87":0.0044,"90":0.0044,"93":0.0054,"98":0.0063,"99":0.0055,"100":0.0055,"101":0.0064,"104":0.0079},
rmsLines = {"2":0.0009,"9":0.009,"19":0.0011,"20":0.0036,"21":0.0046,"22":0.0046,"23":0.0044,"24":0.0045,"25":0.0044,"26":0.0044,"31":0.0042,"52":0.0027,"54":0.004,"56":0.0027,"59":0.0039,"66":0.0038,"70":0.0052,"71":0.0039,"77":0.0047,"85":0.0036,"86":0.0039,"87":0.0044,"90":0.0044,"93":0.0054,"98":0.0063,"99":0.0055,"100":0.0056,"101":0.0064,"104":0.0079},
#rmsLines = {"6-8,39-154":0.007},
# calibration accuracy of the data around each selected line  
calLines         = {"2,9,19,20,21,22,23,24,25,26,31,52,54,56,59,66,70,71,77,85,86,87,90,93,98,99,100,101,104": 0.15},
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
nmol                = {'min':5E13, 'max':9e15, 'nstep':1, 'log_mode':False},
temp                = {'min':10,  'max':40.0,  'nstep':1, 'log_mode':False},
fwhm                = {'min':30.,    'max':80.0,    'nstep':1, 'log_mode':False},
size                = {'min':1.6,    'max':1.6,    'nstep':1, 'log_mode':False},
vlsr                = {'min':245,    'max':255,    'nstep':1, 'log_mode':False},
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
params_1             = {"nmol": 1.0236E15, "fwhm": 60.214, "vlsr": 251.84, "temp": 30.693, "size" : 1.6}


# set the walker and burning
drawNumber         = 2000
cutOff             = drawNumber/2.
ratioAtCutOff      = 0.5

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
lineModel           = userInputs.plotBestModel(moltag = 32083, overSampling=10, tuningBand = 500, telescope = "alma_170m",lines="2,9,19,20,21,22,23,24,25,26,31,52,54,56,59,66,70,71,77,85,86,87,90,93,98,99,100,101,104")

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
fracOfRejWalkers    = "0.25"
trianglePlot        = [myPython+" "+myPythonScript+" "+myDirOutput+" "+myName+" "+fracOfRejWalkers]
subprocess.Popen(trianglePlot, shell=True)
# ==============================================================================
