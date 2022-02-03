import numpy as np
import copy
from LukesExplicitStiffenedEOSAlgorithm.Integrate import Integrate
import os

class Solver():
    def __init__(self, GEO, currentProps, dt_init, cfl_max, tFinal, \
        fluxCalculator, BCs, cfl_criteria, dataSaveDt, reconstruction, \
            sourceTerms, energyMethod) -> None:
        self.fluxCalculator = fluxCalculator
        self.GEO = copy.deepcopy(GEO)
        self.currentProps = copy.deepcopy(currentProps)
        #print(self.currentProps["phi"].tolist())
        self.LargePropsDictionary = {
        }
        
        self.cfl_criteria = cfl_criteria
        self.cfl = cfl_max
        self.dt_fixed = dt_init
        
        t = 0.0
        self.addData(propDictionary = currentProps, GEO = self.GEO, time = t, \
            energyMethod = energyMethod)
        loop = 0
        addedData = False
        #print(currentProps)
        time_tol = 1e-10
        tPlotTol = 1e-10
        tPlot = dataSaveDt
        self.t_list = [t]
        while t < tFinal and abs(t-tFinal) > time_tol:
            if self.cfl_criteria:
                U = self.currentProps["vel_x"]
                if BCs[0][0] == "InFlow":
                    U = np.hstack((np.array([BCs[0][1][2]]), U))
                if BCs[1][0] == "InFlow":
                    U = np.hstack((np.array([BCs[1][1][2]]), U))
                
                uPlusC = np.abs(U) + np.max(self.currentProps["a_L"])
                maxUPlusC = max(np.abs(uPlusC))
                cfl_dt = cfl_max * self.GEO["dx"][0] / maxUPlusC
                dt = cfl_dt
                

            else:
                dt = dt_init
            newProps = Integrate(dt = dt, GEO = copy.deepcopy(self.GEO), \
                                currentProps = self.currentProps, \
                                fluxCalculator = fluxCalculator, BCs = BCs, \
                                reconstruction = reconstruction, \
                                sourceTerms = sourceTerms, \
                                energyMethod = energyMethod).newProps
            t += dt
            print(t)
            loop += 1
            addedData = False
            
            if abs(t - tPlot) < tPlotTol or t >= tPlot:
                print("Writing Data")
                self.addData(propDictionary = copy.deepcopy(newProps), \
                    GEO = self.GEO, time = t, energyMethod = energyMethod)
                addedData = True
                self.t_list.append(t)
                tPlot += dataSaveDt
            self.currentProps = copy.deepcopy(newProps)
        #print(newProps["phi"].tolist())
        if addedData == False:
            self.addData(propDictionary = newProps, GEO = self.GEO, time = t, \
                energyMethod = energyMethod)
            self.t_list.append(t)
        self.writeToTXT()
        

    def addData(self, propDictionary, GEO, time, energyMethod):
        
        ### add averaged properties
        psi = propDictionary["psi"]
        propDictionary["rho_L_psi"] = np.multiply(psi, propDictionary["rho_L"])
        propDictionary["rho"] = np.multiply(propDictionary["rho_L"], psi) \
                                + np.multiply(propDictionary["rho_G"], 1 - psi)
        propDictionary["gamma_m"] = np.multiply(psi, propDictionary["gamma_L"]) + np.multiply(1-psi, propDictionary["gamma_G"])
        #propDictionary["gamma_m"] = (np.divide(psi, propDictionary["gamma_L"] - 1) \
            #+ np.divide(1 - psi, propDictionary["gamma_G"] - 1)) ** -1 + 1
        propDictionary["a"] = np.multiply(psi, propDictionary["a_L"]) \
            + np.multiply(1 - psi, propDictionary["a_G"])
        propDictionary["Ma"] = np.divide(propDictionary["vel_x"], propDictionary["a"])
        if energyMethod == 'Adapted':
            propDictionary["u"] = np.divide(propDictionary["p"], propDictionary["gamma_m"] - 1) \
                + 0.5 * np.multiply(\
                            np.multiply(propDictionary["vel_x"], propDictionary["vel_x"]), \
                        propDictionary["rho"])
        elif energyMethod == "Normal":
            propDictionary["u"] = np.multiply(np.multiply(psi, propDictionary["rho_L"]), propDictionary["u_L"]) \
                + np.multiply(np.multiply(1 - psi, propDictionary["rho_G"]), propDictionary["u_G"]) \
                + 0.5 * np.multiply(\
                            np.multiply(propDictionary["vel_x"], propDictionary["vel_x"]), \
                        propDictionary["rho"])
        propDictionary["h"] = np.multiply(psi, propDictionary["h_L"]) + np.multiply(1-psi, propDictionary["h_G"])
        propDictionary.update(GEO) #merge prop and GEO dictionaries
        self.LargePropsDictionary[str(time)] = propDictionary
    
    def writeToTXT(self):
        dataSavedTimes = self.LargePropsDictionary.keys()
        cwd = os.getcwd()
        fluxCalculator = self.fluxCalculator
        nCells = self.GEO["nCells"]
        for time in dataSavedTimes:
            float_time = float(time)
            variableNames = list(self.LargePropsDictionary[time].keys())
            variableNames.remove("nCells")
            
            if self.cfl_criteria:
                fileName = "dataAt" + str(format(float_time, ".9f")) \
                    + "Using" + fluxCalculator + str(nCells) + "cells" \
                    + str(self.cfl) + "cfl.txt"
            else:
                fileName = "dataAt" + str(format(float_time, ".9f")) \
                    + "Using" + fluxCalculator + str(nCells) + "cells" \
                    + str(self.dt_fixed) + "dt.txt"
            file = open(cwd + "/data/" + fileName, "w")
            file.write("time: " + time + "\n")
            if self.cfl_criteria:
                file.write("cfl: " + str(format(self.cfl, '.9f')) + "\n")
            else:
                file.write("dt: " + str(format(self.dt_fixed, '.9f')) + "\n")
            file.write("FluxCalculator: " + self.fluxCalculator + "\n")
            file.write("nCells: " + str(self.GEO["nCells"]) + "\n")
            file.write("Variables: " + str(len(variableNames)) + "\n")
            file.write(' '.join(variableNames) + "\n")
            for cell in range(self.GEO["nCells"]):
                for name in variableNames:
                    if name == "nCells":
                        pass
                    else:
                        file.write(str(format(self.LargePropsDictionary[time][name][cell], ".9f")) + " ")
                file.write("\n")
            file.close()