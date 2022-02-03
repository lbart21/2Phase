import os
import matplotlib.pyplot as plt
import numpy as np
from LukesExplicitEulerAlgorithm.ProcessData import dataFileToObject
from LukesExplicitEulerAlgorithm.AnalyticalSods import AnalyticalSodsSolver
from LukesExplicitStiffenedEOSAlgorithm.SIUnitsDictionary import SIUnits
from LukesExplicitEulerAlgorithm.ProcessEilmer import ProcessEilmerData

class GenerateSinglePlots():
    def __init__(self, plotVars, PISODataFileName) -> None:
        """
        convert dataFile to object
        loop through plot vars:
            plot dist of var and save file
        """
        data = dataFileToObject(PISODataFileName)
        self.tFinal = data.tFinal
        self.props = data.props
        
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            formattedTime = '{:.3f}'.format(self.tFinal / 1e-6)
            plt.title("Distribution of " + var + " at t = " \
                                                    + formattedTime + "\u03BCs")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(self.props["pos_x"], self.props[var], marker = '.')
            filename = var + " distribution at t = " + formattedTime + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class GenerateWaterfallPlots():
    """
    take dataFiles and convert each to objects which then become elements of list
    loop through plot vars:
        loop through data object list and plot var from each object
        save file
    """
    def __init__(self, PISODataFileNames, plotVars) -> None:
        dataFromFiles = {}
        t_list = []
        for file in PISODataFileNames:
            PISOData = dataFileToObject(filename = file)
            dataFromFiles[str(PISOData.tFinal)] = PISOData.props
            t_list.append(str(PISOData.tFinal))
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Distribution of " + var + " at multiple time values")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for time in t_list:
                formattedTime = '{:.3f}'.format(float(time) / 1e-6)
                
                plt.scatter(dataFromFiles[time]["pos_x"], dataFromFiles[time][var], \
                    label = "Distribution at t = " + formattedTime + "\u03BCs", \
                        marker = ".")
            plt.legend()
            plt.grid()
            filename = var + " distribution at multiple times.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class CompareToAnalytical():
    def __init__(self, nPointsPerRegion, PISODataFileName, plotVars) -> None:
        """
        
        """
        PISOData = dataFileToObject(filename = PISODataFileName)
        # Strip data for Analytical solver
        dx = PISOData.props["pos_x"][1] - PISOData.props["pos_x"][0]
        L = dx * PISOData.nCells
        t = PISOData.tFinal
        p_L = PISOData.props["p"][0]
        p_R = PISOData.props["p"][-1]
        rho_L = PISOData.props["rho"][0]
        rho_R = PISOData.props["rho"][-1]
        vel_x_L = PISOData.props["vel_x"][0]
        vel_x_R = PISOData.props["vel_x"][-1]
        
        # Generate Analytical data
        self.analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Comparison of Distribution of " + var \
                + " to Analytical Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(PISOData.props["pos_x"], PISOData.props[var], \
                marker = '.', label = "PISO data")
            plt.plot(self.analyticalData.finalProps["pos_x"], self.analyticalData.finalProps[var],\
                label = "Analytical Solution")
            plt.legend()
            plt.grid()
            filename = "PISO to Analytical Comparison of " + var + " at t = " + str(t) + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        
class CompareToEilmer():
    def __init__(self, PISODataFileName, EilmerDataNames, plotVars) -> None:
        self.PISOData = dataFileToObject(filename = PISODataFileName)
        self.EilmerData = ProcessEilmerData(dataFiles = EilmerDataNames)
        """
        
        """
        t = self.PISOData.tFinal
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Comparison of Distribution of " + var \
                + " to Eilmer Data at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(self.PISOData.props["pos_x"], self.PISOData.props[var], \
                marker = ',', label = "PISO data")
            plt.scatter(self.EilmerData.Data["pos_x"], self.EilmerData.Data[var], \
                marker = '+', label = "Eilmer data")
            plt.legend()
            plt.grid()
            filename = "PISO to Eilmer  of " + var + " at t = " + str(t) + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class CompareToAnalyticalAndEilmer():
    def __init__(self, PISODataFileName, EilmerDataNames, nPointsPerRegion, plotVars) -> None:
        """
        
        """
        self.PISOData = dataFileToObject(filename = PISODataFileName)
        # Strip data for Analytical solver
        dx = self.PISOData.props["pos_x"][1] - self.PISOData.props["pos_x"][0]
        L = dx * self.PISOData.nCells
        p_L = self.PISOData.props["p"][0]
        p_R = self.PISOData.props["p"][-1]
        rho_L = self.PISOData.props["rho"][0]
        rho_R = self.PISOData.props["rho"][-1]
        vel_x_L = self.PISOData.props["vel_x"][0]
        vel_x_R = self.PISOData.props["vel_x"][-1]
        t = self.PISOData.tFinal
        

        # Generate Analytical data
        self.analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        # Process EilmerData
        self.EilmerData = ProcessEilmerData(dataFiles = EilmerDataNames)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Comparison of Distribution of " + var \
                + " to Eilmer Data and Analytical Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(self.PISOData.props["pos_x"], self.PISOData.props[var], \
                marker = ',', label = "PISO data")
            plt.scatter(self.EilmerData.Data["pos_x"], self.EilmerData.Data[var], \
                marker = '+', label = "Eilmer data")
            plt.plot(self.analyticalData.finalProps["pos_x"], self.analyticalData.finalProps[var])
            plt.legend()
            plt.grid()
            filename = "PISO, Eilmer and Analytical Comparison of " + var + " at t = " + str(t) + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()


class StackMultipleVariablesPISO():
    def __init__(self, PISODataFileName, plotVars) -> None:
        """
        
        """
        self.PISOData = dataFileToObject(filename = PISODataFileName)
        t = self.PISOData.tFinal
        fig = plt.figure(figsize=(15, 5))
        joinedVarNames = ", ".join(plotVars)
        plt.title("Dimensionless Comparison of Distributions of " + joinedVarNames \
                + " at t = " + str(t) + " s")
        plt.xlabel("Position (m)")
        plt.ylabel("Normalised Property (-)", \
            rotation = "horizontal", ha = "right")
        for var in plotVars:
            plt.scatter(self.PISOData.props["pos_x"], np.array(self.PISOData.props[var]) / np.max(np.abs(self.PISOData.props[var])), \
                marker = ',', label = var)
        plt.legend()
        plt.grid()
        filename = "Dimensionless Comparison of " + joinedVarNames + " at t = " + str(t) + ".jpg"
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()


class StackMultipleVaraiblesPISOWithAnalytical():
    def __init__(self, PISODataFileName, plotVars, nPointsPerRegion) -> None:
        """
        
        """
        self.PISOData = dataFileToObject(filename = PISODataFileName)

        # Strip data for Analytical solver
        dx = self.PISOData.props["pos_x"][1] - self.PISOData.props["pos_x"][0]
        L = dx * self.PISOData.nCells
        t = self.PISOData.tFinal
        p_L = self.PISOData.props["p"][0]
        p_R = self.PISOData.props["p"][-1]
        rho_L = self.PISOData.props["rho"][0]
        rho_R = self.PISOData.props["rho"][-1]
        vel_x_L = self.PISOData.props["vel_x"][0]
        vel_x_R = self.PISOData.props["vel_x"][-1]
        

        # Generate Analytical data
        self.analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        fig = plt.figure(figsize=(15, 5))
        joinedVarNames = ", ".join(plotVars)
        plt.title("Dimensionless Comparison of Distributions of " + joinedVarNames \
                + " to Analytical Solution at t = " + str(t) + " s")
        plt.xlabel("Position (m)")
        plt.ylabel("Normalised Property (-)", \
            rotation = "horizontal", ha = "right")
        for var in plotVars:

            plt.scatter(self.PISOData.props["pos_x"], \
                np.array(self.PISOData.props[var]) \
                / np.max(np.abs(self.analyticalData.finalProps[var])), \
                marker = ',', label = "PISO " + var)

            plt.plot(self.analyticalData.finalProps["pos_x"], \
                np.array(self.analyticalData.finalProps[var]) \
                / np.max(np.abs(self.analyticalData.finalProps[var])), label = "Analytical " + var)
        plt.legend()
        plt.grid()
        filename = "Dimensionless Comparison of " + joinedVarNames + " to Analytical Solution at t = " + str(t) + ".jpg"
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

class StackMultipleVariablesAnalytical():
    def __init__(self, PISODataFileName, nPointsPerRegion, plotVars) -> None:
        """
        
        """
        self.PISOData = dataFileToObject(filename = PISODataFileName)
        # Strip data for Analytical solver
        dx = self.PISOData.props["pos_x"][1] - self.PISOData.props["pos_x"][0]
        L = dx * self.PISOData.nCells
        t = self.PISOData.tFinal
        p_L = self.PISOData.props["p"][0]
        p_R = self.PISOData.props["p"][-1]
        rho_L = self.PISOData.props["rho"][0]
        rho_R = self.PISOData.props["rho"][-1]
        vel_x_L = self.PISOData.props["vel_x"][0]
        vel_x_R = self.PISOData.props["vel_x"][-1]
        

        # Generate Analytical data
        self.analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        fig = plt.figure(figsize=(15, 5))
        joinedVarNames = ", ".join(plotVars)
        plt.title("Analytical Solution Dimensionless Comparison of Distributions of " \
            + joinedVarNames + " at t = " + str(t) + " s")
        plt.xlabel("Position (m)")
        plt.ylabel("Normalised Property (-)", \
            rotation = "horizontal", ha = "right")
        for var in plotVars:

            plt.plot(self.analyticalData.finalProps["pos_x"], \
                np.array(self.analyticalData.finalProps[var]) \
                / np.max(np.abs(self.analyticalData.finalProps[var])), label = "Analytical " + var)
       
        plt.legend()
        plt.grid()
        filename = "Analytical Solution Dimensionless Comparison of " + joinedVarNames + " at t = " + str(t) + ".jpg"
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
            
class CompareFluxCalculators():
    def __init__(self, PISOFileNames, plotVars) -> None:
        """
        fluxCalculators must be in same order as PISOFileNames
        """
        self.PISOData = {}
        t = 0
        fluxCalculatorsList = []
        for fileName in PISOFileNames:
            data = dataFileToObject(filename = fileName)
            self.PISOData[data.fluxCalculator] = data.props
            t += data.tFinal
            fluxCalculatorsList.append(data.fluxCalculator)

        t /= len(PISOFileNames)    # get average final time of all sims, 
                                        # shouldn't give too much error as 
                                        # dt is small (ideally)
        joinedFluxNames = ", ".join(fluxCalculatorsList)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Flux Calculators at t = " + str(t) \
                        + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for FC in fluxCalculatorsList:
                plt.scatter(self.PISOData[FC]["pos_x"], self.PISOData[FC][var], \
                                    label = FC, marker = ".")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedFluxNames + " Schemes.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class CompareFluxCalculatorsToEilmer():
    """
    
    """
    def __init__(self, PISOFileNames, EilmerFileNames, plotVars, nPointsPerRegion) -> None:
        t = 0
        PISOData = {}
        fluxCalculatorsList = []
        nCells = None
        for fileName in PISOFileNames:
            data = dataFileToObject(filename = fileName)
            t += data.tFinal
            PISOData[data.fluxCalculator] = data.props
            nCells = data.nCells
            fluxCalculatorsList.append(data.fluxCalculator)
        t /= len(PISOFileNames)
        joinedFluxCalculatorsList = ", ".join(fluxCalculatorsList)
        ### Generate Analytical Data
        dx = PISOData[fluxCalculatorsList[0]]["pos_x"][1] \
                - PISOData[fluxCalculatorsList[0]]["pos_x"][0]
        L = dx * nCells
        p_L = PISOData[fluxCalculatorsList[0]]["p"][0]
        p_R = PISOData[fluxCalculatorsList[0]]["p"][-1]
        rho_L = PISOData[fluxCalculatorsList[0]]["rho"][0]
        rho_R = PISOData[fluxCalculatorsList[0]]["rho"][-1]
        vel_x_L = PISOData[fluxCalculatorsList[0]]["vel_x"][0]
        vel_x_R = PISOData[fluxCalculatorsList[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        EilmerData = ProcessEilmerData(dataFiles = EilmerFileNames)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Flux Calculators, Eilmer and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for FC in fluxCalculatorsList:
                plt.scatter(PISOData[FC]["pos_x"], PISOData[FC][var], \
                    label = "PISO " + FC + " Scheme", marker = ".")
            plt.scatter(EilmerData.Data["pos_x"], EilmerData.Data[var], marker = "+", label = "Eilmer")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedFluxCalculatorsList + " Schemes to Analytical and Eilmer Solutions.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

        
class CompareFluxCalculatorsToAnalytical():
    """
    
    """
    def __init__(self, PISOFileNames, plotVars, nPointsPerRegion) -> None:
        PISOData= {}
        t = 0
        nCells = None
        fluxCalculatorsList = []
        for fileName in PISOFileNames:
            data = dataFileToObject(filename = fileName)
            PISOData[data.fluxCalculator] = data.props
            nCells = data.nCells
            t += data.tFinal
            fluxCalculatorsList.append(data.fluxCalculator)
        t /= len(PISOFileNames)
        joinedFluxCalculatorsList = ", ".join(fluxCalculatorsList)
        ### Generate Analytical Data
        dx = PISOData[fluxCalculatorsList[0]]["pos_x"][1] \
                - PISOData[fluxCalculatorsList[0]]["pos_x"][0]
        L = dx * nCells
        p_L = PISOData[fluxCalculatorsList[0]]["p"][0]
        p_R = PISOData[fluxCalculatorsList[0]]["p"][-1]
        rho_L = PISOData[fluxCalculatorsList[0]]["rho"][0]
        rho_R = PISOData[fluxCalculatorsList[0]]["rho"][-1]
        vel_x_L = PISOData[fluxCalculatorsList[0]]["vel_x"][0]
        vel_x_R = PISOData[fluxCalculatorsList[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Flux Calculators and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for FC in fluxCalculatorsList:
                plt.scatter(PISOData[FC]["pos_x"], PISOData[FC][var], \
                    label = "PISO " + FC + " Scheme", marker = ".")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedFluxCalculatorsList + " Schemes to Analytical Solution.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class CheckGridConvergence():
    def __init__(self, PISOFileNames, plotVars) -> None:
        """
        loop through PISO files, create data objects for each file and make a 
        dictionary element for each object, key being str(nCells)
        """
        nCells_list = []
        t = 0
        data = {}
        for PISOFile in PISOFileNames:
            PISOData = dataFileToObject(filename = PISOFile)
            nCells_list.append(str(PISOData.nCells))
            t += PISOData.tFinal
            data[str(PISOData.nCells)] = PISOData.props
        t /= len(PISOFileNames)
        joinednCellsList = ", ".join(nCells_list)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Cell Numbers at t = " + str(t) \
                        + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for nCells in nCells_list:
                plt.scatter(data[nCells]["pos_x"], data[nCells][var], \
                            label = nCells, marker = ".")
            plt.grid()
            plt.legend()
            filename = var + " distribution at t = " + str(t) \
                    + " Comparison between " + joinednCellsList + " cells.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        
class CompareGridConvergenceToEilmer():
    def __init__(self, PISOFileNames, EilmerFileNames, nPointsPerRegion, plotVars) -> None:
        PISOnCellsList = []
        t = 0
        PISODataDict = {}
        for PISOFile in PISOFileNames:
            PISOData = dataFileToObject(filename = PISOFile)
            PISOnCellsList.append(str(PISOData.nCells))
            t += PISOData.tFinal
            PISODataDict[str(PISOData.nCells)] = PISOData.props
        t /= len(PISOFileNames)
        joinednCellsList = ", ".join(PISOnCellsList)
        ### Generate Analytical Data
        dx = PISODataDict[PISOnCellsList[0]]["pos_x"][1] \
                - PISODataDict[PISOnCellsList[0]]["pos_x"][0]
        L = dx * int(PISOnCellsList[0])
        p_L = PISODataDict[PISOnCellsList[0]]["p"][0]
        p_R = PISODataDict[PISOnCellsList[0]]["p"][-1]
        rho_L = PISODataDict[PISOnCellsList[0]]["rho"][0]
        rho_R = PISODataDict[PISOnCellsList[0]]["rho"][-1]
        vel_x_L = PISODataDict[PISOnCellsList[0]]["vel_x"][0]
        vel_x_R = PISODataDict[PISOnCellsList[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        EilmerData = {}
        EilmernCellsList = []
        for frontInd in range(0, len(EilmerFileNames), 2):
            data = ProcessEilmerData([EilmerFileNames[frontInd], \
                                            EilmerFileNames[frontInd+1]])
            EilmerData[str(data.Data["nCells"])] = data.Data
            EilmernCellsList.append(str(data.Data["nCells"]))
        
        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Cell Numbers for PISO and Eilmer Simulations and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for PISOnCells in PISOnCellsList:
                plt.scatter(PISODataDict[PISOnCells]["pos_x"], PISODataDict[PISOnCells][var], \
                    label = "PISO " + PISOnCells + " cells", marker = ".")
            for EilmernCells in EilmernCellsList:
                plt.scatter(EilmerData[EilmernCells]["pos_x"], EilmerData[EilmernCells][var], \
                    label = "Eilmer " + EilmernCells + " cells", marker = "+")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinednCellsList + " Cells to Analytical Solution and Eilmer.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class CompareGridConvergenceToAnalytical():
    def __init__(self, plotVars, PISODataFiles, nPointsPerRegion) -> None:
        ### Process PISO data
        PISOData = {}
        nCells_list = []
        t = 0
        for file in PISODataFiles:
            data = dataFileToObject(filename = file)
            nCells_list.append(str(data.nCells))
            PISOData[str(data.nCells)] = data.props
            t += data.tFinal
            
        t /= len(PISODataFiles)
        joinednCellsList = ", ".join(nCells_list)
        ### Generate Analytical Data
        dx = PISOData[nCells_list[0]]["pos_x"][1] \
                - PISOData[nCells_list[0]]["pos_x"][0]
        L = dx * int(nCells_list[0])
        p_L = PISOData[nCells_list[0]]["p"][0]
        p_R = PISOData[nCells_list[0]]["p"][-1]
        rho_L = PISOData[nCells_list[0]]["rho"][0]
        rho_R = PISOData[nCells_list[0]]["rho"][-1]
        vel_x_L = PISOData[nCells_list[0]]["vel_x"][0]
        vel_x_R = PISOData[nCells_list[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different Cell Numbers and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for nCells in nCells_list:
                plt.scatter(PISOData[nCells]["pos_x"], PISOData[nCells][var], \
                    label = "PISO " + nCells + " cells", marker = ".")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinednCellsList + " Cells to Analytical Solution.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        

class CheckCFLConvergence():
    def __init__(self, PISOFileNames, plotVars) -> None:
        """
        
        """
        cfl_list = []
        PISOData = {}
        t = 0
        for file in PISOFileNames:
            data = dataFileToObject(filename = file)
            PISOData[str(data.cfl)] = data.props
            cfl_list.append(str(data.cfl))
            t += data.tFinal   
        t /= len(PISOFileNames)     
        joinedCFLList = ", ".join(cfl_list)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different CFL Values at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for cfl in cfl_list:
                plt.scatter(PISOData[cfl]["pos_x"], PISOData[cfl][var], \
                    marker = ".", label = "cfl = " + cfl)
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedCFLList + " CFL Values.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        

class CompareCFLConvergenceToEilmer():
    def __init__(self, PISOFileNames, EilmerFileNames, plotVars, nPointsPerRegion, EilmerCFLList) -> None:
        PISOCFLList = []
        t = 0
        nCells = None
        PISODataDict = {}
        for PISOFile in PISOFileNames:
            PISOData = dataFileToObject(filename = PISOFile)
            PISOCFLList.append(str(PISOData.cfl))
            t += PISOData.tFinal
            nCells = PISOData.nCells
            PISODataDict[str(PISOData.cfl)] = PISOData.props
        t /= len(PISOFileNames)
        joinedCFLList = ", ".join(PISOCFLList)
        ### Generate Analytical Data
        dx = PISODataDict[PISOCFLList[0]]["pos_x"][1] \
                - PISODataDict[PISOCFLList[0]]["pos_x"][0]
        L = dx * nCells
        p_L = PISODataDict[PISOCFLList[0]]["p"][0]
        p_R = PISODataDict[PISOCFLList[0]]["p"][-1]
        rho_L = PISODataDict[PISOCFLList[0]]["rho"][0]
        rho_R = PISODataDict[PISOCFLList[0]]["rho"][-1]
        vel_x_L = PISODataDict[PISOCFLList[0]]["vel_x"][0]
        vel_x_R = PISODataDict[PISOCFLList[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        EilmerData = {}
        
        for frontInd in range(0, len(EilmerFileNames), 2):
            data = ProcessEilmerData([EilmerFileNames[frontInd], \
                                            EilmerFileNames[frontInd+1]])
            EilmerData[str(EilmerCFLList[frontInd // 2])] = data.Data
            
    
        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different CFL Numbers for PISO and Eilmer Simulations and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for PISOCFL in PISOCFLList:
                plt.scatter(PISODataDict[PISOCFL]["pos_x"], PISODataDict[PISOCFL][var], \
                    label = "PISO " + PISOCFL + " CFL", marker = ".")
            for EilmerCFL in EilmerCFLList:
                plt.scatter(EilmerData[EilmerCFL]["pos_x"], EilmerData[EilmerCFL][var], \
                    label = "Eilmer " + EilmerCFL + " CFL", marker = "+")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedCFLList + " CFL Values to Analytical and Eilmer Solutions.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass


class CompareCFLConvergenceToAnalytical():
    def __init__(self, PISOFileNames, plotVars, nPointsPerRegion) -> None:
        cfl_list = []
        PISOData = {}
        t = 0
        nCells = None
        for file in PISOFileNames:
            data = dataFileToObject(filename = file)
            PISOData[str(data.cfl)] = data.props
            cfl_list.append(str(data.cfl))
            nCells = data.nCells
            t += data.tFinal   
        t /= len(PISOFileNames)     
        joinedCFLList = ", ".join(cfl_list)
        ### Generate Analytical Data
        dx = PISOData[cfl_list[0]]["pos_x"][1] \
                - PISOData[cfl_list[0]]["pos_x"][0]
        L = dx * nCells
        p_L = PISOData[cfl_list[0]]["p"][0]
        p_R = PISOData[cfl_list[0]]["p"][-1]
        rho_L = PISOData[cfl_list[0]]["rho"][0]
        rho_R = PISOData[cfl_list[0]]["rho"][-1]
        vel_x_L = PISOData[cfl_list[0]]["vel_x"][0]
        vel_x_R = PISOData[cfl_list[0]]["vel_x"][-1]
        

        analyticalData = AnalyticalSodsSolver(\
                                        leftState = [p_L, rho_L, vel_x_L], \
                                        rightState = [p_R, rho_R, vel_x_R], \
                                        nPointsPerRegion = nPointsPerRegion, \
                                        L = L, t = t)
        for var in plotVars:
            fig = plt.figure(figsize=(15, 7))
            plt.title("Comparison of Distribution of " + var \
                    + " Between Different CFL Numbers and Analytical " \
                        + "Solution at t = " + str(t) + " s")
            plt.ylabel(var + " (" + SIUnits[var] +")", \
                    rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for CFL in cfl_list:
                plt.scatter(PISOData[CFL]["pos_x"], PISOData[CFL][var], \
                    label = "PISO " + CFL + " CFL", marker = ".")
            plt.plot(analyticalData.finalProps["pos_x"], analyticalData.finalProps[var], label = "Analytical Solution")
            plt.grid()
            plt.legend()
            filename = var + " Distribution at t = " + str(t) \
                    + " Comparison Between " + joinedCFLList + " CFL Values to Analytical Solution.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()


