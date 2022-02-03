from syslog import closelog
import numpy as np
import copy
from LukesExplicitStiffenedEOSAlgorithm.EOSElementCoefficients import coefficients
import LukesExplicitStiffenedEOSAlgorithm.LiquidEOS as LiquidEOS
from eilmer.gas import GasModel, GasState
from LukesExplicitEulerAlgorithm.FluxCalculators import AUSMPlusUpBASIC
import math as m


class Integrate():
    def __init__(self, GEO, currentProps, dt, fluxCalculator, reconstruction, sourceTerms, BCs, energyMethod) -> None:
        self.GEO = GEO
        self.currentProps = copy.deepcopy(currentProps)
        self.newProps = {}
        self.dt = dt

        self.GEO["A_c"] = np.hstack((np.array([self.GEO["A_c"][0], self.GEO["A_c"][0]]), \
                                self.GEO["A_c"], np.array([self.GEO["A_c"][-1], self.GEO["A_c"][-1]])))
        self.GEO["A_s"] = np.hstack((np.array([np.nan, np.nan]), \
                                self.GEO["A_s"], np.array([np.nan, np.nan])))
        self.GEO["dV"] = np.hstack((np.array([self.GEO["dV"][0], self.GEO["dV"][0]]), \
                                self.GEO["dV"], np.array([self.GEO["dV"][-1], self.GEO["dV"][-1]])))
        self.GEO["dx"] = np.hstack((self.GEO["dx"][1::-1], \
                                self.GEO["dx"], self.GEO["dx"][:-3:-1]))
        self.currentProps = self.calculateAverageProperties(dictionary = self.currentProps, energyMethod = energyMethod)
        self.currentProps = self.addBoundaryCells(dictionary = self.currentProps, BCs = BCs)
        
        self.solveForConservedProps(\
            a = self.currentProps["a"], \
                liquid_a = self.currentProps["a_L"], \
            vel_x = self.currentProps["vel_x"], \
            rho = self.currentProps["rho"], \
            rho_Liquid = self.currentProps["rho_L"], \
            p = self.currentProps["p"], \
            u = self.currentProps["u"], \
            Ma = self.currentProps["Ma"], \
            A_c = self.GEO["A_c"], \
            phi = self.currentProps["phi"], \
            psi = self.currentProps["psi"], \
            fluxCalculator = fluxCalculator, \
            reconstruction = reconstruction, \
            sourceTerms = sourceTerms, \
            EOS = None, \
            energyMethod = energyMethod)
        
        self.newProps = self.updateNewProps(dictionary = self.newProps)

    def addBoundaryCells(self, dictionary, BCs):
        
        entryNames = dictionary.keys()
        if BCs[0][0] == "OutFlow": #zero gradient
            for name in entryNames:
                if name == "phi":
                    point1 = 3 * dictionary["phi"][0] - 2 * dictionary["phi"][1]
                    point2 = 2 * dictionary["phi"][0] - dictionary["phi"][1]
                    
                    dictionary["phi"] = np.hstack((np.array([point1, point2]), \
                                dictionary["phi"]))
                    psi1 = self.calculatePsi(phi = point1, dx = self.GEO["dx"][0])
                    psi2 = self.calculatePsi(phi = point2, dx = self.GEO["dx"][1])
                    
                    dictionary["psi"] = np.hstack((np.array([psi1, psi2]), \
                                dictionary["psi"]))
                elif name == "psi":
                    pass
                else:
                    dictionary[name] = np.hstack((dictionary[name][1::-1], 
                                        dictionary[name]))
                         
        elif BCs[0][0] == "InFlow": #Fixed value
            gmodel = GasModel('ideal-air-gas-model.lua')
            inletState = GasState(gmodel)
            inletState.p = BCs[0][1][0]
            inletState.rho = BCs[0][1][1]
            inletState.update_thermo_from_rhop()
            inletState.update_sound_speed()
            inletState.update_trans_coeffs()
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass

                elif name == "p":
                        dictionary[name] = [BCs[0][1][0]] * 2 \
                                + dictionary[name] 
                                
                elif name == "rho":
                        dictionary[name] = [BCs[0][1][1]] * 2 \
                                + dictionary[name] 
                               
                elif name == "vel_x":
                        dictionary[name] = [BCs[0][1][2]] * 2 \
                                + dictionary[name] 
                                
                elif name == "T":
                        dictionary[name] = [inletState.T] * 2 \
                                + dictionary[name] 
                                
                elif name == "a":
                        dictionary[name] = [inletState.a] * 2 \
                                + dictionary[name] 
                                
                elif name == "Ma":
                        dictionary[name] = \
                            [BCs[0][1][2] / inletState.a] * 2 \
                                + dictionary[name] 
                                
                elif name == "u":
                        dictionary[name] = [inletState.u] * 2 \
                                + dictionary[name] 
                               
                elif name == "h":
                        dictionary[name] = \
                            [gmodel.enthalpy(inletState)] * 2 \
                                + dictionary[name] 
                               
                elif name == "s":
                        dictionary[name] = [gmodel.entropy(inletState)] * 2 \
                                + dictionary[name] 
                                
                elif name == "Cv":
                        dictionary[name] = [gmodel.Cv(inletState)] * 2\
                                + dictionary[name] 
                
                elif name == "Cp":
                        dictionary[name] = [gmodel.Cp(inletState)] * 2\
                                + dictionary[name]

                elif name == "k":
                        dictionary[name] = [inletState.k] * 2\
                                + dictionary[name]

                elif name == "mu":
                        dictionary[name] = [inletState.mu] * 2\
                                + dictionary[name]

                else: 
                        #For epsilon, mu, assume fixed value 
                        # (doesn't matter since these terms won't be used)
                        dictionary[name] = \
                                  [dictionary[name][0]]*2 + dictionary[name] 
                                
        elif BCs[0][0] == "Wall":
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass

                elif name == "vel_x":
                    dictionary[name] = \
                        (-np.array(dictionary[name][1::-1])).tolist() \
                            + dictionary[name]

                else:
                    dictionary[name] = dictionary[name][1::-1] \
                                            + dictionary[name]
            
        elif BCs[0][0] == "FixedP":
            dictionary["p"] = [BCs[0][1][0]] * 2 + dictionary["p"]
            dictionary["T"] = dictionary["T"][1::-1] \
                                            + dictionary["T"]
            dictionary["vel_x"] = dictionary["vel_x"][1::-1] \
                                            + dictionary["vel_x"]                            
            gmodel = GasModel('ideal-air-gas-model.lua')
            GCellState = GasState(gmodel)
            
            for GhostCell in range(2):
                GCellState.p = dictionary["p"][1 - GhostCell]
                GCellState.T = dictionary["T"][1 - GhostCell]
                GCellState.update_thermo_from_pT()
                GCellState.update_sound_speed()
                GCellState.update_trans_coeffs()
                dictionary["rho"] = [GCellState.rho] + dictionary["rho"]
                dictionary["u"] = [GCellState.u] + dictionary["u"]
                dictionary["a"] = [GCellState.a] + dictionary["a"]
                dictionary["h"] = [gmodel.enthalpy(GCellState)] \
                                                + dictionary["h"]
                dictionary["s"] = [gmodel.entropy(GCellState)] + dictionary["s"]
                dictionary["Cv"] = [gmodel.Cv(GCellState)] + dictionary["Cv"]
                dictionary["Cp"] = [gmodel.Cp(GCellState)] + dictionary["Cp"]
                dictionary["k"] = [GCellState.k] + dictionary["k"]
                dictionary["mu"] = [GCellState.mu] + dictionary["mu"]
                dictionary["Ma"] = \
                    [dictionary["vel_x"][1 - GhostCell] / GCellState.a] \
                            + dictionary["Ma"]

        elif BCs[0][0] == "FromStag":#Calculate P and T from inlet cell velocity
                                     #vel_x in both ghost cells is vel_x in inlet cell
            gamma = dictionary["gamma"]
            vel_x = dictionary["vel_x"][0]
            Ma = dictionary["Ma"][0]
            P = BCs[0][1][0] * \
                (1 + 0.5 * (gamma[0] - 1) * Ma ** 2) ** (gamma[0] / (1 - gamma[0]))
            T = BCs[0][1][1] * (1 + 0.5 * (gamma[0] - 1) * Ma ** 2) ** -1
            gmodel = GasModel('ideal-air-gas-model.lua')
            inletState = GasState(gmodel)
            inletState.p = P
            inletState.T = T
            inletState.update_thermo_from_pT()
            inletState.update_sound_speed()
            inletState.update_trans_coeffs()
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass
                elif name == "p":
                        dictionary[name] = [P] * 2 \
                                + dictionary[name] 
                                
                elif name == "rho":
                        dictionary[name] = [inletState.rho] * 2 \
                                + dictionary[name] 
                               
                elif name == "vel_x":
                        dictionary[name] = [vel_x] * 2 \
                                + dictionary[name] 
                                
                elif name == "T":
                        dictionary[name] = [T] * 2 \
                                + dictionary[name] 
                                
                elif name == "a":
                        dictionary[name] = [inletState.a] * 2 \
                                + dictionary[name] 
                                
                elif name == "Ma":
                        dictionary[name] = \
                            [vel_x / inletState.a] * 2 \
                                + dictionary[name] 
                                
                elif name == "u":
                        dictionary[name] = [inletState.u] * 2 \
                                + dictionary[name] 
                               
                elif name == "h":
                        dictionary[name] = \
                            [gmodel.enthalpy(inletState)] * 2 \
                                + dictionary[name] 
                               
                elif name == "s":
                        dictionary[name] = [gmodel.entropy(inletState)] * 2 \
                                + dictionary[name] 
                                
                elif name == "Cv":
                        dictionary[name] = [gmodel.Cv(inletState)] * 2\
                                + dictionary[name] 
                
                elif name == "Cp":
                        dictionary[name] = [gmodel.Cp(inletState)] * 2\
                                + dictionary[name]
                
                elif name == "k":
                        dictionary[name] = [inletState.k] * 2\
                                + dictionary[name] 
                
                elif name == "mu":
                        dictionary[name] = [inletState.mu] * 2\
                                + dictionary[name]

                else: 
                        #For epsilon, mu, assume fixed value 
                        # (doesn't matter since these terms won't be used)
                        dictionary[name] = \
                                  [dictionary[name][0]]*2 + dictionary[name] 

        if BCs[1][0] == "OutFlow": #Zero gradient
            for name in entryNames:
                if name == "phi":
                    
                    point3 = 2 * dictionary["phi"][-1] - dictionary["phi"][-2]
                    point4 = 3 * dictionary["phi"][-1] - 2 * dictionary["phi"][-2]
                    dictionary["phi"] = np.hstack((dictionary["phi"], np.array([point3, point4])))
                    
                    psi3 = self.calculatePsi(phi = point3, dx = self.GEO["dx"][-2])
                    psi4 = self.calculatePsi(phi = point4, dx = self.GEO["dx"][-1])
                    dictionary["psi"] = np.hstack((dictionary["psi"], np.array([psi3, psi4])))
                elif name == "psi":
                    pass
                else:
                    dictionary[name] = np.hstack((dictionary[name], \
                                                dictionary[name][:-3:-1]))
            
        elif BCs[1][0] == "InFlow": #Fixed value
            gmodel = GasModel('ideal-air-gas-model.lua')
            outletState = GasState(gmodel)
            outletState.p = BCs[1][1][0]
            outletState.rho = BCs[1][1][1]
            outletState.update_thermo_from_rhop()
            outletState.update_sound_speed()
            outletState.update_trans_coeffs()
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass

                elif name == "p":
                        dictionary[name] =  dictionary[name] \
                                                        + [BCs[1][1][0]] * 2
                                
                elif name == "rho":
                        dictionary[name] = dictionary[name] + [BCs[1][1][1]] * 2
                               
                elif name == "vel_x":
                        dictionary[name] = dictionary[name] + [BCs[1][1][2]] * 2
                                
                elif name == "T":
                        dictionary[name] =  \
                                dictionary[name] + [outletState.T] * 2
                                
                elif name == "a":
                        dictionary[name] = \
                                dictionary[name] + [outletState.a] * 2 
                                
                elif name == "Ma":
                        dictionary[name] = dictionary[name] \
                                + [BCs[1][1][2]/outletState.a] * 2
                                
                elif name == "u":
                        dictionary[name] =  \
                            dictionary[name] + [outletState.u] * 2
                               
                elif name == "h":
                        dictionary[name] = dictionary[name] \
                                        + [gmodel.enthalpy(outletState)] * 2
                               
                elif name == "s":
                        dictionary[name] =  dictionary[name] \
                                        + [gmodel.entropy(outletState)] * 2
                                
                elif name == "Cv":
                        dictionary[name] = \
                                dictionary[name] + [gmodel.Cv(outletState)] * 2

                elif name == "Cp":
                        dictionary[name] = \
                                dictionary[name] + [gmodel.Cp(outletState)] * 2

                elif name == "k":
                        dictionary[name] = \
                                dictionary[name] + [outletState.k] * 2

                elif name == "mu":
                        dictionary[name] = \
                                dictionary[name] + [outletState.mu] * 2 

                else: 
                        #For epsilon, mu, assume fixed value 
                        # (doesn't matter since these terms won't be used)
                        dictionary[name] = dictionary[name] \
                            + [dictionary[name][-1]] * 2
            
        elif BCs[1][0] == "Wall":
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass

                elif name == "vel_x":
                    dictionary[name] = dictionary[name] \
                        + (-np.array(dictionary[name][-1:-3:-1])).tolist()

                else:
                    dictionary[name] = dictionary[name] \
                                        + dictionary[name][-1:-3:-1]

        elif BCs[1][0] == "FixedP":
            dictionary["p"] = dictionary["p"] + [BCs[1][1][0]] * 2
            dictionary["T"] = dictionary["T"] = dictionary["T"] \
                                        + dictionary["T"][-1:-3:-1]
            dictionary["vel_x"] = dictionary["vel_x"] = dictionary["vel_x"] \
                                        + dictionary["vel_x"][-1:-3:-1]                            
            gmodel = GasModel('ideal-air-gas-model.lua')
            GCellState = GasState(gmodel)
            
            for i in range(2):
                GCellState.p = dictionary["p"][-2 + GhostCell]
                GCellState.T = dictionary["T"][-2 + GhostCell]
                GCellState.update_thermo_from_pT()
                GCellState.update_sound_speed()
                GCellState.update_trans_coeffs()
                dictionary["rho"] = dictionary["rho"] + [GCellState.rho]
                dictionary["u"] = dictionary["u"] + [GCellState.u]
                dictionary["a"] = dictionary["a"] + [GCellState.a]
                dictionary["h"] = dictionary["h"] \
                                    + [gmodel.enthalpy(GCellState)]
                dictionary["s"] = dictionary["s"] + [gmodel.entropy(GCellState)]
                dictionary["Cv"] = dictionary["Cv"] + [gmodel.Cv(GCellState)]
                dictionary["Cp"] = dictionary["Cp"] + [gmodel.Cp(GCellState)]
                dictionary["k"] = dictionary["k"] + [GCellState.k]
                dictionary["mu"] = dictionary["mu"] + [GCellState.mu]
                dictionary["Ma"] = dictionary["Ma"] \
                    + [dictionary["vel_x"][-2 + GhostCell] / GCellState.a]
            
        elif BCs[1][0] == "FromStag":
            gamma = dictionary["gamma"]
            vel_x = dictionary["vel_x"][-1]
            Ma = dictionary["Ma"][-1]
            P = BCs[1][1][0] * \
                (1 + 0.5 * (gamma[-1] - 1) * Ma ** 2) ** (gamma[-1] / (1 - gamma[-1]))
            T = BCs[1][1][1] * (1 + 0.5 * (gamma[-1] - 1) * Ma ** 2) ** -1
            gmodel = GasModel('ideal-air-gas-model.lua')
            outletState = GasState(gmodel)
            outletState.p = P
            outletState.T = T
            outletState.update_thermo_from_pT()
            outletState.update_sound_speed()
            outletState.update_trans_coeffs()
            for name in entryNames:
                if name == "R" or name == "gamma":
                    pass
                elif name == "p":
                        dictionary[name] =  dictionary[name] + [P] * 2
                                
                elif name == "rho":
                        dictionary[name] = \
                            dictionary[name] + [outletState.rho] * 2
                               
                elif name == "vel_x":
                        dictionary[name] = dictionary[name] + [vel_x] * 2
                                
                elif name == "T":
                        dictionary[name] =  \
                                dictionary[name] + [T] * 2
                                
                elif name == "a":
                        dictionary[name] = \
                                dictionary[name] + [outletState.a] * 2 
                                
                elif name == "Ma":
                        dictionary[name] = \
                             dictionary[name] + [vel_x / outletState.a] * 2
                                
                elif name == "u":
                        dictionary[name] =  \
                            dictionary[name] + [outletState.u] * 2
                               
                elif name == "h":
                        dictionary[name] = dictionary[name] \
                                        + [gmodel.enthalpy(outletState)] * 2
                               
                elif name == "s":
                        dictionary[name] =  dictionary[name] \
                                        + [gmodel.entropy(outletState)] * 2
                                
                elif name == "Cv":
                        dictionary[name] = \
                                dictionary[name] + [gmodel.Cv(outletState)] * 2

                elif name == "Cp":
                        dictionary[name] = \
                                dictionary[name] + [gmodel.Cp(outletState)] * 2
                
                elif name == "k":
                        dictionary[name] = \
                                dictionary[name] + [outletState.k] * 2

                elif name == "mu":
                        dictionary[name] = \
                                dictionary[name] + [outletState.mu] * 2

                else: 
                        #For epsilon, mu, assume fixed value 
                        # (doesn't matter since these terms won't be used)
                        dictionary[name] = dictionary[name] \
                            + [dictionary[name][-1]] * 2
        return dictionary
    
    def calculateAverageProperties(self, dictionary, energyMethod):
        psi = dictionary["psi"]
        dictionary["rho"] = np.multiply(dictionary["rho_L"], psi) \
                                + np.multiply(dictionary["rho_G"], 1 - psi)
        dictionary["gamma_m"] = np.multiply(psi, dictionary["gamma_L"]) + np.multiply(1 - psi, dictionary["gamma_G"])
        #dictionary["gamma_m"] = (np.divide(psi, dictionary["gamma_L"] - 1) \
            #+ np.divide(1 - psi, dictionary["gamma_G"] - 1)) ** -1 + 1
        dictionary["a"] = np.multiply(psi, dictionary["a_L"]) \
            + np.multiply(1 - psi, dictionary["a_G"])
        dictionary["Ma"] = np.divide(dictionary["vel_x"], dictionary["a"])
        if energyMethod == 'Adapted':
            dictionary["u"] = np.divide(dictionary["p"], dictionary["gamma_m"] - 1) \
                + 0.5 * np.multiply(\
                            np.multiply(dictionary["vel_x"], dictionary["vel_x"]), \
                        dictionary["rho"])
        elif energyMethod == "Normal":
            dictionary["u"] = np.multiply(np.multiply(psi, dictionary["rho_L"]), dictionary["u_L"]) \
                + np.multiply(np.multiply(1 - psi, dictionary["rho_G"]), dictionary["u_G"]) \
                + 0.5 * np.multiply(\
                            np.multiply(dictionary["vel_x"], dictionary["vel_x"]), \
                        dictionary["rho"])
        
        #print(dictionary["u"])
        return dictionary
    
    def solveForConservedProps(self, a, liquid_a, vel_x, rho, rho_Liquid, p, u, Ma, A_c, phi, psi, fluxCalculator, reconstruction, sourceTerms, EOS, energyMethod):
        nCells = self.GEO["nCells"]
        dV = self.GEO["dV"]
        dt = self.dt
        dx = self.GEO["dx"]

        massFlux_e = np.full(nCells, np.nan)
        massFlux_w = np.full(nCells, np.nan)
        liquidMassFlux_e = np.full(nCells, np.nan)
        liquidMassFlux_w = np.full(nCells, np.nan)
        momentumFlux_w = np.full(nCells, np.nan)
        momentumFlux_e = np.full(nCells, np.nan)
        enthalpyFlux_e = np.full(nCells, np.nan)
        enthalpyFlux_w = np.full(nCells, np.nan)
        p_w = np.full(nCells, np.nan)
        p_e = np.full(nCells, np.nan)
        phi_w = np.full(nCells, np.nan)
        phi_e = np.full(nCells, np.nan)
        psi_w = np.full(nCells, np.nan)
        psi_e = np.full(nCells, np.nan)
        interface_velocity_w = np.full(nCells, np.nan)
        
        A_w = np.hstack((np.array([0.5 * (3*A_c[2] - A_c[3])]), 0.5 * (A_c[2:-3] + A_c[3:-2])))
        A_e = np.hstack((0.5 * (A_c[2:-3] + A_c[3:-2]), np.array([0.5 * (3*A_c[-3] - A_c[-4])])))
        density_new = np.full(nCells, np.nan)
        momentum_new = np.full(nCells, np.nan)
        totalInternalEnergy_new = np.full(nCells, np.nan)
        for cell in range(2, nCells + 2):
            ### Reconstruction
            if reconstruction == "Copy":
                    rho_L_w,        rho_R_w,      rho_L_e,      rho_R_e = \
                    rho[cell-1],    rho[cell],    rho[cell],    rho[cell + 1]

                    a_L_w,          a_R_w,        a_L_e,        a_R_e = \
                    a[cell-1],      a[cell],      a[cell],      a[cell + 1]

                    liquid_a_L_w,          liquid_a_R_w,        liquid_a_L_e,        liquid_a_R_e = \
                    psi[cell-1]*liquid_a[cell-1],      psi[cell]*liquid_a[cell],      psi[cell]*liquid_a[cell],      psi[cell+1]*liquid_a[cell + 1]

                    p_L_w,          p_R_w,        p_L_e,        p_R_e = \
                    p[cell-1],      p[cell],      p[cell],      p[cell + 1]

                    vel_x_L_w,      vel_x_R_w,    vel_x_L_e,    vel_x_R_e = \
                    vel_x[cell-1],  vel_x[cell],  vel_x[cell],  vel_x[cell + 1]

                    massFlux_L_w, massFlux_R_w, massFlux_L_e, massFlux_R_e = \
                                vel_x_L_w * rho_L_w, vel_x_R_w * rho_R_w, \
                                vel_x_L_e * rho_L_e, vel_x_R_e * rho_R_e
                    
                    H_L_w, H_R_w, H_L_e, H_R_e = \
                                        (u[cell-1] + p[cell-1])/rho[cell-1], \
                                        (u[cell] + p[cell])/rho[cell], \
                                        (u[cell] + p[cell])/rho[cell], \
                                        (u[cell+1] + p[cell+1])/rho[cell+1]
                    
                    Ma_L_w, Ma_R_w, Ma_L_e, Ma_R_e = \
                    Ma[cell-1],  Ma[cell],  Ma[cell],  Ma[cell + 1]

                    liquid_Ma_L_w, liquid_Ma_R_w, liquid_Ma_L_e, liquid_Ma_R_e = \
                    vel_x[cell-1]/liquid_a[cell-1],  vel_x[cell]/liquid_a[cell],  vel_x[cell]/liquid_a[cell],  vel_x[cell + 1]/liquid_a[cell+1]

                    rho_Liquid_L_w,         rho_Liquid_R_w, \
                    rho_Liquid_L_e,         rho_Liquid_R_e = \
                    rho_Liquid[cell-1],     rho_Liquid[cell], \
                    rho_Liquid[cell],       rho_Liquid[cell + 1]
            if fluxCalculator == "AUSMPlusUpBASIC":
                a_interface_w = 0.5 * (a_L_w + a_R_w)
                a_interface_e = 0.5 * (a_L_e + a_R_e)
                liquid_a_interface_w = 0.5 * (liquid_a_L_w + liquid_a_R_w)
                liquid_a_interface_e = 0.5 * (liquid_a_L_e + liquid_a_R_e)
                
                """
                a_interface_w       = AUSMPlusUpBASIC().interfaceSoundSpeed(\
                    a_crit          = max(abs(vel_x_L_w), abs(vel_x_R_w)), \
                    vel_x_L         = vel_x_L_w,    vel_x_R     = vel_x_R_w)
                a_interface_e       = AUSMPlusUpBASIC().interfaceSoundSpeed(\
                    a_crit          = max(abs(vel_x_L_e), abs(vel_x_R_e)), \
                    vel_x_L         = vel_x_L_e,    vel_x_R     = vel_x_R_e)
                """
                Ma_interface_w      = AUSMPlusUpBASIC().interfaceMachNumber(\
                    order = 4, \
                    rho_L           = rho_L_w,      rho_R       = rho_R_w,\
                    a_L             = a_L_w,        a_R         = a_R_w, \
                    p_L             = p_L_w,        p_R         = p_R_w, \
                    vel_x_L         = vel_x_L_w,    vel_x_R     = vel_x_R_w, \
                    Ma_L            = vel_x_L_w/a_interface_w,       Ma_R        = vel_x_R_w/a_interface_w)
                Ma_interface_e      = AUSMPlusUpBASIC().interfaceMachNumber(\
                    order = 4, \
                    rho_L           = rho_L_e,      rho_R       = rho_R_e, \
                    a_L             = a_L_e,        a_R         = a_R_e, \
                    p_L             = p_L_e,        p_R         = p_R_e, \
                    vel_x_L         = vel_x_L_e,    vel_x_R     = vel_x_R_e, \
                    Ma_L            = vel_x_L_e/a_interface_e,       Ma_R        = vel_x_R_e/a_interface_e)
                if abs(liquid_a_interface_w) < 1e-9:
                    liquid_Ma_interface_w = 0.0
                else:
                    liquid_Ma_interface_w      = AUSMPlusUpBASIC().interfaceMachNumber(\
                    order = 4, \
                    rho_L           = rho_Liquid_L_w,      rho_R       = rho_Liquid_R_w,\
                    a_L             = liquid_a_L_w,        a_R         = liquid_a_R_w, \
                    p_L             = p_L_w,        p_R         = p_R_w, \
                    vel_x_L         = vel_x_L_w,    vel_x_R     = vel_x_R_w, \
                    Ma_L            = vel_x_L_w/liquid_a_interface_w,       Ma_R        = vel_x_R_w/liquid_a_interface_w)
                if abs(liquid_a_interface_e) < 1e-9:
                    liquid_Ma_interface_e = 0.0
                else:
                    liquid_Ma_interface_e     = AUSMPlusUpBASIC().interfaceMachNumber(\
                    order = 4, \
                    rho_L           = rho_Liquid_L_e,      rho_R       = rho_Liquid_R_e, \
                    a_L             = liquid_a_L_e,        a_R         = liquid_a_R_e, \
                    p_L             = p_L_e,        p_R         = p_R_e, \
                    vel_x_L         = vel_x_L_e,    vel_x_R     = vel_x_R_e, \
                    Ma_L            = vel_x_L_e/liquid_a_interface_e,       Ma_R        = vel_x_R_e/liquid_a_interface_e)
                massFlux_w[cell - 2] = AUSMPlusUpBASIC().interfaceMassFlux(\
                    a_interface     = a_interface_w, \
                    Ma_interface    = Ma_interface_w, \
                    rho_L           = rho_L_w,      rho_R       = rho_R_w)
                massFlux_e[cell - 2] = AUSMPlusUpBASIC().interfaceMassFlux(\
                    a_interface     = a_interface_e, \
                    Ma_interface    = Ma_interface_e, \
                    rho_L           = rho_L_e,      rho_R       = rho_R_e)
                liquidMassFlux_w[cell-2] = AUSMPlusUpBASIC().interfaceMassFlux(\
                    a_interface     = liquid_a_interface_w, \
                    Ma_interface    = liquid_Ma_interface_w, \
                    rho_L           = rho_Liquid_L_w, \
                    rho_R           = rho_Liquid_R_w)
                liquidMassFlux_e[cell-2] = AUSMPlusUpBASIC().interfaceMassFlux(\
                    a_interface     = liquid_a_interface_e, \
                    Ma_interface    = liquid_Ma_interface_e, \
                    rho_L           = rho_Liquid_L_e, \
                    rho_R           = rho_Liquid_R_e)
                momentumFlux_w[cell - 2] = AUSMPlusUpBASIC().interfaceMomentumFlux(\
                    a_interface     = a_interface_w, \
                    Ma_interface    = Ma_interface_w, \
                    rho_L           = rho_L_w,      rho_R       = rho_R_w, \
                    vel_x_L         = vel_x_L_w,    vel_x_R     = vel_x_R_w)
                momentumFlux_e[cell - 2] = AUSMPlusUpBASIC().interfaceMomentumFlux(\
                    a_interface     = a_interface_e, \
                    Ma_interface    = Ma_interface_e, \
                    rho_L           = rho_L_e,      rho_R       = rho_R_e, \
                    vel_x_L         = vel_x_L_e,    vel_x_R     = vel_x_R_e)
                
                enthalpyFlux_w[cell - 2] = AUSMPlusUpBASIC().interfaceEnthalpyFlux(\
                    a_interface     = a_interface_w, \
                    Ma_interface    = Ma_interface_w, \
                    rho_L           = rho_L_w,      rho_R       = rho_R_w, \
                    H_L             = H_L_w,        H_R         = H_R_w)
                enthalpyFlux_e[cell - 2] = AUSMPlusUpBASIC().interfaceEnthalpyFlux(\
                    a_interface     = a_interface_e, \
                    Ma_interface    = Ma_interface_e, \
                    rho_L           = rho_L_e,      rho_R       = rho_R_e, \
                    H_L             = H_L_e,        H_R         = H_R_e)
                
                p_w[cell - 2]       = AUSMPlusUpBASIC().interfacePressure(\
                    rho_L           = rho_L_w,      rho_R       = rho_R_w, \
                    a_L             = a_L_w,        a_R         = a_R_w, \
                    p_L             = p_L_w,        p_R         = p_R_w, \
                    vel_x_L         = vel_x_L_w,    vel_x_R     = vel_x_R_w, \
                    Ma_L            = vel_x_L_w/a_interface_w,       Ma_R        = vel_x_R_w/a_interface_w, \
                    a_crit          = max(abs(vel_x_L_w), abs(vel_x_R_w)))
                p_e[cell - 2]       = AUSMPlusUpBASIC().interfacePressure(\
                    rho_L           = rho_L_e,      rho_R       = rho_R_e, \
                    a_L             = a_L_e,        a_R         = a_R_e, \
                    p_L             = p_L_e,        p_R         = p_R_e, \
                    vel_x_L         = vel_x_L_e,    vel_x_R     = vel_x_R_e, \
                    Ma_L            = vel_x_L_e/a_interface_e,       Ma_R        = vel_x_R_e/a_interface_e, \
                    a_crit          = max(abs(vel_x_L_e), abs(vel_x_R_e)))
                phi_w[cell-2] = 0.5 * (phi[cell - 1] + phi[cell])
                phi_e[cell-2] = 0.5 * (phi[cell] + phi[cell + 1])   
                psi_w[cell-2] = self.calculatePsi(phi = phi_w[cell-2], dx = dx[cell])
                psi_e[cell-2] = self.calculatePsi(phi = phi_e[cell-2], dx = dx[cell])
                interface_velocity_w[cell-2] = massFlux_w[cell - 2] / (0.5 * (rho[cell-1] +   rho[cell]))
        
        ### Mass
        rho_old_no_BCs = rho[2:-2]
        
        TotalMassFlux_w = np.multiply(A_w, massFlux_w)
        TotalMassFlux_e = np.multiply(A_e, massFlux_e)
        """
        print(TotalMassFlux_w[125])
        print(TotalMassFlux_e[125])
        print(TotalMassFlux_w[126])
        print(TotalMassFlux_e[126])
        """
        rho_new_no_BCs = rho_old_no_BCs + dt * np.divide(TotalMassFlux_w - TotalMassFlux_e, dV[2:-2])
        
        
        ### Momentum
        TotalxMomentumFlux_w = np.multiply(A_w, momentumFlux_w)
        TotalxMomentumFlux_e = np.multiply(A_e, momentumFlux_e)
        PressureForceDifferential = np.multiply(A_c[2:-2], p_e - p_w)
        rhovel_x_old_no_BCs = np.multiply(rho[2:-2], vel_x[2:-2])
        rhovel_x_new_no_BCs = rhovel_x_old_no_BCs + dt * np.divide(TotalxMomentumFlux_w - TotalxMomentumFlux_e - PressureForceDifferential, dV[2:-2])

        
        ### Liquid mass
        
        rho_L_psi_old_no_BCs = np.multiply(rho_Liquid[2:-2], psi[2:-2])
        TotalLiquidMassFlux_w = np.multiply(A_w, liquidMassFlux_w)
        TotalLiquidMassFlux_e = np.multiply(A_e, liquidMassFlux_e)
        """
        print(TotalLiquidMassFlux_w[125])
        print(TotalLiquidMassFlux_e[125])
        print(TotalLiquidMassFlux_w[126])
        print(TotalLiquidMassFlux_e[126])
        """
        rho_L_psi_no_BCs = rho_L_psi_old_no_BCs + dt * np.divide(np.multiply(TotalLiquidMassFlux_w, psi_w) - np.multiply(TotalLiquidMassFlux_e, psi_e), dV[2:-2])
        for ind, value in enumerate(rho_L_psi_no_BCs):
            if np.abs(value) < 1e-6:
                rho_L_psi_no_BCs[ind] = 0.0

        ### Signed distance
        phi_old_no_BCs = phi[2:-2]
        phi_W_old = phi[1:-3]
        phi_P_old = phi[2:-2]
        dx_WP = dx[2:-2]
        
        phi_no_BCs = phi_old_no_BCs + dt * np.divide(np.multiply((interface_velocity_w), phi_W_old - phi_P_old), dx_WP)
        new_phi1 = 3 * phi_no_BCs[0] - 2 * phi_no_BCs[1]
        new_phi2 = 2 * phi_no_BCs[0] - phi_no_BCs[1]
        new_phi3 = 2 * phi_no_BCs[-1] - phi_no_BCs[-2]
        new_phi4 = 3 * phi_no_BCs[-1] - 2 * phi_no_BCs[-2]
        print(interface_velocity_w.tolist())
        phi_with_BCs = np.hstack((np.array([new_phi1, new_phi2]), phi_no_BCs, np.array([new_phi3, new_phi4])))
        
        steadyStatePhi = self.SteadyStatePhi(phi = phi_with_BCs)
        steadyStatePsi = np.full(nCells + 4, np.nan)
        for cell in range(nCells + 4):
            steadyStatePsi[cell] = self.calculatePsi(phi = steadyStatePhi[cell], dx = dx[cell])

        ### Energy
        if energyMethod == "Normal":
            # Original Method
            rhoU_old_no_BCs = u[2:-2]
            TotalEnthalpyFlux_w = np.multiply(A_w, enthalpyFlux_w)
            TotalEnthalpyFlux_e = np.multiply(A_e, enthalpyFlux_e)
            rhoU_new_no_BCs = rhoU_old_no_BCs + dt * np.divide(TotalEnthalpyFlux_w - TotalEnthalpyFlux_e, dV[2:-2])
        
        elif energyMethod == "Adapted":
            # Their method
            p_L_inf = coefficients["water"]["p_inf"]
            gamma_L = self.currentProps["gamma_L"][2:-2]
            gamma_L_dummy = copy.deepcopy(gamma_L)
            for cell in range(2, nCells + 2):
                if steadyStatePsi[cell] != 0.0 or steadyStatePsi[cell] != 1.0 or steadyStatePsi[cell-1] != 0.0 or steadyStatePsi[cell-1] != 1.0:
                    gamma_L_dummy[cell-2] = 0.0
            dx_PW = 0.5 * (dx[1:-3] + dx[2:-2])
            e_prime_old_no_BCs = u[2:-2]
            dpsidt = (steadyStatePsi[2:-2] - psi[2:-2])/dt
            dpsivel_xdx = np.divide(np.multiply(psi[2:-2], vel_x[2:-2]) - np.multiply(psi[1:-3], vel_x[1:-3]), dx_PW)
            TotalEnthalpyFlux_w = np.multiply(A_w, enthalpyFlux_w)
            TotalEnthalpyFlux_e = np.multiply(A_e, enthalpyFlux_e)
            
            e_prime_no_BCs = e_prime_old_no_BCs + dt * np.divide(TotalEnthalpyFlux_w - TotalEnthalpyFlux_e, dV[2:-2]) \
                - dt * p_L_inf * np.multiply(np.divide(gamma_L_dummy, gamma_L - 1), dpsidt + dpsivel_xdx)
            
        #self.newProps["rho_L_psi"] = rho_L_psi_no_BCs
        #print(TotalMassFlux_w.tolist())
        """
        print('test')
        print(rho_old_no_BCs.tolist())
        print(rho_new_no_BCs.tolist())
        #print(rhovel_x_old_no_BCs.tolist())
        #print(rhovel_x_new_no_BCs.tolist())
        print(rho_L_psi_old_no_BCs.tolist())
        print(rho_L_psi_no_BCs.tolist())
        #print(phi_no_BCs.tolist())
        #print(np.shape(steadyStatePhi))
        #print(np.shape(steadyStatePsi))
        #print(steadyStatePhi.tolist())
        print(psi[2:-2].tolist())
        print(steadyStatePsi[2:-2].tolist())
        
        if energyMethod == "Normal":
            print(rhoU_old_no_BCs.tolist())
            print(rhoU_new_no_BCs.tolist())
        elif energyMethod == "Adapted":
            print(e_prime_old_no_BCs.tolist())
            print(e_prime_no_BCs.tolist())
        """
        ### Reconstruct variables
        if energyMethod == "Normal":

            self.reconstructPrimaryVariables(energyMethod = energyMethod, \
                rho = rho_new_no_BCs, \
                rho_LPsi = rho_L_psi_no_BCs, \
                xMomentum = rhovel_x_new_no_BCs, \
                phi = steadyStatePhi[2:-2], \
                energy = rhoU_new_no_BCs)
        elif energyMethod == "Adapted":
            self.reconstructPrimaryVariables(energyMethod = energyMethod, \
                rho = rho_new_no_BCs, \
                rho_LPsi = rho_L_psi_no_BCs, \
                xMomentum = rhovel_x_new_no_BCs, \
                phi = steadyStatePhi[2:-2], \
                energy = e_prime_no_BCs)
        
    def SteadyStatePhi(self, phi):
        nCells = self.GEO["nCells"]
        dx = self.GEO["dx"]
        new_phi = phi.copy()
        closestIndexToInterface = np.argmin(np.abs(new_phi))
        closestCellToInterface = new_phi[closestIndexToInterface]
        #print(closestCellToInterface)
        for cell in range(2, nCells + 2):
                if cell < closestIndexToInterface:
                        new_phi[cell] = closestCellToInterface - np.sum(dx[cell:closestIndexToInterface])
                elif cell > closestIndexToInterface:
                        new_phi[cell] = closestCellToInterface + np.sum(dx[closestIndexToInterface:cell])
                else:
                        new_phi[cell] = closestCellToInterface
        new_phi[0] = 3 * new_phi[2] - 2 * new_phi[3]
        new_phi[1] = 2 * new_phi[2] - new_phi[3]
        new_phi[-2] = 2 * new_phi[-3] - new_phi[-4]
        new_phi[-1] = 3 * new_phi[-3] - 2 * new_phi[-4]
        return new_phi

    def calculatePsi(self, phi, dx):
        alpha = 2.5 * dx
        if phi > alpha:
            return 1
        elif phi < - alpha:
            return 0
        else:
            return 0.5 * (1 + phi/alpha + 1 / m.pi * m.sin(m.pi * phi / alpha))
    
    def reconstructPrimaryVariables(self, energyMethod, rho, rho_LPsi, \
                                            xMomentum, phi, energy):
        nCells = self.GEO["nCells"]
        dx = self.GEO["dx"]
        rho_L = np.full(nCells, np.nan)
        rho_G = np.full(nCells, np.nan)
        u_L = np.full(nCells, np.nan)
        u_G = np.full(nCells, np.nan)
        vel_x = np.full(nCells, np.nan)
        psi = np.full(nCells, np.nan)
        p = np.full(nCells, np.nan)

        b_L = coefficients["water"]["b"]
        q_L = coefficients["water"]["q"]
        p_inf_L = coefficients["water"]["p_inf"]
        
        if energyMethod == "Normal":
            for cell in range(nCells):
                #print(cell)
                psi_val = self.calculatePsi(phi = phi[cell], dx = dx[cell + 2])
                if abs(psi_val) < 1e-6:
                    psi[cell] = 0.0
                elif abs(psi_val - 1.0) < 1e-6:
                    psi[cell] = 1.0
                else:
                    psi[cell] = psi_val
                gamma_L = self.currentProps["gamma_L"][cell]
                gamma_G = self.currentProps["gamma_G"][cell]
                #print(psi[cell], rho_LPsi[cell])
                if psi[cell] == 0.0:
                    ### rho_L
                    rho_L[cell] = self.currentProps["rho_L"][cell + 2]
                    ### rho_G
                    rho_G[cell] = (rho[cell] - rho_LPsi[cell]) / (1-psi[cell])
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### u_L
                    u_L[cell] = ((gamma_L - 1) * q_L / (rho_G[cell] * (gamma_G - 1) * (1 / rho_L[cell] - b_L)) \
                        + gamma_L * p_inf_L / (rho_G[cell] * (gamma_G - 1)) \
                            + (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) / (rho_G[cell] * (1 - psi[cell]))) \
                                / ((gamma_L - 1) / (rho_G[cell] * (gamma_G - 1) * (1 / rho_L[cell] - b_L)) \
                                    + rho_LPsi[cell] / (rho_G[cell] * (1 - psi[cell])))
                    ### u_G
                    u_G[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell] - rho_LPsi[cell] * u_L[cell]) \
                        / (rho_G[cell] * (1-psi[cell]))
                    ### P
                    p[cell] = (gamma_G - 1) * rho_G[cell] * u_G[cell]
                    continue
                elif psi[cell] == 1.0:
                    ### rho_L
                    rho_L[cell] = rho_LPsi[cell] / psi[cell]
                    ### rho_G
                    rho_G[cell] = self.currentProps["rho_G"][cell + 2]
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### u_G
                    u_G[cell] = (((energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) / rho_LPsi[cell]) \
                        - q_L - (1 / rho_L[cell] - b_L) * gamma_L * p_inf_L / (gamma_L - 1)) \
                            / ((1 - psi[cell]) * rho_G[cell] / rho_LPsi[cell] \
                                + rho_G[cell] * (1 / rho_L[cell] - b_L) * (gamma_G - 1) / (gamma_L - 1))
                    ### u_L
                    u_L[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell] - (1 - psi[cell]) * rho_G[cell] * u_G[cell]) \
                        / rho_LPsi[cell]
                    ### P
                    p[cell] = (gamma_L - 1) * (u_L[cell] - q_L) / (1 / rho_L[cell] - b_L) - gamma_L * p_inf_L
                    continue
                else:
                    ### rho_L
                    rho_L[cell] = rho_LPsi[cell] / psi[cell]
                    #print(rho_L[cell])
                    ### rho_G
                    rho_G[cell] = (rho[cell] - rho_LPsi[cell]) / (1-psi[cell])
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### u_G
                    u_G[cell] = (((energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) / rho_LPsi[cell]) \
                        - q_L - (1 / rho_L[cell] - b_L) * gamma_L * p_inf_L / (gamma_L - 1)) \
                            / ((1 - psi[cell]) * rho_G[cell] / rho_LPsi[cell] \
                                + rho_G[cell] * (1 / rho_L[cell] - b_L) * (gamma_G - 1) / (gamma_L - 1))
                    ### u_L
                    u_L[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell] - (1 - psi[cell]) * rho_G[cell] * u_G[cell]) \
                        / rho_LPsi[cell]
                    ### P
                    p[cell] = (gamma_G - 1) * rho_G[cell] * u_G[cell]
                    #print(psi[cell], rho_L[cell], rho_G[cell], vel_x[cell], u_G[cell], u_L[cell], p[cell])
                    continue
                
        elif energyMethod == "Adapted":
            for cell in range(nCells):
                #print(cell)
                psi_val = self.calculatePsi(phi = phi[cell], dx = dx[cell + 2])
                if abs(psi_val) < 1e-6:
                    psi[cell] = 0.0
                elif abs(psi_val - 1.0) < 1e-6:
                    psi[cell] = 1.0
                else:
                    psi[cell] = psi_val
                gamma_L = self.currentProps["gamma_L"][cell]
                gamma_G = self.currentProps["gamma_G"][cell]
                gamma_m = psi[cell] * gamma_L + (1-psi[cell]) * gamma_G
                #gamma_m = (psi[cell] / (gamma_L - 1) + (1-psi[cell])/(gamma_G - 1)) ** -1 + 1
                #print(psi[cell], rho_LPsi[cell])
                if psi[cell] == 0.0:
                    ### rho_L
                    rho_L[cell] = self.currentProps["rho_L"][cell + 2]
                    ### rho_G
                    rho_G[cell] = (rho[cell] - rho_LPsi[cell]) / (1-psi[cell])
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### p
                    p[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) * (gamma_m - 1)
                    ### u_L
                    u_L[cell] = (p[cell] + gamma_L * p_inf_L) * (1 / rho_L[cell] - b_L) / (gamma_L - 1) + q_L
                    ### u_G
                    u_G[cell] = p[cell] / (rho_G[cell] * (gamma_G - 1))
                elif psi[cell] == 1.0:
                    ### rho_L
                    rho_L[cell] = rho_LPsi[cell] / psi[cell]
                    ### rho_G
                    rho_G[cell] = self.currentProps["rho_G"][cell + 2]
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### p
                    p[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) * (gamma_m - 1)
                    ### u_L
                    u_L[cell] = (p[cell] + gamma_L * p_inf_L) * (1 / rho_L[cell] - b_L) / (gamma_L - 1) + q_L
                    ### u_G
                    u_G[cell] = p[cell] / (rho_G[cell] * (gamma_G - 1))
                else:
                    ### rho_L
                    rho_L[cell] = rho_LPsi[cell] / psi[cell]
                    ### rho_G
                    rho_G[cell] = (rho[cell] - rho_LPsi[cell]) / (1-psi[cell])
                    ### vel_x
                    vel_x[cell] = xMomentum[cell] / rho[cell]
                    ### p
                    p[cell] = (energy[cell] - 0.5 * xMomentum[cell] * vel_x[cell]) * (gamma_m - 1)
                    ### u_L
                    u_L[cell] = (p[cell] + gamma_L * p_inf_L) * (1 / rho_L[cell] - b_L) / (gamma_L - 1) + q_L
                    ### u_G
                    u_G[cell] = p[cell] / (rho_G[cell] * (gamma_G - 1))
        #print("next")
        #print(self.currentProps["psi"][128])   
        #print(psi[127])
        #print(self.currentProps["phi"][127])
        #print(phi[127])
        #print(rho_LPsi[128])
        #print(rho[128])
        #print(self.currentProps["rho_G"][2:-2].tolist())
        #print(rho_G.tolist())
        #print(self.currentProps["rho_L"][2:-2].tolist())
        #print(rho_L.tolist())
        #print(self.currentProps["p"][127])
        #print(p[127])
        #print(self.currentProps["u_L"][127])
        #print(u_L[127])
        #print(self.currentProps["u_G"][127])
        #print(u_G[127])
        #print(self.currentProps["vel_x"][127])
        #print(vel_x[127])
        
        self.newProps["psi"] = psi
        self.newProps["phi"] = phi
        self.newProps["rho_G"] = rho_G
        self.newProps["rho_L"] = rho_L
        self.newProps["p"] = p
        self.newProps["u_L"] = u_L
        self.newProps["u_G"] = u_G
        self.newProps["vel_x"] = vel_x   
            
    def updateNewProps(self, dictionary):
        nCells = self.GEO["nCells"]

        T_G = np.full(nCells, np.nan)
        T_L = np.full(nCells, np.nan)
        h_G = np.full(nCells, np.nan)
        h_L = np.full(nCells, np.nan)
        Cp_G = np.full(nCells, np.nan)
        Cp_L = np.full(nCells, np.nan)
        Cv_G = np.full(nCells, np.nan)
        Cv_L = np.full(nCells, np.nan)
        R_G = np.full(nCells, np.nan)
        gamma_G = np.full(nCells, np.nan)
        gamma_L = np.full(nCells, np.nan)
        a_G = np.full(nCells, np.nan)
        a_L = np.full(nCells, np.nan)


        p = dictionary["p"]
        #print(p)
        rho_L = dictionary["rho_L"]
        rho_G = dictionary["rho_G"]
                
        for cell in range(nCells):
            ### Fill gas props
            gmodel = GasModel("ideal-air-gas-model.lua")
            cellState_G = GasState(gmodel=gmodel)
            #print(cell)
            #print(p[cell], rho_G[cell])
            #print(self.currentProps["p"][cell+2], self.currentProps["rho_G"][cell+2])
            cellState_G.p = p[cell]
            cellState_G.rho = rho_G[cell]
            cellState_G.update_thermo_from_rhop()
            cellState_G.update_sound_speed()
            T_G[cell] = cellState_G.T
            h_G[cell] = gmodel.enthalpy(cellState_G)
            Cp_G[cell] = cellState_G.Cp
            Cv_G[cell] = cellState_G.Cv
            R_G[cell] = gmodel.R(cellState_G)
            gamma_G[cell] = gmodel.gamma(cellState_G)
            a_G[cell] = cellState_G.a

            ### Fill liquid props
            coeffs = coefficients["water"]
            liquidProps = {}
            liquidProps["rho"] = rho_L[cell]
            liquidProps["p"] = p[cell]
            liquidProps = LiquidEOS.update_from_rhoP(state = liquidProps, coefficients = coeffs)
            T_L[cell] = liquidProps["T"]
            h_L[cell] = liquidProps["h"]
            Cp_L[cell] = liquidProps["Cp"]
            Cv_L[cell] = liquidProps["Cv"]
            gamma_L[cell] = liquidProps["gamma"]
            a_L[cell] = liquidProps["a"]

        dictionary["T_G"] = T_G
        dictionary["T_L"] = T_L
        dictionary["h_G"] = h_G
        dictionary["h_L"] = h_L
        dictionary["Cp_G"] = Cp_G
        dictionary["Cp_L"] = Cp_L
        dictionary["Cv_G"] = Cv_G
        dictionary["Cv_L"] = Cv_L
        dictionary["R_G"] = R_G
        dictionary["gamma_G"] = gamma_G
        dictionary["gamma_L"] = gamma_L
        dictionary["a_G"] = a_G
        dictionary["a_L"] = a_L
        return dictionary