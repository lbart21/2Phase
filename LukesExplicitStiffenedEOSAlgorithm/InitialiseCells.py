import numpy as np
import math as m
from eilmer.gas import GasModel, GasState
from LukesExplicitStiffenedEOSAlgorithm.EOSElementCoefficients import coefficients
import LukesExplicitStiffenedEOSAlgorithm.LiquidEOS as LiquidEOS
class InitialiseCells():
    def __init__(self, Geometry, pBounds, vel_xBounds, rhoBounds, interfaceIndex) -> None:
        """
        Geometry = [str, [data]]
            str == "StraightPipe" \
                -> data = D, L, nCells
            str == "LinearNozzle" \
                -> data = D_left, D_right, L, nCells
            str == "SinusoidalNozzle" \
                -> data = D_left, D_right, L, nCells
            str == "QuadraticNozzle" \
                -> data = D_left, D_right, L, nCells
            str == "LinearDoubleNozzle" \
                -> data = D_left, D_throat, D_right, \
                    L_1, L_2, nCells_1, nCells_2
            str == "SinusoidalDoubleNozzle" \
                -> data = D_left, D_throat, D_right, \
                    L_1, L_2, nCells_1, nCells_2
            str == "QuadraticDoubleNozzle" \
                -> data = D_left, D_throat, D_right, \
                    L_1, L_2, nCells_1, nCells_2
        """
        self.GEO = {
        }
        self.props = {
        }
        self.fillGEO(Geometry = Geometry)
        self.fillBasicProps(pBounds              = pBounds, \
                        vel_xBounds         = vel_xBounds, \
                        rhoBounds           = rhoBounds, \
                        interfaceIndex   = interfaceIndex)
        self.fillRestOfPropsFromRHOPVEL_X()

    def fillGEO(self, Geometry):
        if Geometry[0] == "StraightPipe":
            D, L, nCells = Geometry[1]
            self.GEO["nCells"] = nCells
            self.GEO["dx"] = np.full(nCells, L / nCells)
            self.GEO["A_c"] = np.full(nCells, 0.25 * m.pi * D ** 2)
            self.GEO["dV"] = np.full(nCells, 0.25 * m.pi * D ** 2 * L / nCells)
            self.GEO["A_s"] = np.full(nCells, m.pi * D * L / nCells)
            self.GEO["pos_x"] = np.linspace(\
                            0.5 * L / nCells, L * (1 - 0.5/nCells), nCells)

        elif Geometry[0] == "LinearNozzle":
            D_left, D_right, L, nCells = Geometry[1]
            self.GEO["nCells"] = nCells
            self.GEO["dx"] = np.full(nCells, L / nCells)
            self.GEO["pos_x"] = np.linspace(\
                            0.5 * L / nCells, L * (1 - 0.5/nCells), nCells)
            self.GEO["A_c"] = m.pi * (0.5 * D_left + 0.5 * (D_right - D_left) / L * self.GEO["pos_x"]) ** 2
            self.GEO["dV"] = np.multiply(self.GEO["A_c"], self.GEO["dx"])
            self.GEO["A_s"] = 2 * m.pi * np.multiply(self.GEO["dx"], 0.5 * D_left + 0.5 * (D_right - D_left) / L * self.GEO["pos_x"])

        elif Geometry[0] == "SinusoidalNozzle":
            D_left, D_right, L, nCells = Geometry[1]
            self.GEO["nCells"] = nCells
            self.GEO["dx"] = np.full(nCells, L / nCells)
            self.GEO["pos_x"] = np.linspace(\
                            0.5 * L / nCells, L * (1 - 0.5/nCells), nCells)
            self.GEO["A_c"] = m.pi * (0.25 * (D_left + D_right) + 0.25 * (D_left - D_right) * np.cos(self.GEO["pos_x"] * m.pi / L)) ** 2
            self.GEO["dV"] = np.multiply(self.GEO["A_c"], self.GEO["dx"])
            self.GEO["A_s"] = 2 * m.pi * np.multiply(self.GEO["dx"], 0.25 * (D_left + D_right) + 0.25 * (D_left - D_right) * np.cos(self.GEO["pos_x"] * m.pi / L))

        elif Geometry[0] == "QuadraticNozzle":
            D_left, D_right, L, nCells = Geometry[1]
            self.GEO["nCells"] = nCells
            self.GEO["dx"] = np.full(nCells, L / nCells)
            self.GEO["pos_x"] = np.linspace(\
                            0.5 * L / nCells, L * (1 - 0.5/nCells), nCells)
            self.GEO["A_c"] = m.pi * (0.5 * D_right + 0.5 * (D_left - D_right) * (1 - self.GEO["pos_x"] / L) ** 2) ** 2
            self.GEO["dV"] = np.multiply(self.GEO["A_c"], self.GEO["dx"])
            self.GEO["A_s"] = 2 * m.pi * np.multiply(self.GEO["dx"], 0.5 * D_right + 0.5 * (D_left - D_right) * (1 - self.GEO["pos_x"] / L) ** 2)

        elif Geometry[0] == "LinearDoubleNozzle":
            D_left, D_throat, D_right, L_1, L_2, nCells_1, nCells_2 = Geometry[1]
            
            dx_1 = np.full(nCells_1, L_1 / nCells_1)
            pos_x_1 = np.linspace(\
                            0.5 * L_1 / nCells_1, L_1 * (1 - 0.5/nCells_1), nCells_1)
            A_c_1 = m.pi * (0.5 * D_left + 0.5 * (D_throat - D_left) / L_1 * pos_x_1) ** 2
            dV_1 = np.multiply(dx_1, A_c_1)
            A_s_1 = 2 * m.pi * np.multiply(dx_1, 0.5 * D_left + 0.5 * (D_throat - D_left) / L_1 * pos_x_1)

            dx_2 = np.full(nCells_2, L_2 / nCells_2)
            pos_x_2 = np.linspace(\
                            0.5 * L_2 / nCells_2, L_2 * (1 - 0.5/nCells_2), nCells_2)
            A_c_2 = m.pi * (0.5 * D_throat + 0.5 * (D_right - D_throat) / L_2 * pos_x_2) ** 2
            dV_2 = np.multiply(dx_2, A_c_2)
            A_s_2 = 2 * m.pi * np.multiply(dx_2, 0.5 * D_throat + 0.5 * (D_right - D_throat) / L_2 * pos_x_2)

            self.GEO["nCells"] = nCells_1 + nCells_2
            self.GEO["dx"] = np.hstack((dx_1, dx_2))
            self.GEO["pos_x"] = np.hstack((pos_x_1, L_1 + pos_x_2))
            self.GEO["A_c"] = np.hstack((A_c_1, A_c_2))
            self.GEO["dV"] = np.hstack((dV_1, dV_2))
            self.GEO["A_s"] = np.hstack((A_s_1, A_s_2))

        elif Geometry[0] == "SinusoidalDoubleNozzle":
            D_left, D_throat, D_right, L_1, L_2, nCells_1, nCells_2 = Geometry[1]
            
            dx_1 = np.full(nCells_1, L_1 / nCells_1)
            pos_x_1 = np.linspace(\
                            0.5 * L_1 / nCells_1, L_1 * (1 - 0.5/nCells_1), nCells_1)
            A_c_1 = m.pi * (0.25 * (D_left + D_throat) + 0.25 * (D_left - D_throat) * np.cos(pos_x_1 * m.pi / L_1)) ** 2
            dV_1 = np.multiply(dx_1, A_c_1)
            A_s_1 = 2 * m.pi * np.multiply(dx_1, 0.25 * (D_left + D_throat) + 0.25 * (D_left - D_throat) * np.cos(pos_x_1 * m.pi / L_1))

            dx_2 = np.full(nCells_2, L_2 / nCells_2)
            pos_x_2 = np.linspace(\
                            0.5 * L_2 / nCells_2, L_2 * (1 - 0.5/nCells_2), nCells_2)
            A_c_2 = m.pi * (0.25 * (D_throat + D_right) + 0.25 * (D_throat - D_right) * np.cos(pos_x_2 * m.pi / L_2)) ** 2
            dV_2 = np.multiply(dx_2, A_c_2)
            A_s_2 = 2 * m.pi * np.multiply(dx_2, 0.25 * (D_throat + D_right) + 0.25 * (D_throat - D_right) * np.cos(pos_x_2 * m.pi / L_2))
            
            self.GEO["nCells"] = nCells_1 + nCells_2
            self.GEO["dx"] = np.hstack((dx_1, dx_2))
            self.GEO["pos_x"] = np.hstack((pos_x_1, L_1 + pos_x_2))
            self.GEO["A_c"] = np.hstack((A_c_1, A_c_2))
            self.GEO["dV"] = np.hstack((dV_1, dV_2))
            self.GEO["A_s"] = np.hstack((A_s_1, A_s_2))

        elif Geometry[0] == "QuadraticDoubleNozzle":
            D_left, D_throat, D_right, L_1, L_2, nCells_1, nCells_2 = Geometry[1]
            
            dx_1 = np.full(nCells_1, L_1 / nCells_1)
            pos_x_1 = np.linspace(\
                            0.5 * L_1 / nCells_1, L_1 * (1 - 0.5/nCells_1), nCells_1)
            A_c_1 = m.pi * (0.5 * D_throat + 0.5 * (D_left - D_throat) * (1 - pos_x_1 / L_1) ** 2) ** 2
            dV_1 = np.multiply(dx_1, A_c_1)
            A_s_1 = 2 * m.pi * np.multiply(dx_1, 0.5 * D_throat + 0.5 * (D_left - D_throat) * (1 - pos_x_1 / L_1) ** 2)

            dx_2 = np.full(nCells_2, L_2 / nCells_2)
            pos_x_2 = np.linspace(\
                            0.5 * L_2 / nCells_2, L_2 * (1 - 0.5/nCells_2), nCells_2)
            A_c_2 = m.pi * (0.5 * D_right + 0.5 * (D_throat - D_right) * (1 - pos_x_2 / L_2) ** 2) ** 2
            dV_2 = np.multiply(dx_2, A_c_2)
            A_s_2 = 2 * m.pi * np.multiply(dx_2, 0.5 * D_right + 0.5 * (D_throat - D_right) * (1 - pos_x_2 / L_2) ** 2)            

            self.GEO["nCells"] = nCells_1 + nCells_2
            self.GEO["dx"] = np.hstack((dx_1, dx_2))
            self.GEO["pos_x"] = np.hstack((pos_x_1, L_1 + pos_x_2))
            self.GEO["A_c"] = np.hstack((A_c_1, A_c_2))
            self.GEO["dV"] = np.hstack((dV_1, dV_2))
            self.GEO["A_s"] = np.hstack((A_s_1, A_s_2))
        

    def fillBasicProps(self, pBounds, vel_xBounds, rhoBounds, interfaceIndex):
        nCells = self.GEO["nCells"]
        self.props["psi"] = np.full(nCells, np.nan)
        self.props["phi"] = np.full(nCells, np.nan)
        self.props["p"] = np.full(nCells, np.nan)
        self.props["vel_x"] = np.full(nCells, np.nan)
        interfacePosition = 0.5 * (self.GEO["pos_x"][interfaceIndex - 1] + self.GEO["pos_x"][interfaceIndex])
        for index, pos_x in enumerate(self.GEO["pos_x"]):
            dx = self.GEO["dx"][index]
            signedDist = pos_x - interfacePosition
            alpha = 2.5 * dx
            if signedDist > alpha:
                psi = 1
            elif signedDist < -alpha:
                psi = 0
            else:
                psi = 0.5 * (1 + (signedDist / alpha + 1 / m.pi * m.sin(m.pi * signedDist / alpha)))
            if abs(psi) <1e-6:
                self.props["psi"][index] = 0.0
            elif abs(psi - 1.0) < 1e-6:
                self.props["psi"][index] = 1.0
            else:
                self.props["psi"][index] = psi
            self.props["phi"][index] = signedDist
            if index <= interfaceIndex:
                p = pBounds[0]
                vel_x = vel_xBounds[0]
            else:
                p = pBounds[1]
                vel_x = vel_xBounds[1]
            self.props["p"][index] = p
            self.props["vel_x"][index] = vel_x
        rho_L = np.full(nCells, rhoBounds[1])
        rho_G = np.full(nCells, rhoBounds[0])
        self.props["rho_G"] = rho_G
        self.props["rho_L"] = rho_L
    
    def fillRestOfPropsFromRHOPVEL_X(self):
        nCells = self.GEO["nCells"]
        T_G = np.full(nCells, np.nan)
        T_L = np.full(nCells, np.nan)
        u_G = np.full(nCells, np.nan)
        u_L = np.full(nCells, np.nan)
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

        for cell in range(nCells):
            ### Fill gas props
            gmodel = GasModel("ideal-air-gas-model.lua")
            cellState_G = GasState(gmodel=gmodel)
            p = self.props["p"][cell]
            rho_G = self.props["rho_G"][cell]
            cellState_G.p = p
            cellState_G.rho = rho_G
            cellState_G.update_thermo_from_rhop()
            cellState_G.update_sound_speed()
            T_G[cell] = cellState_G.T
            u_G[cell] = cellState_G.u
            h_G[cell] = gmodel.enthalpy(cellState_G)
            Cp_G[cell] = cellState_G.Cp
            Cv_G[cell] = cellState_G.Cv
            R_G[cell] = gmodel.R(cellState_G)
            gamma_G[cell] = gmodel.gamma(cellState_G)
            a_G[cell] = cellState_G.a

            ### Fill liquid props
            coeffs = coefficients["water"]
            liquidProps = {}
            liquidProps["rho"] = self.props["rho_L"][cell]
            liquidProps["p"] = self.props["p"][cell]
            liquidProps = LiquidEOS.update_from_rhoP(state = liquidProps, coefficients = coeffs)
            T_L[cell] = liquidProps["T"]
            u_L[cell] = liquidProps["u"]
            h_L[cell] = liquidProps["h"]
            Cp_L[cell] = liquidProps["Cp"]
            Cv_L[cell] = liquidProps["Cv"]
            gamma_L[cell] = liquidProps["gamma"]
            a_L[cell] = liquidProps["a"]
            
        self.props["T_G"] = T_G
        self.props["T_L"] = T_L
        self.props["u_G"] = u_G
        self.props["u_L"] = u_L
        self.props["h_G"] = h_G
        self.props["h_L"] = h_L
        self.props["Cp_G"] = Cp_G
        self.props["Cp_L"] = Cp_L
        self.props["Cv_G"] = Cv_G
        self.props["Cv_L"] = Cv_L
        self.props["R_G"] = R_G
        self.props["gamma_G"] = gamma_G
        self.props["gamma_L"] = gamma_L
        self.props["a_G"] = a_G
        self.props["a_L"] = a_L

