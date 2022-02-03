class AUSM():
    def __init__(self) -> None:
        pass
    
    def interfaceVelocity(self, a_L, a_R, vel_x_L, vel_x_R):
        if abs(vel_x_L) <= a_L:
            vel_x_L_P = (vel_x_L + a_L) ** 2 / (4 * a_L)
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2
        if abs(vel_x_R) <= a_R:
            vel_x_R_M = - (vel_x_R - a_R) ** 2 / (4 * a_R)
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2

        return vel_x_L_P + vel_x_R_M
    
    def interfacePressure(self, a_L, a_R, vel_x_L, vel_x_R, p_L, p_R):
        if abs(vel_x_L) <= a_L:
            p_L_P = p_L * (vel_x_L + a_L) ** 2 * (2 - vel_x_L/a_L) / (4 * a_L ** 2)
            
        else:
            p_L_P = p_L * (vel_x_L + abs(vel_x_L)) / (2 * vel_x_L)
            
        if abs(vel_x_R) <= a_R:
            p_R_M = p_R * (vel_x_R - a_R) ** 2 * (2 + vel_x_R/a_R) / (4 * a_R ** 2)
            
        else:
            p_R_M = p_R * (vel_x_R - abs(vel_x_R)) / (2 * vel_x_R)
            
        return p_L_P + p_R_M
        
    def interfaceMassFlux(self, a_L, a_R, vel_x_L, vel_x_R, rho_L, rho_R):
        vel_x_interface = self.interfaceVelocity(a_L = a_L, a_R = a_R, \
                                        vel_x_L = vel_x_L, vel_x_R = vel_x_R)

        return 0.5 * ( vel_x_interface * (rho_L + rho_R) \
                                    - abs(vel_x_interface)  * (rho_R - rho_L) )
    
    def interfaceMomentumFlux(self):
        pass

    def interfaceEnthalpyFlux(self):
        pass

class AUSMD():
    def __init__(self) -> None:
        pass
    
    def interfaceVelocity(self, vel_x_L, vel_x_R, a_L, a_R, rho_L, rho_R, p_L, p_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2
            
        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2
            
        return vel_x_L_P + vel_x_R_M

    def interfacePressure(self, vel_x_L, vel_x_R, p_L, p_R, a_L, a_R):
        a_m = max(a_L, a_R)

        if abs(vel_x_L) <= a_m:
            p_L_P = p_L * (vel_x_L + a_m) ** 2 * (2 - vel_x_L/a_m) / (4 * a_m ** 2)
            
        else:
            p_L_P = p_L * (vel_x_L + abs(vel_x_L)) / (2 * vel_x_L)
            
        if abs(vel_x_R) <= a_m:
            p_R_M = p_R * (vel_x_R - a_m) ** 2 * (2 + vel_x_R/a_m) / (4 * a_m ** 2)
            
        else:
            p_R_M = p_R * (vel_x_R - abs(vel_x_R)) / (2 * vel_x_R)
            
        return p_L_P + p_R_M
        
    def interfaceMassFlux(self, rho_L, rho_R, a_L, a_R, p_L, p_R, vel_x_L, vel_x_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) \
                + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2

        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) \
                + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2

        return vel_x_L_P * rho_L + vel_x_R_M * rho_R

    def interfaceMomentumFlux(self, rho_L, a_L, p_L, vel_x_L, \
                                        rho_R, a_R, p_R, vel_x_R):
        massFlux = self.interfaceMassFlux(rho_L = rho_L, rho_R = rho_R, \
            a_L = a_L, a_R = a_R, \
                p_L = p_L, p_R = p_R, \
                    vel_x_L = vel_x_L, vel_x_R = vel_x_R)
        return 0.5 * ((massFlux)*(vel_x_L + vel_x_R) - abs(massFlux) * (vel_x_R - vel_x_L))

    def interfaceEnthalpyFlux(self):
        pass

class AUSMV():
    def __init__(self) -> None:
        pass
    
    def vel_x_L_plus(self, a_L, a_R, vel_x_L, rho_L, rho_R, p_L, p_R):
        a_m = max(a_L, a_R)
        alpha_L = (2 * (p_L/rho_L)) / ( (p_L / rho_L)+ (p_R / rho_R) )
        if abs(vel_x_L) <= a_m:
            vel_x_L_plus = alpha_L * ((vel_x_L + a_m) ** 2 / (4 * a_m)) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
        else:
            vel_x_L_plus = (vel_x_L + abs(vel_x_L)) / 2
        return vel_x_L_plus

    def vel_x_R_minus(self, a_L, a_R, vel_x_R, rho_L, rho_R, p_L, p_R):
        a_m = max(a_L, a_R)
        alpha_R = (2 * (p_R/rho_R)) / ( (p_L / rho_L)+ (p_R / rho_R) )
        if abs(vel_x_R) <= a_m:
            vel_x_R_minus = -alpha_R * ((vel_x_R - a_m) ** 2 / (4 * a_m)) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
        else:
            vel_x_R_minus = (vel_x_R - abs(vel_x_R)) / 2
        return vel_x_R_minus

    def interfaceVelocity(self, vel_x_L, vel_x_R, a_L, a_R, rho_L, rho_R, p_L, p_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2
            
        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2
            
        return vel_x_L_P + vel_x_R_M

    def interfacePressure(self, vel_x_L, vel_x_R, p_L, p_R, a_L, a_R):
        a_m = max(a_L, a_R)

        if abs(vel_x_L) <= a_m:
            p_L_P = p_L * (vel_x_L + a_m) ** 2 * (2 - vel_x_L/a_m) / (4 * a_m ** 2)
            
        else:
            p_L_P = p_L * (vel_x_L + abs(vel_x_L)) / (2 * vel_x_L)
            
        if abs(vel_x_R) <= a_m:
            p_R_M = p_R * (vel_x_R - a_m) ** 2 * (2 + vel_x_R/a_m) / (4 * a_m ** 2)
            
        else:
            p_R_M = p_R * (vel_x_R - abs(vel_x_R)) / (2 * vel_x_R)
            
        return p_L_P + p_R_M
        
    def interfaceMassFlux(self, rho_L, rho_R, a_L, a_R, p_L, p_R, vel_x_L, vel_x_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2

        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2

        return vel_x_L_P * rho_L + vel_x_R_M * rho_R
        
    def interfaceMomentumFlux(self, rho_L, a_L, p_L, vel_x_L, \
                            rho_R, a_R, p_R, vel_x_R, massFlux_L, massFlux_R):
        vel_x_L_P = self.vel_x_L_plus(a_L, a_R, vel_x_L, rho_L, rho_R, p_L, p_R)
        vel_x_R_M = self.vel_x_R_minus(a_L, a_R, vel_x_R, rho_L, rho_R, p_L, p_R)
        return vel_x_L_P * massFlux_L + vel_x_R_M * massFlux_R

    def interfaceEnthalpyFlux(self):
        pass

class AUSMDV():
    def __init__(self) -> None:
        pass

    def vel_x_L_plus(self, rho_L, a_L, p_L, vel_x_L, rho_R, a_R, p_R):
        a_m = max(a_L, a_R)
        alpha_L = (2 * (p_L/rho_L)) / ( (p_L / rho_L)+ (p_R / rho_R) )
        if abs(vel_x_L) <= a_m:
            vel_x_L_plus = alpha_L * ((vel_x_L + a_m) ** 2 / (4 * a_m)) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
        else:
            vel_x_L_plus = (vel_x_L + abs(vel_x_L)) / 2
        return vel_x_L_plus

    def vel_x_R_minus(self, rho_L, a_L, p_L, rho_R, a_R, p_R, vel_x_R):
        a_m = max(a_L, a_R)
        alpha_R = (2 * (p_R/rho_R)) / ( (p_L / rho_L)+ (p_R / rho_R) )
        if abs(vel_x_R) <= a_m:
            vel_x_R_minus = -alpha_R * ((vel_x_R - a_m) ** 2 / (4 * a_m)) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
        else:
            vel_x_R_minus = (vel_x_R - abs(vel_x_R)) / 2
        return vel_x_R_minus
        
    def interfaceVelocity(self, rho_L, a_L, p_L, vel_x_L, rho_R, a_R, p_R, vel_x_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2
            
        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2
            
        return vel_x_L_P + vel_x_R_M

    def interfacePressure(self, a_L, p_L, vel_x_L, a_R, p_R, vel_x_R):
        a_m = max(a_L, a_R)

        if abs(vel_x_L) <= a_m:
            p_L_P = p_L * (vel_x_L + a_m) ** 2 * (2 - vel_x_L/a_m) / (4 * a_m ** 2)
            
        else:
            p_L_P = p_L * (vel_x_L + abs(vel_x_L)) / (2 * vel_x_L)
            
        if abs(vel_x_R) <= a_m:
            p_R_M = p_R * (vel_x_R - a_m) ** 2 * (2 + vel_x_R/a_m) / (4 * a_m ** 2)
            
        else:
            p_R_M = p_R * (vel_x_R - abs(vel_x_R)) / (2 * vel_x_R)
            
        return p_L_P + p_R_M

    def interfaceMassFlux(self, rho_L, rho_R, a_L, a_R, p_L, p_R, vel_x_L, vel_x_R):
        a_m = max(a_L, a_R)
        alpha_L = 2 * (p_L/rho_L) / ((p_L/rho_L) + (p_R/rho_R))
        alpha_R = 2 * (p_R/rho_R) / ((p_L/rho_L) + (p_R/rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_P = alpha_L * (vel_x_L + a_m) ** 2 / (4 * a_m) + (1-alpha_L) * (vel_x_L + abs(vel_x_L)) / 2
            
        else:
            vel_x_L_P = (vel_x_L + abs(vel_x_L)) / 2

        if abs(vel_x_R) <= a_m:
            vel_x_R_M = - alpha_R * (vel_x_R - a_m) ** 2 / (4 * a_m) + (1-alpha_R) * (vel_x_R - abs(vel_x_R)) / 2
            
        else:
            vel_x_R_M = (vel_x_R - abs(vel_x_R)) / 2

        return vel_x_L_P * rho_L + vel_x_R_M * rho_R

    def interfaceMomentumFlux(self, rho_L, a_L, vel_x_L, p_L, \
                                    rho_R, a_R, vel_x_R, p_R, \
                                    massFlux_L, massFlux_R):
        K = 10
        s = min(0.5, abs(p_R - p_L) / min(p_L, p_R))
        AUSMVInterfaceMomentumFlux = AUSMV().interfaceMomentumFlux(\
                    rho_L = rho_L, a_L = a_L, p_L = p_L, vel_x_L = vel_x_L, \
                    rho_R = rho_R, a_R = a_R, p_R = p_R, vel_x_R = vel_x_R, \
                    massFlux_L = massFlux_L, massFlux_R = massFlux_R)
        AUSMDInterfaceMomentumFlux = AUSMD().interfaceMomentumFlux(\
                    rho_L = rho_L, a_L = a_L, p_L = p_L, vel_x_L = vel_x_L, \
                    rho_R = rho_R, a_R = a_R, p_R = p_R, vel_x_R = vel_x_R)

        return 0.5 * (1 + s) * AUSMVInterfaceMomentumFlux \
                    + 0.5 * (1 - s) * AUSMDInterfaceMomentumFlux

    def interfaceEnthalpyFlux(self, rho_L, a_L, p_L, vel_x_L, H_L, rho_R, a_R, p_R, vel_x_R, H_R):
        massFlux = self.interfaceMassFlux(rho_L = rho_L, a_L = a_L, p_L = p_L, vel_x_L = vel_x_L, \
                        rho_R = rho_R,  a_R = a_R,  p_R = p_R,  vel_x_R = vel_x_R)
        return 0.5 * (massFlux * (H_L + H_R) - abs(massFlux) * (H_R - H_L))

class AUSMPlusUpBASIC():
    def __init__(self) -> None:
        self.beta = 1/8
        self.alpha = 3/16
        self.Ku = 0.75
        self.Kp = 0.25
        self.sigma = 1.0
        pass
    
    def MachPolynomial(self, order, P_M, Ma):
        if order == 1:
            if P_M == "P":
                Ma_m_PM =  0.5 * (Ma + abs(Ma))
            elif P_M == "M":
                Ma_m_PM = 0.5 * (Ma - abs(Ma))
        elif order == 2:
            if P_M == "P":
                Ma_m_PM = 0.25 * (Ma + 1) ** 2
            elif P_M == "M":
                Ma_m_PM = - 0.25 * (Ma - 1) ** 2
        elif order == 4:
            if abs(Ma) >= 1:
                if P_M == "P":
                    Ma_m_PM = self.MachPolynomial(order = 1, P_M = "P", Ma = Ma)
                elif P_M == "M":
                    Ma_m_PM = self.MachPolynomial(order = 1, P_M = "M", Ma = Ma)
            else:
                if P_M == "P":
                    Ma_m_PM = self.MachPolynomial(order = 2, P_M = "P", Ma = Ma) \
                        * (1 - 16 * self.beta * self.MachPolynomial(order = 2, P_M = "M", Ma = Ma))
                elif P_M == "M":
                    Ma_m_PM = self.MachPolynomial(order = 2, P_M = "M", Ma = Ma) \
                        * (1 + 16 * self.beta * self.MachPolynomial(order = 2, P_M = "P", Ma = Ma))
        return Ma_m_PM

    def Ma_0(self, vel_x_L, vel_x_R, a_interface):
        M_inf = 0
        Ma_bar_squared = (vel_x_L ** 2 + vel_x_R ** 2)/(2 * a_interface ** 2)
        return (min(1, max(Ma_bar_squared, M_inf)))

    def interfacePressure(self, rho_L, a_L, p_L, vel_x_L, Ma_L, rho_R, a_R, p_R, vel_x_R, Ma_R, a_crit):
        P_5_plus = self.PressurePolynomialOrder5(Ma = Ma_L, P_M = "P")
        P_5_minus = self.PressurePolynomialOrder5(Ma = Ma_R, P_M = "M")
        #a_L_tilder = self.a_L_tilder(a_crit = a_crit, vel_x_L = vel_x_L)
        #a_R_tilder = self.a_R_tilder(a_crit = a_crit, vel_x_R = vel_x_R)
        a_interface = 0.5 * (a_L + a_R)
        #a_interface = self.interfaceSoundSpeed(a_crit = max(abs(vel_x_L), abs(vel_x_R)), vel_x_L = vel_x_L, vel_x_R = vel_x_R)
        P_interface = P_5_plus * p_L + P_5_minus * p_R \
        - self.Ku * P_5_minus * P_5_plus * (rho_L + rho_R) \
            * (vel_x_R - vel_x_L) * a_interface
        return P_interface

    def PressurePolynomialOrder5(self, Ma, P_M):
        if abs(Ma) >= 1:
            if P_M == "P":
                P_5_PM = self.MachPolynomial(order = 1, P_M = "P", Ma = Ma) / Ma
            elif P_M == "M":
                P_5_PM = self.MachPolynomial(order = 1, P_M = "M", Ma = Ma) / Ma
        else:
            if P_M == "P":
                P_5_PM = self.MachPolynomial(order = 2, P_M = "P", Ma = Ma) \
                    * (2 - Ma - 16 * self.alpha * Ma * self.MachPolynomial(order = 2, P_M = "M", Ma = Ma))
            elif P_M == "M":
                P_5_PM = self.MachPolynomial(order = 2, P_M = "M", Ma = Ma) \
                    * (-2 - Ma + 16 * self.alpha * Ma * self.MachPolynomial(order = 2, P_M = "P", Ma = Ma))
        return P_5_PM
    
    def interfaceSoundSpeed(self, a_crit, vel_x_L, vel_x_R):
        a_L_tilder = self.a_L_tilder(a_crit = a_crit, vel_x_L = vel_x_L)
        a_R_tilder = self.a_R_tilder(a_crit = a_crit, vel_x_R = vel_x_R)
        return min(a_L_tilder, a_R_tilder)
    
    def interfaceMachNumber(self, order, rho_L, a_L, p_L, vel_x_L, Ma_L, rho_R, a_R, p_R, vel_x_R, Ma_R):
        rho_interface = 0.5 * (rho_L + rho_R)
        Ma_m_plus_L = self.MachPolynomial(order = order, P_M = "P", Ma = Ma_L)
        Ma_m_minus_R = self.MachPolynomial(order = order, P_M = "M", Ma = Ma_R)
        a_interface = 0.5 * (a_L + a_R)
        """
        a_interface = self.interfaceSoundSpeed(\
                                a_crit = max(abs(vel_x_L), abs(vel_x_R)), \
                                vel_x_L = vel_x_L, vel_x_R = vel_x_R)
        """
        Ma_bar_squared = (vel_x_L ** 2 + vel_x_R ** 2) / (2 * a_interface ** 2)
        M_p = - self.Kp * max(1 - self.sigma * Ma_bar_squared, 0) * (p_R - p_L) / (rho_interface * a_interface ** 2)
        Ma_interface = Ma_m_plus_L + Ma_m_minus_R + M_p
        return Ma_interface
        
    def a_R_tilder(self, a_crit, vel_x_R):
        return a_crit ** 2 / max(a_crit, - vel_x_R)
    
    def a_L_tilder(self, a_crit, vel_x_L):
        return a_crit ** 2 / max(a_crit, vel_x_L)

    def interfaceMassFlux(self, a_interface, Ma_interface, rho_L, rho_R):
        if Ma_interface >= 0:
            interfaceMassFlux =  a_interface * Ma_interface * rho_L
        else:
            interfaceMassFlux =  a_interface * Ma_interface * rho_R
        return interfaceMassFlux

    def interfaceMomentumFlux(self, a_interface, Ma_interface, rho_L, vel_x_L, rho_R, vel_x_R):
        interfaceMassFlux = self.interfaceMassFlux(a_interface = a_interface, Ma_interface = Ma_interface, rho_L = rho_L, rho_R = rho_R)
        if interfaceMassFlux >= 0:
            interfaceMomentumFlux = interfaceMassFlux * vel_x_L
        else:
            interfaceMomentumFlux = interfaceMassFlux * vel_x_R
        return interfaceMomentumFlux

    def interfaceEnthalpyFlux(self, a_interface, Ma_interface, rho_L, H_L, rho_R, H_R):
        interfaceMassFlux = self.interfaceMassFlux(a_interface = a_interface, Ma_interface = Ma_interface, rho_L = rho_L, rho_R = rho_R)
        if interfaceMassFlux >= 0:
            interfaceEnthalpyFlux = interfaceMassFlux * H_L
        else:
            interfaceEnthalpyFlux = interfaceMassFlux * H_R
        return interfaceEnthalpyFlux

class AUSMPlusUpFULL():
    def __init__(self) -> None:
        self.beta = 1/8
        self.Ku = 0.75
        self.Kp = 0.25
        self.sigma = 1.0
        
    
    def MachPolynomial(self, order, P_M, Ma):
        if order == 1:
            if P_M == "P":
                Ma_m_PM =  0.5 * (Ma + abs(Ma))
            elif P_M == "M":
                Ma_m_PM = 0.5 * (Ma - abs(Ma))
        elif order == 2:
            if P_M == "P":
                Ma_m_PM = 0.25 * (Ma + 1) ** 2
            elif P_M == "M":
                Ma_m_PM = - 0.25 * (Ma - 1) ** 2
        elif order == 4:
            if abs(Ma) >= 1:
                if P_M == "P":
                    Ma_m_PM = self.MachPolynomial(order = 1, P_M = "P", Ma = Ma)
                elif P_M == "M":
                    Ma_m_PM = self.MachPolynomial(order = 1, P_M = "M", Ma = Ma)
            else:
                if P_M == "P":
                    Ma_m_PM = self.MachPolynomial(order = 2, P_M = "P", Ma = Ma) \
                        * (1 - 16 * self.beta * self.MachPolynomial(order = 2, P_M = "M", Ma = Ma))
                elif P_M == "M":
                    Ma_m_PM = self.MachPolynomial(order = 2, P_M = "M", Ma = Ma) \
                        * (1 + 16 * self.beta * self.MachPolynomial(order = 2, P_M = "P", Ma = Ma))
        return Ma_m_PM

    def fa(self, Ma_0):
        return Ma_0 * (2 - Ma_0) 

    def Ma_0(self, vel_x_L, vel_x_R, a_interface):
        M_inf = 0
        Ma_bar_squared = (vel_x_L ** 2 + vel_x_R ** 2)/(2 * a_interface ** 2)
        Ma_0 = min(1, max(Ma_bar_squared, M_inf))
        self.alpha = 3/16 * (-4 + 5 * self.fa(Ma_0 = Ma_0) ** 2)
        return Ma_0

    def interfacePressure(self, rho_L, p_L, vel_x_L, Ma_L, rho_R, p_R, vel_x_R, Ma_R, a_crit):
        P_5_plus = self.PressurePolynomialOrder5(Ma = Ma_L, P_M = "P")
        P_5_minus = self.PressurePolynomialOrder5(Ma = Ma_R, P_M = "M")
        a_L_tilder = self.a_L_tilder(a_crit = a_crit, vel_x_L = vel_x_L)
        a_R_tilder = self.a_R_tilder(a_crit = a_crit, vel_x_R = vel_x_R)
        a_interface = self.interfaceSoundSpeed(a_L_tilder = a_L_tilder, a_R_tilder = a_R_tilder)
        Ma_0 = self.Ma_0()
        f_a = self.fa(Ma_0 = Ma_0)
        P_interface = P_5_plus * p_L + P_5_minus * p_R \
        - self.Ku * P_5_minus * P_5_plus * (rho_L + rho_R) \
            * (vel_x_R - vel_x_L) * a_interface * f_a
        return P_interface

    def PressurePolynomialOrder5(self, Ma, P_M):
        if abs(Ma) >= 1:
            if P_M == "P":
                P_5_PM = self.MachPolynomial(order = 1, P_M = "P", Ma = Ma) / Ma
            elif P_M == "M":
                P_5_PM = self.MachPolynomial(order = 1, P_M = "M", Ma = Ma) / Ma
        else:
            if P_M == "P":
                P_5_PM = self.MachPolynomial(order = 2, P_M = "P", Ma = Ma) \
                    * (2 - Ma - 16 * self.alpha * Ma * self.MachPolynomial(order = 2, P_M = "M", Ma = Ma))
            elif P_M == "M":
                P_5_PM = self.MachPolynomial(order = 2, P_M = "M", Ma = Ma) \
                    * (-2 - Ma + 16 * self.alpha * Ma * self.MachPolynomial(order = 2, P_M = "P", Ma = Ma))
        return P_5_PM
    
    def interfaceSoundSpeed(self, a_L_tilder, a_R_tilder):
        return min(a_L_tilder, a_R_tilder)
    
    def interfaceMachNumber(self, order, rho_L, p_L, vel_x_L, Ma_L, rho_R, p_R, vel_x_R, Ma_R, a_crit):
        rho_interface = 0.5 * (rho_L + rho_R)
        Ma_m_plus_L = self.MachPolynomial(order = order, P_M = "P", Ma = Ma_L)
        Ma_m_minus_R = self.MachPolynomial(order = order, P_M = "M", Ma = Ma_R)
        a_L_tilder = self.a_L_tilder(a_crit = a_crit, vel_x_L = vel_x_L)
        a_R_tilder = self.a_R_tilder(a_crit = a_crit, vel_x_R = vel_x_R)
        a_interface = self.interfaceSoundSpeed(a_L_tilder, a_R_tilder)
        Ma_0 = self.Ma_0(vel_x_L = vel_x_L, vel_x_R = vel_x_R, a_interface = a_interface)
        fa = self.fa(Ma_0 = Ma_0)
        Ma_bar_squared = (vel_x_L ** 2 + vel_x_R ** 2) / (2 * a_interface ** 2)
        M_p = - self.Kp / fa * max(1 - self.sigma * Ma_bar_squared, 0) * (p_R - p_L) / (rho_interface * a_interface ** 2)
        Ma_interface = Ma_m_plus_L + Ma_m_minus_R + M_p
        return Ma_interface
        
    def a_R_tilder(self, a_crit, vel_x_R):
        return a_crit ** 2 / max(a_crit, - vel_x_R)
    
    def a_L_tilder(self, a_crit, vel_x_L):
        return a_crit ** 2 / max(a_crit, vel_x_L)

    def interfaceMassFlux(self, a_interface, Ma_interface, rho_L, rho_R):
        if Ma_interface >= 0:
            interfaceMassFlux =  a_interface * Ma_interface * rho_L
        else:
            interfaceMassFlux =  a_interface * Ma_interface * rho_R
        return interfaceMassFlux

    def interfaceMomentumFlux(self, a_interface, Ma_interface, rho_L, vel_x_L, rho_R, vel_x_R):
        interfaceMassFlux = self.interfaceMassFlux(a_interface = a_interface, Ma_interface = Ma_interface, rho_L = rho_L, vel_x_L = vel_x_L, rho_R = rho_R, vel_x_R = vel_x_R)
        if interfaceMassFlux >= 0:
            interfaceMomentumFlux = interfaceMassFlux * vel_x_L
        else:
            interfaceMomentumFlux = interfaceMassFlux * vel_x_R
        return interfaceMomentumFlux

    def interfaceEnthalpyFlux(self, a_interface, Ma_interface, rho_L, H_L, rho_R, H_R):
        interfaceMassFlux = self.interfaceMassFlux(a_interface = a_interface, Ma_interface = Ma_interface, rho_L = rho_L, rho_R = rho_R, H_L = H_L, H_R = H_R)
        if interfaceMassFlux >= 0:
            interfaceEnthalpyFlux = interfaceMassFlux * H_L
        else:
            interfaceEnthalpyFlux = interfaceMassFlux * H_R
        return interfaceEnthalpyFlux