
from tkinter import W
import scipy.optimize as opt
import math
import numpy as np

class Design:
    def __init__(self):
        self.design_parameters = {}
        self.constraints = {}

    def evaluate_design(self,params):
        #Power out(aerodynamic power, after all drivetrain losses) is just denoted by p_{X}, battery/source power is denoted by p_{X}_bat
        #Design parameters
        mac = params[0]
        b = params[1]
        cl_cruise = params[2]

        #Constants
        rho = 1.225
        g=9.8066

        #Design constants
        cd_cruise = 0.014
        cd_max_lift = 0.014
        cd_0 = 0.012
        e=0.95

        cd_fuse = 0.4

        cl_safe = 1.2
        cl_stall = 1.4

        m_fuse=0.9
        r_fuse=0.06

        foil_thickness=0.12

        t_frac = 0.2

        sig_wings = 0.203
        rho_shear = 22.87

        max_load_factor = 3

        cell_dim = 0.125
        max_cell_power = 3.55

        n_drivetrain = 0.95*0.9*0.8
        p_avail = 200*n_drivetrain

        def mass_builtup(b,mac):
            s = mac*b
            m_LE = 0.15*s*0.58*2
            m_spar = ((0.0254*0.25)*(b*mac*foil_thickness)*100)+(2*b*0.0254*0.25*0.002*1750)
            m_ribs = 0.85*s*0.5*foil_thickness*(100*0.025)
            m_TE = b*1e-6*200
            m_skin = s*0.75*2*0.036584
            return sum([m_LE,m_spar,m_ribs,m_TE,m_skin])

        #Derived
        #Helpers
        def velocity(cl):
            return ((2*weight)/(rho*s*cl))**0.5
        def drag(n,velocity):
            q=0.5*rho*(velocity**2)
            d_parasitic_wing = q*(s*cd_cruise)
            d_parasitic_tail = q*(s*t_frac*cd_0)
            d_parasitic_fuse = q*(cd_fuse*a_fuse)
            d_parasitic = d_parasitic_tail+d_parasitic_wing+d_parasitic_fuse
            d_induced = ((weight*n)**2)/(0.5*rho*(velocity**2)*s*math.pi*e*AR)
            d=d_parasitic+d_induced
            return [d,d_parasitic_wing,d_parasitic_tail,d_parasitic_fuse,d_induced]
        #Geometry
        s=b*mac
        AR = b/mac
        a_fuse = math.pi*(r_fuse**2)
        #Mass
        #m_wing = (s*sig_wings)+(s*mac*foil_thickness*0.5*rho_shear)
        m_wing = 2*mass_builtup(b,mac)
        #m_tail = (s*sig_wings*t_frac)+(s*(t_frac**1.5)*mac*foil_thickness*0.5*rho_shear)#Calculated assuming fiberglass and foam core
        m_tail = mass_builtup(b*(t_frac**0.5),mac*(t_frac**0.5))
        cell_rows = math.floor((0.85*mac)/cell_dim) #Fore 10% chord unusable for high curvature, rear 15% reserved for control surfaces
        cell_cols = math.floor((b-0.12)/cell_dim)
        aileron_cell_loss = math.ceil(((b/2)*0.3)/cell_dim)
        num_cells = cell_rows*cell_cols-(2*aileron_cell_loss)
        m_cells = num_cells*0.007
        m=m_wing+m_tail+m_fuse+m_cells
        weight=m*g
            #Cell power
        p_out_cells = 3.55*num_cells

        #Aero
        v_cruise=velocity(cl_cruise)
        turn_load_factor_cruise = min(cl_safe/cl_cruise, max_load_factor)
        d_cruise_list = drag(1,v_cruise)
        d_cruise = d_cruise_list[0]
        p_cruise = d_cruise*v_cruise
        p_cruise_bat = p_cruise/n_drivetrain
        radius_cruise = (v_cruise**2)/(g*((turn_load_factor_cruise**2)-1)**0.5)

        v_takeoff = velocity(cl_safe)
        v_stall = velocity(cl_stall)

        n_takeoff = 0.5
        d_takeoff=(m*(v_takeoff**3))/(3*n_takeoff*p_avail)

        p_excess_cruise = p_avail-p_cruise
        climb_cruise = p_excess_cruise/weight
        p_excess_cells_cruise = p_out_cells*n_drivetrain - p_cruise_bat
        climb_solar_cruise = p_excess_cells_cruise/weight
        p_excess_cells_bat = p_out_cells - p_cruise_bat
        

        self.design_parameters["Span m"] = b
        self.design_parameters["MAC m"] = mac
        self.design_parameters["AR"] = AR
        self.design_parameters["Cl_cruise"] = cl_cruise
        self.design_parameters["Mass kg"] = m
        self.design_parameters["Mass wing builtup construction kg"] = m_wing
        self.design_parameters["Mass wing foamcore construction kg"] = (s*sig_wings)+(s*mac*foil_thickness*0.5*rho_shear)
        self.design_parameters["Wing Area Density kg/m^2"] = m_wing/s
        self.design_parameters["Takeoff Dist"] = d_takeoff
        self.design_parameters["n rows"] = cell_rows
        self.design_parameters["n cols"] = cell_cols
        self.design_parameters["n cells"] = num_cells


        self.design_parameters["Cruise Drag N"] = d_cruise
        self.design_parameters["Cruise Induced Drag N"] = d_cruise_list[4]
        self.design_parameters["Cruise Parasitic Wing Drag N"] = d_cruise_list[1]
        self.design_parameters["Cruise Parasitic Tail Drag N"] = d_cruise_list[2]
        self.design_parameters["Cruise Parasitic Fuse Drag N"] = d_cruise_list[3]
        self.design_parameters["Cruise L/D"] = m*g/d_cruise
        self.design_parameters["Cruise Velocity m/s"] = v_cruise
        self.design_parameters["Stall Velocity m/s"] = v_stall
        self.design_parameters["Cruise Radius m"] = radius_cruise
        self.design_parameters["Cruise Power W"] = p_cruise_bat
        self.design_parameters["Cell Power Margin Bat W"] = p_excess_cells_bat
        self.design_parameters["Cell Power Margin Out W"] = p_excess_cells_cruise
        self.design_parameters["Cruise Solar Climb m/s"] = climb_solar_cruise
        self.design_parameters["Cruise Burst Climb m/s"] = climb_cruise

        maximum_radius = 10
        self.constraints["Con1: Cruise Radius"] = maximum_radius - radius_cruise
        min_stall_margin = 1.3
        self.constraints["Con2: Stall Margin"] = (v_cruise/v_stall) - min_stall_margin
        self.constraints["Con3: Even Cell Col Count"] = -1*(cell_cols%2)
        max_cells = 32
        self.constraints["Con4: Cell Count"] = max_cells - num_cells

    def objective(self, params):
        self.evaluate_design(params)
        return self.design_parameters["Cruise Solar Climb m/s"]*-1

    def con1(self, params):
        self.evaluate_design(params)
        return self.constraints['Con1: Cruise Radius']

    def con2(self,params):
        self.evaluate_design(params)
        return self.constraints["Con2: Stall Margin"]

    def con3(self,params):
        self.evaluate_design(params)
        return self.constraints["Con3: Even Cell Col Count"]

    def con4(self,params):
        self.evaluate_design(params)
        return self.constraints["Con4: Cell Count"]

class landing_gear:
    def __init__(self):
        self.design_parameters = {}
        self.constraints = {}

    def evaluate_design(self,params):
        """
        params = 1D arraylike:
        0 = t m - cap thickness
        1 = w m - width
        2 = L m - length
        3 = theta rad - angle with ground
        4 = c m - core thickness
        5 = lam {0} - taper ratio
        """
        t = params[0]
        w = params[1]
        L = params[2]
        theta = params[3]
        c = params[4]
        lam = params[5]

        E_cap = 100e9 #assume 46.5 GPa modulus for unidirectional bulk E glass composite
        #All properties of 54%SiO2-15%Al2O3-12%CaO E glass fiber (https://www.azom.com/properties.aspx?ArticleID=764)
        E_core = 0.3e9 #Radial compressive modulus of medium density balsa
        v = 7.78*math.sin(10*math.pi/180)+0.3 #0.3 m/s for downdrafts
        m_tot = 6.48

        z_2 = ((2*t)+c)/2
        t_tot = z_2*2
        z_1 = (c/2)
        EI = 2*((E_cap*(2*w/3)*((z_2**3)-(z_1**3)))+(E_core*(w*(c**3)/12))) #I is doubled to account for two landing gear pylons

        a = EI/w
        del_x = w*(lam-1)/(lam+1)
        root = del_x+w
        tip = -del_x+w
        T = 2*del_x/L
        R=w+del_x
        negt = L*T-R
        c2 = (L/R)/2
        c3 = (negt/(R**2))/6
        c4 = (2*T*negt/(R**3))/24
        c5 = (6*(T**2)*negt/(R**4))/120

        #k = 3*(EI)/(L**3) 
        k=a/(c2*(L**2)+c3*(L**3)+c4*(L**4)+c5*(L**5))
        F_max = v*((m_tot*k)**0.5)
        z_accel_max = ((1/math.cos(theta))*F_max)/m_tot
        F_c = F_max*(math.tan(theta))
        D_max = F_max/k
        rho_min = (a*root)/(F_max*L)
        max_z = (c/2)+t
        sigma_max = ((E_cap*max_z)/rho_min)+(F_c/(2*t*w))
        m_caps = 2*t*w*L*1550
        m_core = c*w*L*100
        m_struts = 2*(m_caps+m_core)
        m_fairings = 2*((t_tot*2*L*2*0.58)+((t_tot/0.5)*t_tot*L*30))
        m = m_struts+m_fairings
        F_com_crit = math.pi**2*EI*(1/(4*L**2))

        self.design_parameters["mass"] = m
        self.design_parameters["core mass"] = 2*m_core
        self.design_parameters["cap mass"] = 2*m_caps
        self.design_parameters["fairing mass"] = m_fairings
        self.design_parameters["cap thickness"] = t
        self.design_parameters["width"] = w
        self.design_parameters["root"] = root
        self.design_parameters["tip"] = tip
        self.design_parameters["taper ratio"] = lam
        self.design_parameters["length"] = L
        self.design_parameters["theta"] = theta
        self.design_parameters["core thickness"] = c
        self.design_parameters["total thickness"] = t_tot
        self.design_parameters["aspect ratio"] = w/t_tot

        
        self.design_parameters["EI"] = EI
        self.design_parameters["k"] = k
        self.design_parameters["F_max"] = F_max
        self.design_parameters["F_c"] = F_c
        self.design_parameters["D_max"] = D_max
        self.design_parameters["Z_max"] = D_max * math.cos(theta)
        self.design_parameters["Z_height"] = L * math.sin(theta)
        self.design_parameters["z_accel_max"] = z_accel_max
        self.design_parameters["sigma_max"] = sigma_max

        self.constraints["Con1: Max Stress"] = 1200e6 - sigma_max*3 #E glass compressive strength 3500 MPa
        l_boom = 0.745
        h_min_tipback = math.tan(15*(math.pi/180))*l_boom
        self.constraints["Con2: Max Deflection"] = min((0.25*L)-D_max,
                                                       (L*math.sin(theta))-(0.05+(D_max*math.cos(theta))),
                                                       (L*math.sin(theta)+0.05)-h_min_tipback,#tipback
                                                       tip-0.02#Min tip width constraint
                                                       )
        max_neg_load_factor = 3
        self.constraints["Con3: Max z accel"] = max_neg_load_factor*9.8066 - z_accel_max
        self.constraints["Con4: Buckling"] = F_com_crit-F_c

    def objective(self,params):
        self.evaluate_design(params)
        return self.design_parameters["mass"]

    def con1_max_stress(self,params):
        self.evaluate_design(params)
        return self.constraints["Con1: Max Stress"]

    def con2_max_deflection(self,params):
        self.evaluate_design(params)
        return self.constraints["Con2: Max Deflection"]

    def con3_max_z_accel(self,params):
        self.evaluate_design(params)
        return self.constraints["Con3: Max z accel"]

    def con4_buckling(self,params):
        self.evaluate_design(params)
        return self.constraints["Con4: Buckling"]




