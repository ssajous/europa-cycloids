import math


# Numerical parameters
nr = 200 # number of grid points
relaxation_parameter=.01 # used in nonlinear loop.
maxiter=300

# Define physical constants and parameters
# Physical constants
seconds_in_year = 3.1558e7
R=8.314e-3     # in kJ/mol/K
# Boundary conditions and internal heating
H=0 # internal heating rate.
Tb=270
Ts=100
Ro = 2.52e5             # outer radius of ice shell (m)
Ri = 2.52e5-12e3         # inner radius of ice shell (m)
Rc = 1.60e5             # core radius (m)
# Elastic and Viscous propertiespara
E = 5e9        # shear modulus of ice (Pa)
nu = 0.3      # Poisson ratio of ice (-)
beta_w = 4e-10 # Compressibility of water (1/Pa)
alpha_l = 1e-4 # coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
rho_i=900      # density of ice (kg/m^3)
rho_w=1000     # density of water (kg/m^3)
Q=40           # activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e15       # basal viscosity (Pa-s)
mu = lambda T: mub * math.exp(Q * (Tb-T) / R / Tb / T) # function to evaluate viscosity in Pa-s given T
g = 0.113      # used to plot a failure curve
tau = 3e6 # tensile strength, Pa
# Thermal properties
Cp = 2100 #heat capacity of ice J/kg/K
Lf = 334*1000 # latent heat of fusion (J/kg)
kappa = 1e-6# m/s/s
k=kappa*rho_i*Cp
