import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve
import math
import pandas as pd
from numba import jit

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

dtmax = 1e3*seconds_in_year
dtmin = seconds_in_year


# def mu(T):
#     mub * math.exp(Q * (Tb - T) / R / Tb / T)


def interpolate_solution(new_grid_r, grid_r, T_last, sigma_r_last, sigma_t_last, er_last, et_last, Tb):

    interp_r = np.copy(grid_r)
    interp_r[0] = new_grid_r[0]

    tmp = T_last
    tmp[0] = Tb

    tlast = np.interp(new_grid_r, interp_r, tmp)
    sigmarlast = np.interp(new_grid_r, interp_r, sigma_r_last)
    sigmatlast = np.interp(new_grid_r, interp_r, sigma_t_last)
    erlast = np.interp(new_grid_r, interp_r, er_last)
    etlast = np.interp(new_grid_r, interp_r, et_last)

    return tlast, sigmarlast, sigmatlast, erlast, etlast


def solve_temperature_shell(grid_r, T_last, Tb, Ts, k, rho_i, Cp, H, dt):
    # Solve the energy equation in spherical coordinates
    # VECTOR inputs
    # grid_r - node locations
    # T_last - last timestep solution
    #
    # SCALAR inputs
    # Tb - basal temperature
    # Ts - surface temperature
    # k - thermal conductivity
    # rho - density
    # Cp - heat capacity
    # H - (constant) internal heating
    # dt - timestep

    nr = len(grid_r)

    L = np.zeros((nr, nr))
    R = np.zeros(nr)

    for i in range(nr):
        r = grid_r[i]

        if i == 0:
            drm = grid_r[i + 1] - grid_r[i]
        else:
            drm = grid_r[i] - grid_r[i - 1]

        if i == nr - 1:
            drp = drm
        else:
            drp = grid_r[i + 1] - grid_r[i]

        rA = r + drp / 2
        rB = r - drm / 2
        kA = k  # thermal conductivities
        kB = k
        dr = rA - rB
        coef_plus = -kA * rA ** 2 / r ** 2 / drp / dr
        coef_center = rho_i * Cp / dt + kA * rA ** 2 / r ** 2 / drp / dr + kB * rB ** 2 / r ** 2 / drm / dr
        coef_minus = -kB * rB ** 2 / r ** 2 / drm / dr

        if i == 0:
            L[i, i] = coef_center
            #             L(i,i+1) = coef_plus-coef_minus
            #             R(i) = R(i) - 2*Tb*coef_minus
            R[i] = coef_center * Tb
        elif i == nr - 1:
            L[i, i] = coef_center
            #             L(i,i-1) = coef_minus-coef_plus
            #             R(i) = R(i) - 2*Ts*coef_plus
            R[i] = coef_center * Ts
        else:
            L[i, i] = coef_center
            L[i, i - 1] = coef_minus
            L[i, i + 1] = coef_plus
            R[i] = rho_i * Cp / dt * T_last[i] + H

    T = np.linalg.solve(L, R)

    return T


def solve_stress_viscoelastic_shell(grid_r, mu, sigma_r_last, alpha_dTdotdr, Pex, E, nu, dt):
    # Solver for viscoelastic stresses
    # Subject to basal pressure boundary condition given by Pex
    # surface boundary condition is sigma_r = 0 (free surface)
    #
    # VECTOR inputs
    # grid_r contains nodal coordinates
    # mu contains viscosities at each node
    # sigma_r_last contains previous stress (sigma_r) values
    # alpha_dTdotdr is the alpha_l * d/dt(dT/dr)
    #
    # SCALAR inputs
    # Pex excess pressure
    # E shear modulus
    # nu poisson ratio
    # dt timestep
    #
    #

    nr = len(grid_r)

    M1 = np.zeros((nr, nr))  # coefficients on (dsigma/dt)
    M2 = np.zeros((nr, nr))  # coefficients on (sigma_r)
    R = np.zeros(nr)

    for i in range(nr):
        if i == 0:
            drm = grid_r[i + 1] - grid_r[i]
        else:
            drm = grid_r[i] - grid_r[i - 1]

        if i == nr - 1:
            drp = grid_r[i] - grid_r[i - 1]
        else:
            drp = grid_r[i + 1] - grid_r[i]

        rA = grid_r[i] + drp / 2  # half a cell +
        rB = grid_r[i] - drm / 2  # half a cell -
        drc = rA - rB
        this_mu = mu[i]  # viscosity

        # M1 - coefficients of dsigma_r/dt
        const1 = (3 - 3 * nu) / (2 * E)
        const2 = (1 - nu) / E  # coefficient on d(r/2 d/dr)/dr term
        coef_a = rA / 2 / drp / drc
        coef_b = rB / 2 / drm / drc
        coef_plus = 1 / dt * (const1 / (drp + drm) + const2 * coef_a)
        coef_center = 1 / dt * (-const2 * coef_a - const2 * coef_b)
        coef_minus = 1 / dt * (-const1 / (drp + drm) + const2 * coef_b)
        if i == 0:
            pass
            #                 M1(i,i)   = coef_center
            #                 M1(i,i+1) = coef_plus-coef_minus
            #                 R(i) = R(i) - 2*coef_minus*P
        elif i == nr - 1:
            M1[i, i - 1] = coef_minus - coef_plus
            M1[i, i] = coef_center
            R[i] = R[i] - 2 * coef_plus * 0
        else:
            M1[i, i - 1] = coef_minus
            M1[i, i] = coef_center
            M1[i, i + 1] = coef_plus

        # M2 - coefficients of sigma_r
        if i == 0:
            mu_B = mu[i]
        else:
            mu_B = np.exp(np.mean(np.log([mu[i - 1], mu[i]])))  # viscosity halfway between nodes i,i+1

        if i == nr - 1:
            mu_A = mu[i]
        else:
            mu_A = np.exp(np.mean(np.log([mu[i], mu[i + 1]])))  # viscosity halfway between nodes i,i+1

        coef_plus = 1 / (4 * this_mu) / (drp + drm) + rA / 12 / mu_A / drp / drc
        coef_center = -rA / 12 / mu_A / drp / drc - rB / 12 / mu_B / drm / drc
        coef_minus = -1 / (4 * this_mu) / (drp + drm) + rB / 12 / mu_B / drm / drc
        if i == 0:
            M2[i, i] = coef_center
            M2[i, i + 1] = coef_plus - coef_minus
            R[i] = R[i] - 2 * coef_minus * Pex
        elif i == nr - 1:
            M2[i, i - 1] = coef_minus - coef_plus
            M2[i, i] = coef_center
            R[i] = R[i] - 2 * coef_plus * 0  # surface sigma_r = 0
        else:
            M2[i, i - 1] = coef_minus
            M2[i, i] = coef_center
            M2[i, i + 1] = coef_plus

        if i == 0:
            R[i] = R[i] - coef_minus * 2 * Pex
            M2[i, i + 1] = M2[i, i + 1] - coef_minus
        elif i == nr - 1:
            M2[i, i - 1] = M2[i, i - 1] - coef_plus
            R[i] = R[i] - coef_plus * 0  # no change because sigma_r = 0 at surface

        R[i] = R[i] - alpha_dTdotdr[i]
        #         R(i) = R(i)+alpha_l*(Tdot(i+1)-Tdot(i))/2/drc # this term
        #         includes the coupling to the energy equation - Tdot needs
        #         to be updated

    LHS = (M1 + M2)
    R1term = M1 @ sigma_r_last  # this represents terms involving dsigma/dr at previous timestep

    RHS = (R + R1term)

    LHS[0, :] = 0
    LHS[0, 0] = abs(LHS[1, 1])
    RHS[0] = Pex * LHS[0, 0]
    LHS[-1, :] = 0
    LHS[-1, -1] = abs(LHS[-2, -2])
    RHS[-1] = LHS[-1, -1] * 0

    LHS = sps.csr_matrix(LHS)
    sigma_r = spsolve(LHS, RHS)

    # 4. calculate the tangential stress sigma_t
    # first, calculate dsr/dr
    dsrdr = np.zeros(sigma_r.shape)
    for i in range(1, nr - 1):
        dr = grid_r[i + 1] - grid_r[i - 1]
        dsrdr[i] = (sigma_r[i + 1] - sigma_r[i - 1]) / dr

    sigma_g = Pex - (sigma_r[1] - Pex)
    dsrdr[0] = (sigma_r[1] - sigma_g) / 2 / (grid_r[1] - grid_r[0])  # special formula using ghost value
    sigma_g = 0 - (sigma_r[nr - 2] - 0)
    dsrdr[-1] = (sigma_g - sigma_r[-2]) / 2 / (grid_r[-1] - grid_r[-2])

    sigma_t = sigma_r + (grid_r.T / 2) * dsrdr

    # deviatoric stresses, from Hillier and Squyres (1991) equations A8-9
    sigma_tD = (grid_r / 6) * dsrdr
    sigma_rD = (-grid_r / 3) * dsrdr

    return sigma_r, sigma_t, sigma_rD, sigma_tD


def calculate_stress_curve_at_time(startTime, endTime, gridr):
    # 1. Calculate the amount of basal freeze-on and advance the mesh
    # 2. Solve the heat equation using an implicit method
    # 3. Solve for sigma_r
    # 4. Calculate sigma_t
    # 5. Calculate the radial displacements u(r)

    # initialize solution vectors (IC)
    grid_r = np.copy(gridr)
    sigma_r_last = np.zeros(nr)  # initial stresses
    sigma_t_last = np.zeros(nr)  # initial stresses
    T_last = np.linspace(Tb, Ts, nr)
    er_last = np.zeros(nr)  # strains
    et_last = np.zeros(nr)
    ur_last = np.zeros(nr)  # displacement
    z_last = 0  # total amount of thickening
    figPex_last = 0  # initial overpressure
    Pex_last = 0
    time = startTime

    results = None

    while time < endTime:
        # 1. Calculate basal freeze-on and interpolate old solution onto new mesh
        # calculate heat flux
        Tg = Tb - (T_last[1] - Tb)
        dTdr_b_last = (T_last[1] - Tg) / 2 / (grid_r[1] - grid_r[0])
        qb = -k * dTdr_b_last

        # determine the timestep
        dt = dtmax
        if np.abs(qb / Lf / rho_i * dt) > (grid_r[1] - grid_r[0]) / 2:
            dt = (grid_r[1] - grid_r[0]) / 2 / (qb / Lf / rho_i)

        if dt < dtmin:
            dt = dtmin

        # thickening would be dx/dt = qb/(L*rho_i)
        delta_rb = dt * qb / Lf / rho_i
        z = z_last + delta_rb

        # calculate new ocean pressure (Manga and Wang 2007, equation 5)
        Pex_pred = Pex_last + 3 * Ri ** 2 / beta_w / (Ri ** 3 - Rc ** 3) * (
                delta_rb * (rho_w - rho_i) / rho_w - ur_last[0])  # ur_last because we don't yet know the uplift

        # re-mesh onto new grid
        new_grid_r = np.linspace(Ri - z, Ro, nr)
        T_last, sigma_r_last, sigma_t_last, er_last, et_last = interpolate_solution(new_grid_r, grid_r, T_last,
                                                                                    sigma_r_last, sigma_t_last, er_last,
                                                                                    et_last, Tb)

        grid_r = new_grid_r  # end interpolation step

        # 2. form discrete operators and solve the heat equation
        T = solve_temperature_shell(grid_r, T_last, Tb, Ts, k, rho_i, Cp, H, dt)

        # Compute d(Tdot)/dr
        Tdot = (T - T_last) / dt
        dTdr_b = (T[1] - T[0]) / (grid_r[1] - grid_r[0])
        Tdot[0] = delta_rb * dTdr_b / dt  # this is an approximation to the Eulerian cooling rate at the ocean-ice interface

        dTdotdr = np.zeros(nr)
        for i in range(1, nr - 1):
            dTdotdr[i] = (Tdot[i + 1] - Tdot[i - 1]) / (grid_r[i + 1] - grid_r[i - 1])

        dTdotdr[0] = (Tdot[1] - Tdot[0]) / (grid_r[1] - grid_r[0])
        dTdotdr[-1] = (Tdot[-1] - Tdot[-2]) / (grid_r[-1] - grid_r[-2])

        # 3. Nonlinear loop over pressure.
        # because the ocean pressure depends on the uplift, we make a guess
        # (above). Using this guess, we calculate stresses, strains, and
        # displacements. Then we re-calculate the pressure using the new value
        # of radial displacement. We continue until the pressure used in the
        # calculations has converged to the pressure consistent with the
        # calculated displacement
        for iterator in range(maxiter):
            if iterator > 0:
                Pex = Pex + relaxation_parameter * (Pex_post - Pex)
            else:
                Pex = Pex_pred

            # calculate viscosity at each node
            mu_node = np.zeros(nr)
            for i in range(nr):
                mu_node[i] = mu(T[i])

            sigma_r, sigma_t, sigma_rD, sigma_tD = solve_stress_viscoelastic_shell(grid_r, mu_node, sigma_r_last,
                                                                                   alpha_l * dTdotdr, -Pex, E, nu, dt)

            # 5. Calculate the strains
            dT = T - T_last
            dT[0] = delta_rb * dTdr_b

            dsigma_t = sigma_t - sigma_t_last
            dsigma_r = sigma_r - sigma_r_last

            de_t = 1 / E * (dsigma_t - nu * (dsigma_t + dsigma_r)) + alpha_l * dT + dt / 2 * (
                    sigma_tD / mu_node)  # change in tangential strain
            de_r = 1 / E * (dsigma_r - 2 * nu * dsigma_t) + alpha_l * dT + dt / 2 * (
                    sigma_rD / mu_node)  # HS91 equations A5-6
            er = er_last + de_r
            et = et_last + de_t
            ur = grid_r * et  # radial displacement

            # re-calculate excess pressure using new uplift
            Pex_post = 3 * Ri ** 2 / beta_w / (Ri ** 3 - Rc ** 3) * (
                    z * (rho_w - rho_i) / rho_w - ur[0])  # ur_last because we don't yet know the uplift

            # check for convergence
            if abs(Pex_post - Pex) / abs(Pex) < 1e-2:
                print(
                    f'dt={dt / seconds_in_year} yr, time={(time + dt) / seconds_in_year / 1e6} Myr, Pex_post {Pex_post} Pex {Pex}, converged in {iterator} iterations')
                break
            elif iterator == maxiter - 1:
                raise Exception('Nonlinear loop failed to converge')

        # 6. advance to next time step and plot (if needed)
        sigma_r_last = sigma_r
        sigma_t_last = sigma_t
        T_last = T
        er_last = er
        et_last = et
        z_last = z
        ur_last = ur
        Pex_last = Pex

        time += dt

        result = np.concatenate([(Ro - grid_r).reshape(grid_r.shape[0], 1),
                                    sigma_r.reshape(sigma_r.shape[0], 1),
                                    sigma_t.reshape(sigma_r.shape[0], 1),
                                    er.reshape(et.shape[0], 1),
                                    et.reshape(et.shape[0], 1),
                                    T.reshape(et.shape[0], 1),
                                    ur.reshape(et.shape[0], 1),
                                    np.ones((grid_r.shape[0], 1)) * time,
                                    np.ones((grid_r.shape[0], 1)) * z
                                 ], axis=1)
        if results is None:
            results = result
        else:
            results = np.concatenate([results, result])

    return pd.DataFrame(results, columns=['r', 'sigma_r', 'stress', 'strain', 'et',
                                          'temp', 'displacement', 'time', 'thickness'])


def find_data_at_peak_stress(data, stressThreshold):
    groups = data.groupby('time')
    stressData = groups.filter(lambda x: x.stress.max() >= stressThreshold)

    if stressData.shape[0] == 0:
        return None

    peakedStress = stressData.loc[stressData['time'] == stressData['time'].min()]
    return peakedStress

