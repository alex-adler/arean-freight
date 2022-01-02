import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

# Globals
g0 = 9.81  # Standard gravity

# Mars properties
rho = 0.02  # Density kg/m3
a = 240  # Speed of sound m/s
g = 3.721   # Gravity m/s2
gamma = 1.31  # Adiabatic gas constant

# Shock Table
M_max = 5
M_step = 0.1
beta_max = 60
beta_step = 1


def calc_CL(L, rho, V, S):
    return (2*L)/(rho*V**2*S)


def calc_V(L, rho, C_L, S):
    return ((2*L)/(rho*C_L*S))**.5


def PG(M):
    return 1/(1-M**2)


def LinearSupersonic(M):
    return 1/(M**2-1)


def Ackeret(theta, M):
    return 2 * theta/np.sqrt(M**2-1)


def SupersonicFlatPlate(alpha, M):
    # Sectional lift coeff, Sectional wave drag coeff
    return (4 * alpha/np.sqrt(M**2-1), 4 * alpha**2/np.sqrt(M**2-1))


def PlotFuelWeightISP(F_N, flight_time_s):
    # Test Values
    Isp = np.linspace(500, 1000, 100)

    w = (F_N/(g0*Isp))*flight_time_s*g

    fig, ax = plt.subplots()
    ax.plot(Isp, w)
    ax.plot([0, Isp[-1]], [1000000, 1000000])
    plt.ylabel("Fuel Weight (N)")
    plt.xlabel("Specific Impulse (s)")
    plt.xlim([0, Isp[-1]])
    plt.ylim([0, w[0]])
    plt.show()


def Oblique(beta, M, gamma):
    beta = np.radians(beta)
    return np.degrees(np.arctan((2*(M**2*np.sin(beta)**2-1))/(np.tan(beta)*(M**2*(gamma+np.cos(2*beta))+2))))


def GetObliqueShockAngle(M, theta, gamma, beta_max=60, beta_step=.1):
    theta_left = 0
    beta_left = 0
    shock_table = []

    # Create the theta-beta plot
    beta_array = np.arange(1, beta_max, beta_step)
    for beta in beta_array:
        shock_table.append((beta, Oblique(beta, M, gamma)))

    # Look for theta and interpolate when passed
    for value in shock_table:
        if value[1] < theta:
            theta_left = value[1]
            beta_left = value[0]
        else:
            return np.interp(theta, [theta_left, value[1]], [beta_left, value[0]])


# Returns (M2, shock angle, p1/p0, p01/p02, T1/T0, rho1/rho0)
def CalcObliqueShock(M1, theta, gamma):
    beta_rad = np.radians(GetObliqueShockAngle(M1, theta, gamma))
    theta_rad = np.radians(theta)
    M1n2 = M1**2*np.sin(beta_rad)**2
    M2 = np.sqrt(((gamma-1)*M1n2+2) /
                 ((2*gamma*M1n2-(gamma-1))*np.sin(beta_rad-theta_rad)**2))
    p2_p1 = (2*gamma*M1n2-(gamma-1))/(gamma+1)
    rho2_rho1 = ((gamma+1)*M1n2) / \
        ((gamma-1)*M1n2+2)
    T2_T1 = ((2*gamma*M1n2-(gamma-1))*((gamma-1)*M1n2+2)) / \
        ((gamma+1)**2*M1n2)
    p02_p01 = rho2_rho1**(gamma/(gamma-1))*(p2_p1)**(-1/(gamma-1))
    return M2, p2_p1, p02_p01, rho2_rho1, T2_T1


print(CalcObliqueShock(2, 5, gamma))

# Design requirements
range_km = 10000
payload_kg = 30000

# Engine properties (nuclear thermal)
Isp_s = 800
F_N = 250000

# Estimated properties
flight_time_s = 5000
We_Wto = 0.5  # Empty weight fraction
cruise_alpha_rad = np.radians(2.5)


# Vehicle properties
# Weight
W_PL_N = payload_kg * g   # Payload weight
W_C_N = 0  # Crew Weight
W_F_N = (F_N/(g0*Isp_s))*flight_time_s*g
MTOW_N = W_PL_N + W_C_N+W_F_N
MTOW_N *= 1+We_Wto

print(f"MTOW = {MTOW_N:.4} N, Payload Weight = {W_PL_N:.4} N, Fuel Weight = {W_F_N:.4} N ({W_F_N/(1000*g):.4} tonnes), Total Mass = {MTOW_N/(1000*g):.4} tonnes")

print(f"M = {calc_V(MTOW_N, rho, 0.1, 1000)/a:.4}")
print(f"C_L = {calc_CL(MTOW_N, rho, 4*a, 1000):.4}")

print(np.degrees(Oblique(np.radians(16.0965164), 4, 1.31)))
