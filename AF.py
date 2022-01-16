import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from pyatmos import coesa76

# Globals
g0 = 9.81  # Standard gravity

# Mars properties
rho = 0.02  # Density kg/m3
a = 240  # Speed of sound m/s
g = 3.721  # Gravity m/s2
gamma = 1.31  # Adiabatic gas constant
surface_pressure = 610

# Shock Table
beta_max = 60
beta_step = 1


def calc_CL(L, rho, V, S):
    return (2 * L) / (rho * V ** 2 * S)


def calc_V(L, rho, C_L, S):
    return ((2 * L) / (rho * C_L * S)) ** 0.5


def DesignVehicle():
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
    W_PL_N = payload_kg * g  # Payload weight
    W_C_N = 0  # Crew Weight
    W_F_N = (F_N / (g0 * Isp_s)) * flight_time_s * g
    MTOW_N = W_PL_N + W_C_N + W_F_N
    MTOW_N *= 1 + We_Wto

    print(
        f"MTOW = {MTOW_N:.4} N, Payload Weight = {W_PL_N:.4} N, Fuel Weight = {W_F_N:.4} N ({W_F_N/(1000*g):.4} tonnes), Total Mass = {MTOW_N/(1000*g):.4} tonnes"
    )

    print(f"M = {calc_V(MTOW_N, rho, 0.1, 1000)/a:.4}")
    print(f"C_L = {calc_CL(MTOW_N, rho, 4*a, 1000):.4}")


def PG(M):
    return 1 / (1 - M ** 2)


def LinearSupersonic(M):
    return 1 / (M ** 2 - 1)


def Ackeret(theta_deg, M):
    return 2 * np.radians(theta_deg) / np.sqrt(M ** 2 - 1)


def SupersonicFlatPlate(alpha, M):
    # Sectional lift coeff, Sectional wave drag coeff
    return (
        4 * np.radians(alpha) / np.sqrt(M ** 2 - 1),
        4 * np.radians(alpha) ** 2 / np.sqrt(M ** 2 - 1),
    )


def PlotFuelWeightISP(F_N, flight_time_s):
    # Test Values
    Isp = np.linspace(500, 1000, 100)

    w = (F_N / (g0 * Isp)) * flight_time_s * g

    fig, ax = plt.subplots()
    ax.plot(Isp, w)
    ax.plot([0, Isp[-1]], [1000000, 1000000])
    plt.ylabel("Fuel Weight (N)")
    plt.xlabel("Specific Impulse (s)")
    plt.xlim([0, Isp[-1]])
    plt.ylim([0, w[0]])
    plt.show()


# Calculates the deflection required for a given Mach number and shock angle beta
def Oblique(beta, M, gamma):
    beta = np.radians(beta)
    return np.degrees(
        np.arctan(
            (2 * (M ** 2 * np.sin(beta) ** 2 - 1))
            / (np.tan(beta) * (M ** 2 * (gamma + np.cos(2 * beta)) + 2))
        )
    )


def GetObliqueShockAngle(M, theta_deg, gamma, beta_max=90, beta_step=0.1):
    theta_left = 0
    beta_left = 0
    shock_table = []

    # Create the theta-beta plot
    beta_array = np.arange(1, beta_max, beta_step)
    for beta in beta_array:
        shock_table.append((beta, Oblique(beta, M, gamma)))

    # Look for theta and interpolate when passed
    for value in shock_table:
        if value[1] < theta_deg:
            theta_left = value[1]
            beta_left = value[0]
        else:
            return np.interp(theta_deg, [theta_left, value[1]], [beta_left, value[0]])

    print("Normal shock")
    return 90


# Returns (M2, shock angle, p2/p1, p02/p01, rho2/rho1, T2/T1)
def CalcObliqueShock(M1, theta_deg, gamma):
    beta_rad = np.radians(GetObliqueShockAngle(M1, theta_deg, gamma))
    theta_rad = np.radians(theta_deg)
    M1n2 = M1 ** 2 * np.sin(beta_rad) ** 2
    M2 = np.sqrt(
        ((gamma - 1) * M1n2 + 2)
        / ((2 * gamma * M1n2 - (gamma - 1)) * np.sin(beta_rad - theta_rad) ** 2)
    )
    p2_p1 = (2 * gamma * M1n2 - (gamma - 1)) / (gamma + 1)
    rho2_rho1 = ((gamma + 1) * M1n2) / ((gamma - 1) * M1n2 + 2)
    T2_T1 = ((2 * gamma * M1n2 - (gamma - 1)) * ((gamma - 1) * M1n2 + 2)) / (
        (gamma + 1) ** 2 * M1n2
    )
    p02_p01 = rho2_rho1 ** (gamma / (gamma - 1)) * \
        (p2_p1) ** (-1 / (gamma - 1))
    return M2, np.degrees(beta_rad), p2_p1, p02_p01, rho2_rho1, T2_T1


# Get coordinates needed to draw a wing on a plot
def GetWingCoords(h, alpha):
    alpha = np.radians(alpha)
    x = [1 - np.cos(alpha), 1]
    y = [h + np.sin(alpha), h]
    return (x, y)


# Get the coordinates needed to draw a shock wave
def GetShockCoords(x1, y1, M, theta_deg, wing_m, wing_c, gamma):
    shock = CalcObliqueShock(M, theta_deg, gamma)

    if y1 == 0:
        if shock[1] > 89.9:
            # Normal shock
            x2 = x1
            y2 = wing_m * x2 + wing_c
            if x2 > 1:
                y2 = 10
        else:
            # Bounce up
            shock_m = np.tan(np.radians(shock[1]))
            shock_c = -shock_m * x1

            x2 = (shock_c - wing_c) / (wing_m - shock_m)
            # If the shock will miss the wing, draw it off the canvas
            if x2 > 1:
                x2 = 2
            y2 = shock_m * x2 + shock_c
    else:
        # Bounce down
        y2 = 0
        x2 = x1 + (y1 - y2) / (np.tan(np.radians(shock[1])))

    return ([x1, x2], [y1, y2], shock)


def CalcPrandtlMeyer(M, gamma):
    M2 = M**2
    return np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(np.sqrt(((gamma - 1) / (gamma + 1)) * (M2 - 1))) - np.arctan(np.sqrt(M2 - 1))


def InvPrandtlMeyer(nu_rad, gamma):
    nuMax = (np.pi / 2) * (np.sqrt((gamma + 1) / (gamma - 1)) - 1)
    # Deflection greater than max turning angle
    if (nu_rad > nuMax):
        return 0

    mach_max = 18
    mach_step = .01
    mach_left = 0
    nu_left = 0
    PM_table = []

    # Create the mach - Prandtl Meyer Table
    mach_array = np.arange(1, mach_max, mach_step)
    for mach in mach_array:
        PM_table.append((mach, CalcPrandtlMeyer(mach, gamma)))

    for value in PM_table:
        if value[1] < nu_rad:
            mach_left = value[0]
            nu_left = value[1]
        else:
            return np.interp(nu_rad, [nu_left, value[1]], [mach_left, value[0]])

    print("Invalid Prandtl Meyer angle")


def CalcExpansion(M1, theta_rad, gamma):
    nu_M2 = theta_rad + CalcPrandtlMeyer(M1, gamma)
    M2 = InvPrandtlMeyer(nu_M2, gamma)
    T2_T1 = 0
    p2_p1 = 0
    rho2_rho1 = 0
    # If the expansion was successful
    if M2 != 0:
        numerator = 1 + (gamma - 1) / 2 * M1**2
        denominator = 1 + (gamma - 1) / 2 * M2**2
        T2_T1 = numerator / denominator
        p2_p1 = np.power((numerator / denominator), (gamma / (gamma - 1)))
        rho2_rho1 = np.power((numerator / denominator), (1 / (gamma - 1)))

    return (M2, T2_T1, p2_p1, rho2_rho1)


# Draw the plot showing how the wing interacts
def DrawSupersonicWing(M, alpha, h, gamma):

    M_inf = M

    fig, ax = plt.subplots()

    # Draw wing
    wing = GetWingCoords(h, alpha)
    wing_m = (wing[1][0] - wing[1][1]) / (wing[0][0] - wing[0][1])
    wing_c = wing[1][0] - wing[0][0] * wing_m
    ax.plot(wing[0], wing[1])

    x = wing[0][0]
    y = wing[1][0]

    pressure_ratio = [(x, y, 1)]

    max_shocks = 20
    shock_count = 0
    # Loop until flow is no longer supersonic or it has reached the end of the wing
    while M > 1:
        shock = GetShockCoords(x, y, M, alpha, wing_m, wing_c, gamma)
        ax.plot(shock[0], shock[1])
        x = shock[0][1]
        y = shock[1][1]
        M = shock[2][0]
        pressure_ratio.append((x, y, shock[2][2]))
        # print(f"Mach Number is now {M}")
        if x > 1 and y != 0:
            break
        shock_count += 1
        if shock_count > max_shocks:
            print("Too many shocks")
            break

    # Calculate lift and drag
    x_last = pressure_ratio[0][0]
    y_last = pressure_ratio[0][1]
    current_pressure = 1
    lift = 0
    drag = 0
    for p in pressure_ratio:
        if p[1] != 0:
            x = p[0]
            y = p[1]
            if x > 1:
                x = 1
                y = h

            lift += (x - x_last) * current_pressure
            drag += (y_last - y) * current_pressure
            # print(
            #     f"Adding {((x-x_last) * current_pressure)/(0.5*gamma*M_inf**2)} of lift coeff")

            x_last = x
            y_last = y
        current_pressure *= p[2]

    alpha_rad = np.radians(alpha)

    # Upper surface
    expansion = CalcExpansion(M_inf, alpha_rad, gamma)
    lift -= (1-pressure_ratio[0][0])*expansion[2]
    drag -= (pressure_ratio[0][1]-h)*expansion[2]

    # Non-dimensionalise
    c_l = lift / (0.5 * gamma * M_inf**2)
    c_d = drag / (0.5 * gamma * M_inf**2)

    unbounded = SupersonicFlatPlate(alpha, M_inf)

    print(f"Lift Coefficient = {c_l:.4f}, Wave Drag Coefficient = {c_d:.4f}")
    print(f"Unbounded flow c_l= {unbounded[0]:.4f}, c_d = {unbounded[1]:.4f}")
    print(
        f"Gound effect increases c_l by factor of {c_l/unbounded[0]:.2f} and c_d by a factor {c_d/unbounded[1]:.2f}"
    )
    print(f"Final pressure ratio is {current_pressure}")

    # Plot formatting
    plt.title(f"Wing section at alpha = {alpha}, M = {M_inf}, and h = {h} c")
    plt.ylabel(r"$\frac{y}{c}$")
    plt.xlabel(r"$\frac{x}{c}$")
    # plt.xlim([-0.5, 1.5])
    # plt.ylim([0, 1.5 * h])
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    # plt.show()


# Plot altitude vs pressure for Earth and Mars atmospheres
def Atmosphere():
    # plt.figure(figsize=(4.5, 8))
    fig, ax = plt.subplots(figsize=(4.5, 8))

    # Earth atmosphere
    h = np.linspace(0, 100, 100)
    atmo = coesa76(h)
    ax.plot(atmo.P / 1000, h, c="royalblue")
    plt.xlim([0, atmo.P[0] / 1000])

    # Mars Atmosphere
    h_mars = np.linspace(0, 100, 100)
    p_mars = 0.699 * np.exp(-0.00009 * h_mars * 1000)
    ax.plot(p_mars, h_mars, c="coral")
    plt.ylim([0, 100])
    # plt.ylabel("Altitude (km)")
    # plt.xlabel("Pressure (kPa)")
    ax.set_xlabel("Pressure (kPa)")
    ax.set_ylabel("Altitude (km)")

    fig.tight_layout()
    plt.savefig("./big_atmo.svg", transparent=True, bbox_inches="tight")


# Plot Mars pressure at Earth equivalent altitude
def AtmosphereZoomedIn():
    # Earth atmosphere
    h = np.linspace(0, 100, 100)
    atmo = coesa76(h)

    # Mars Atmosphere
    h_mars = np.linspace(0, 100, 100)
    p_mars = 0.699 * np.exp(-0.00009 * h_mars * 1000)

    fig, ax1 = plt.subplots(figsize=(4.5, 8))
    # Earth
    color = "royalblue"
    ax1.set_xlabel("Pressure (kPa)")
    ax1.set_ylabel("Earth Altitude (km)", color=color)
    ax1.plot(atmo.P / 1000, h, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.set_ylim([34, 100])

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # Mars
    color = "coral"
    ax2.set_ylabel("Mars Altitude (km)", color=color)
    ax2.plot(p_mars, h_mars, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    plt.xlim([0, 0.7])
    plt.ylim([0, 100])
    fig.tight_layout()
    plt.savefig("./small_atmo.svg", transparent=True, bbox_inches="tight")


# DesignVehicle()
DrawSupersonicWing(3, 10, 0.02, gamma)
# Atmosphere()
# AtmosphereZoomedIn()
