from collections import namedtuple
import math
AirProperties = namedtuple("air_p", ["Density", "Dynamic_Viscosity", "Thermal_Conductivity", "Specific_Heat"])
import pandas as pd

def calc_re_tube(density, velocity, diameter, viscosity):# calcultes reynolds number for flow inside a tube
    density, velocity, diameter, viscosity = float(density), float(velocity), float(diameter), float(viscosity)
    if viscosity == 0:
        raise ValueError("Viscosity cannot be zero.")
    reynolds_number = (density * velocity * diameter) / viscosity
    return reynolds_number

def calc_re_plane(density, velocity, glassDimeter,od, viscosity): # calculates reynolds number for flow over the tube
    density, velocity, viscosity = float(density), float(velocity), float(viscosity)
    if viscosity == 0:
        raise ValueError("Viscosity cannot be zero.")
    reynolds_number = (density * velocity * 2*(glassDimeter-od)) / viscosity
    return reynolds_number

def calc_pr(viscosity, specific_heat, thermal_conductivity):# calculates prandtl number
    viscosity, specific_heat, thermal_conductivity = float(viscosity), float(specific_heat), float(thermal_conductivity)
    if thermal_conductivity == 0:
        raise ValueError("Thermal conductivity cannot be zero.")
    prandtl_number = (viscosity * specific_heat) / thermal_conductivity
    return prandtl_number

def calc_nu_tube(re, pr):# calculates nusselt number for flow inside a tube using the correlation for turbulent flow 0.23*Re^0.8*Pr^0.4 or n=4.36 for laminar flow
    re, pr = float(re), float(pr)
    if re > 5000:
        nu = 0.023 * (re ** 0.8) * (pr ** 0.4)
    elif 0 < re < 5000:
        nu = 4.36
    else:
        raise ValueError("Re must be greater than zero")
    return nu

def calc_nu_plane(re, pr):# calculates nusselt number for flow over the tube using flow over plate correlation of 0.037*Re^0.8*Pr^0.43 or 0.66*Re^0.5*Pr^0.33 for laminar flow
    re, pr = float(re), float(pr)
    if re > 5000:
        nu = 0.037 * (re ** 0.8) * (pr ** 0.43)
    elif 0 < re < 5000:
        nu = 0.66 * (re ** 0.5) * (pr ** 0.33)
    else:
        raise ValueError("Re must be greater than zero")
    return nu

def calc_q_overall(id, od, nu_in, k, nu_out, t_led, t_air,D_glass):# calculates heat transfer rate
    id, od, nu_in, k, nu_out, t_led, t_air = map(float, (id, od, nu_in, k, nu_out, t_led, t_air))
    h_in = nu_in * k/id
    q_in = h_in * id * 3.14 * (t_led - t_air)
    h_out = nu_out * k/(2*(D_glass-od))
    q_out = h_out * od * 3.14 * (t_led - t_air)
    q_total = q_in + q_out
    return q_total

def air_properties_calculator(temperature):# calculates air properties

    temperature = float(temperature)
    temperature_k = temperature + 273.15
    R = 287.05  # Specific gas constant for air, J/(kg·K)
    P = 101325  # Atmospheric pressure, Pa
    density = P / (R * temperature_k)
    C1 = 1.458e-6  # kg/(m·s·K^0.5)
    S = 110.4  # Sutherland's constant, K
    dynamic_viscosity = C1 * temperature_k ** 1.5 / (temperature_k + S)
    thermal_conductivity = 0.0241 + 7.38e-5 * temperature  # W/(m·K)
    specific_heat = 1005 + 0.1 * temperature  # J/(kg·K)
    return AirProperties(
        Density=density,
        Dynamic_Viscosity=dynamic_viscosity,
        Thermal_Conductivity=thermal_conductivity,
        Specific_Heat=specific_heat
    )


def colebrook(epsilon, diameter, reynolds_number, tolerance=1e-6): # calculates friction factor using colebrook equation or Darcy–Weisbach if laminar
    
    if reynolds_number < 2300:
       friction_factor=float(64/reynolds_number)
       return friction_factor

    # Initial guess for friction factor
    friction_factor = 0.02

    for _ in range(100):  # Limit iterations to ensure no infinite loop
        # Colebrook equation rearranged for iterative solution
        colebrook_lhs = -2 * math.log10((epsilon / (3.7 * diameter)) + (2.51 / (reynolds_number * math.sqrt(friction_factor))))
        new_friction_factor = 1 / (colebrook_lhs ** 2)

        # Check for convergence
        if abs(new_friction_factor - friction_factor) < tolerance:
            return new_friction_factor

        friction_factor = new_friction_factor

    raise RuntimeError("Failed to converge to a solution for the friction factor.")

def pressure_drop(friction_tube,friction_shell,id,od,d_glass,density,velocity):# calculates pressure drop in the tube and shell using darcy weisbach equation
    
    pressure_drop_tube = (friction_factor_tube * 3 * density * velocity ** 2) / (id)
    pressure_drop_shell = (friction_factor_shell * 3 * density * velocity ** 2) / (2*(d_glass-od))
    pressure_drop = pressure_drop_tube + pressure_drop_shell
    return pressure_drop



    
ID = float(0.012)  # inner diameter of LED tube
OD = float(0.015)  # outside diameter of LED tube
glassDiameter = float(0.02)  # diameter of glass tube
V = float(input("Enter the velocity of air inside LED tube (m/s): "))
T_led = float(input("Enter the temperature of LED tube (°C): "))
T_air = float(25 + (T_led - 25) / 2)  # temperature of air


air = air_properties_calculator(T_air)
re_tube = calc_re_tube(air.Density, V, ID, air.Dynamic_Viscosity)
pr = calc_pr(air.Dynamic_Viscosity, air.Specific_Heat, air.Thermal_Conductivity)
re_plane = calc_re_plane(air.Density, V,glassDiameter,OD, air.Dynamic_Viscosity)
nu_in = calc_nu_tube(re_tube, pr)
nu_out = calc_nu_plane(re_plane, pr)
q_total = calc_q_overall(ID, OD, nu_in, air.Thermal_Conductivity, nu_out, T_led, T_air,glassDiameter)
friction_factor_tube = colebrook(0.001, ID, re_tube)
friction_factor_shell = colebrook(0.001, 2*(glassDiameter-OD), re_plane)

pressure_drop = pressure_drop(friction_factor_tube,friction_factor_shell,ID,OD,glassDiameter,air.Density,V)
print("the pressure drop along the led tube and shell is:   " f"{pressure_drop:.2f} Pa")
pressure_drop_mmh2o = pressure_drop * 0.10197
print("or:  "f"{pressure_drop_mmh2o:.2f} mmH2O")
print("the heat transfer due to convection in the tube and shell is "f"{q_total:.2f} w/m")
cfm=(V*ID**2*3.14/4)*3600*1.7
print("the volumetric flow is:  "f"{cfm:.2f} cfm")
just_input=input("press enter to exit")

