# Load required tudatpy modules
import numpy as np
from matplotlib import pyplot as plt
from tudatpy.interface import spice
from tudatpy.astro import time_conversion, element_conversion
from tudatpy.math import interpolators
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup, environment, propagation_setup, estimation, estimation_setup
from tudatpy.numerical_simulation.estimation_setup import observation
from tudatpy.numerical_simulation.environment import Tle
from datetime import datetime, timedelta
import matplotlib.dates as mdatesx  
from itertools import zip_longest
from tudatpy.util import result2array
import math

spice.load_standard_kernels()


def tle_epoch_to_datetime(epoch_str: str) -> datetime:
    """
    Convert a TLE epoch string (YYDDD.FFFFFFFF) to a Python datetime in UTC.
    
    Parameters:
    -----------
    epoch_str : str
        TLE epoch in the format 'YYDDD.FFFFFFFFF', e.g., '25147.20450743'
    
    Returns:
    --------
    datetime (UTC)
    """
    # Parse year (YY)
    yy = int(epoch_str[:2])
    # TLE years: 00-56 => 2000-2056, 57-99 => 1957-1999
    year = 2000 + yy if yy < 57 else 1900 + yy
    
    # Parse day-of-year (DDD)
    day_of_year = int(epoch_str[2:5])
    
    # Fractional part of the day
    frac_day = float("0." + epoch_str.split('.', 1)[1])
    
    # Compute the base date (January 1 of year) + day offset
    base_date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
    # Add the fractional day as a timedelta
    full_date = base_date + timedelta(days=frac_day)

    # Round to the nearest second
    # Convert to seconds since UNIX epoch for rounding
    epoch_origin = datetime(1970, 1, 1)
    total_seconds = (full_date - epoch_origin).total_seconds()
    rounded_seconds = round(total_seconds)
    rounded_date = epoch_origin + timedelta(seconds=rounded_seconds)

    return rounded_date


def sat_prop(
        simulation_start_epoch,
        simulation_end_epoch,
        target,
        target_TLE = None,
        target_initial_state = None,
        fixed_step_size = 60.0
        ):
    """
    Propogates the target object

    Parameters:
    simulation_start_epoch: float
        Simulation start epoch, in seconds from J2000
    simulation_end_epoch: float
        Simulation end eopch, in sends from J2000
    target: str
        Name of celestial body being observed
    target_TLE: tudatpy.kernel.numerical_simulation.environment.Tle object
        TLE of target object. Run TLE values through environment.Tle
    target_initial_state: numpy array
        Contains cartesian position and velocity of target's initial state

    Output:
    states_array: numpy array
        7 column array with time, position and velocities, in m and m/s
    """
    # Create bodies
    bodies_to_create = ['Sun','Earth','Moon']
    body_settings = environment_setup.get_default_body_settings(bodies_to_create,'Earth','J2000')
    body_settings.add_empty_settings(target)
    bodies_to_propagate = [target]
    central_bodies = ["Earth"]

    # # CREATE ACCELERATION SETTINGS
    # # Aerodynamic drag
    ref_area = np.pi * 0.5**2
    drag_coef = 2.0
    aero_coef_settings = environment_setup.aerodynamic_coefficients.constant(ref_area,[drag_coef,0.0,0.0])

    # Solar radiation pressure
    rad_coef = 1.2
    occulting_bodies_dict = dict()
    occulting_bodies_dict["Sun"] = ["Earth"]
    vehicle_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        ref_area, rad_coef, occulting_bodies_dict)
    
    body_settings.get(target).aerodynamic_coefficient_settings = aero_coef_settings
    body_settings.get(target).radiation_pressure_target_settings = vehicle_target_settings
    
    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.get(target).mass = 10

    accelerations_settings_target = dict(
        Sun=[
            propagation_setup.acceleration.radiation_pressure(),
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Earth=[
            propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
            propagation_setup.acceleration.aerodynamic()
        ],
        Moon=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
    )

    # # Create global accelerations settings dictionary.
    acceleration_settings = {target: accelerations_settings_target}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        acceleration_settings,
        bodies_to_propagate,
        central_bodies)
    
    # Get initial states
    if target_initial_state is None:
        target_ephemeris = environment.TleEphemeris("Earth", "J2000", target_TLE, False)
        initial_state = target_ephemeris.cartesian_state(simulation_start_epoch)
    else:
        initial_state = target_initial_state

    # Create termination settings
    termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)

    # Create numerical integrator settings
    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        fixed_step_size, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
    )

    # Create propagation settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        simulation_start_epoch,
        integrator_settings,
        termination_condition,
    )

    # Create simulation object and propagate the dynamics
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(
        bodies, propagator_settings
    )

    # Extract the resulting state and dependent variable history and convert it to an ndarray
    states = dynamics_simulator.propagation_results.state_history
    states_array = result2array(states)

    return states_array

def compute_tle_checksum(tle_line):
    """
    Calculates and returns the checksum value for a TLE line
    """
    checksum = 0
    for c in tle_line[:68]:  # Only the first 68 characters are used
        if c.isdigit():
            checksum += int(c)
        elif c == '-':
            checksum += 1
    return str(checksum % 10)   


def sat_prop_two_body(
    simulation_start_epoch,
    simulation_end_epoch,
    target,
    target_initial_state,
    fixed_step_size = 60.0
):
    """
    Same as sat prop, but does a simple 2 body propogation with no perturbations involved
    """
    bodies_to_create = ['Earth']
    body_settings = environment_setup.get_default_body_settings(bodies_to_create, 'Earth', 'J2000')
    body_settings.add_empty_settings(target)

    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.get(target).mass = 10  # arbitrary nonzero mass

    # Only point-mass gravity
    acceleration_settings = {
        target: {
            "Earth": [propagation_setup.acceleration.point_mass_gravity()]
        }
    }

    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, [target], ["Earth"]
    )

    # Integrator & propagation
    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        fixed_step_size, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
    )

    termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)

    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"], acceleration_models, [target], target_initial_state,
        simulation_start_epoch, integrator_settings, termination_condition
    )

        # Define parameter settings (just initial state for STM)
    parameter_settings = estimation_setup.parameter.initial_states(propagator_settings, bodies)

    parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, bodies)

    # Solve variational equations
    variational_solver = numerical_simulation.create_variational_equations_solver(
        bodies,
        propagator_settings,
        parameters_to_estimate,
        simulate_dynamics_on_creation=True
    )

    state_history = variational_solver.state_history
    stm_history = variational_solver.state_transition_matrix_history

    return state_history, stm_history



def lla_to_ECEF(lat_rad,lon_rad,alt_m,flattening,equitorial_radius):
    """
    Converts from lat,lon,alt (rad,rad,m) to ecef coordinates (m)
    """
    # eccentricity squared
    e2 = 2 * flattening - flattening**2
    # prime vertical radius of curvature
    N_convert = equitorial_radius / np.sqrt(1 - e2 * np.sin(lat_rad)**2)

    # Compute ECEF (ITRF) position
    x = (N_convert + alt_m) * np.cos(lat_rad) * np.cos(lon_rad)
    y = (N_convert + alt_m) * np.cos(lat_rad) * np.sin(lon_rad)
    z = ((1 - e2) * N_convert + alt_m) * np.sin(lat_rad)
    r_ecef = np.array([x, y, z])  # shape (3,)

    return(r_ecef)


def state_eci_to_radec(x_eci,gs_eci):
    """
    Maps 6 dimensional cartesian coordinates in ECI to right ascension,declination
    """
    pos_topo = x_eci[0:3] - gs_eci
    range = np.sqrt(pos_topo[0]**2 + pos_topo[1]**2 + pos_topo[2]**2)
    ra = math.atan2(pos_topo[1],pos_topo[0])
    dec = math.asin(pos_topo[2]/range)

    return ra,dec


def sat_prop_with_stm(sim_start, sim_end, target, initial_state, fixed_step_size = 60.0):
    """
    Similar to sat_prop, but outputs both the state as well as the STM"""
    
    bodies_to_create = ['Sun','Earth','Moon']
    body_settings = environment_setup.get_default_body_settings(bodies_to_create,'Earth','J2000')
    body_settings.add_empty_settings(target)
    bodies_to_propagate = [target]
    central_bodies = ["Earth"]

    # CREATE ACCELERATION SETTINGS
    # Aerodynamic drag
    ref_area = np.pi * 0.5**2
    drag_coef = 2.0
    aero_coef_settings = environment_setup.aerodynamic_coefficients.constant(ref_area,[drag_coef,0.0,0.0])

    # Solar radiation pressure
    rad_coef = 1.2
    occulting_bodies_dict = dict()
    occulting_bodies_dict["Sun"] = ["Earth"]
    vehicle_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        ref_area, rad_coef, occulting_bodies_dict)
    
    body_settings.get(target).aerodynamic_coefficient_settings = aero_coef_settings
    body_settings.get(target).radiation_pressure_target_settings = vehicle_target_settings
    
    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.get(target).mass = 10

    accelerations_settings_target = dict(
        Sun=[
            propagation_setup.acceleration.radiation_pressure(),
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Earth=[
            propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
            propagation_setup.acceleration.aerodynamic()
        ],
        Moon=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
    )

    # Create global accelerations settings dictionary.
    acceleration_settings = {target: accelerations_settings_target}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        acceleration_settings,
        bodies_to_propagate,
        central_bodies)
    

    # Create termination settings
    termination_condition = propagation_setup.propagator.time_termination(sim_end)

    # Create numerical integrator settings
    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        fixed_step_size, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
    )

    # Create propagation settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        sim_start,
        integrator_settings,
        termination_condition,
    )

    # Define parameter settings (just initial state for STM)
    parameter_settings = estimation_setup.parameter.initial_states(propagator_settings, bodies)

    parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, bodies)

    # Solve variational equations
    variational_solver = numerical_simulation.create_variational_equations_solver(
        bodies,
        propagator_settings,
        parameters_to_estimate,
        simulate_dynamics_on_creation=True
    )

    state_history = variational_solver.state_history
    stm_history = variational_solver.state_transition_matrix_history

    return state_history, stm_history

