### MCTS-Based Sensor Tasking in Earth Orbit (TUDAT, Python)
### Simplified Version: 1 Sensor, 1 Object, Ground-Based Observations

# Load standard modules
import numpy as np
from matplotlib import pyplot as plt

# Load tudatpy modules
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup, propagation_setup
from tudatpy.astro import element_conversion, time_conversion
from tudatpy import constants
from tudatpy.util import result2array
from tudatpy.astro.time_conversion import DateTime

# --- INITIALIZATION ---

# Load spice kernels
spice.load_standard_kernels()

# Define simulation start and end epochs
simulation_start_epoch = time_conversion.calendar_date_to_julian_day_since_epoch(2025, 1, 1)
simulation_end_epoch = simulation_start_epoch + 0.01  # ~15 minutes

# Create body settings for Earth and satellite
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create=['Earth'],
    global_frame_origin='Earth',
    global_frame_orientation='J2000'
)

# Create system of bodies
bodies = environment_setup.create_system_of_bodies(body_settings)

# Add ground station
station_altitude = 0.0
station_latitude = np.deg2rad(19.7)   # Haleakala latitude
station_longitude = np.deg2rad(-156.3)

bodies.get_body('Earth').ground_station_model.add_station(
    'Haleakala',
    environment_setup.ground_station.ground_station_settings(
        'Haleakala', [station_longitude, station_latitude, station_altitude]
    )
)

# Define satellite
satellite_initial_state = [7000.0e3, 0, 0, 0, 7.5e3, 1.0e3]  # (x, y, z, vx, vy, vz)

# Create satellite body
bodies.create_empty_body('Satellite')

# Create acceleration settings
acceleration_settings = {
    'Satellite': {
        'Earth': [environment_setup.acceleration.point_mass_gravity()]
    }
}

# Create propagation settings
propagator_settings = numerical_simulation.propagation_setup.propagator_settings.propagator_settings(
    bodies=bodies,
    central_bodies=['Earth'],
    acceleration_models=acceleration_settings,
    bodies_to_propagate=['Satellite'],
    initial_states=np.array(satellite_initial_state),
    termination_settings=numerical_simulation.propagation_setup.propagator_settings.propagation_termination_settings(
        simulation_end_epoch
    )
)

# Create the dynamics simulator
simulator = numerical_simulation.create_dynamics_simulator(
    bodies=bodies,
    propagator_settings=propagator_settings
)

# Extract ephemeris
state_history = simulator.state_history

# Interpolator for states
state_interpolator = interpolators.create_one_dimensional_vector_interpolator(
    state_history,
    interpolators.lagrange_interpolation(8)
)

# --- MCTS ELEMENTS ---

class Node:
    def __init__(self, belief, depth):
        self.belief = belief  # mean + covariance
        self.children = []
        self.visits = 0
        self.total_reward = 0.0
        self.depth = depth

    def add_child(self, child):
        self.children.append(child)

    def get_average_reward(self):
        if self.visits == 0:
            return 0
        return self.total_reward / self.visits

# Observation model: get RA/Dec from satellite relative to station

def simulate_observation(state, station_coords):
    rel_position = state[:3] - station_coords
    r = np.linalg.norm(rel_position)
    ra = np.arctan2(rel_position[1], rel_position[0])
    dec = np.arcsin(rel_position[2]/r)
    return np.array([ra, dec])

# Unscented Kalman update (simplified)

def unscented_update(belief_mean, belief_cov, measurement, R):
    H = np.zeros((2,6))
    H[0,0] = 1.0
    H[1,1] = 1.0
    S = H @ belief_cov @ H.T + R
    K = belief_cov @ H.T @ np.linalg.inv(S)
    innovation = measurement - (H @ belief_mean)
    new_mean = belief_mean + K @ innovation
    new_cov = (np.eye(6) - K @ H) @ belief_cov
    return new_mean, new_cov

# Reward function: negative trace of position covariance

def reward_function(old_cov, new_cov):
    return - (np.trace(new_cov[:3,:3]) - np.trace(old_cov[:3,:3]))

# --- MCTS CORE FUNCTIONS ---


def mcts_search(root, iterations, station_coords, measurement_noise_cov):
    for _ in range(iterations):
        node = root
        path = [node]

        # Selection
        while node.children:
            node = max(node.children, key=lambda n: n.get_average_reward() + np.sqrt(2*np.log(node.visits + 1)/(n.visits+1)))
            path.append(node)

        # Expansion
        current_belief = node.belief

        # Simulate an observation
        true_state = state_interpolator.evaluate(node.depth * 60)  # 60 seconds per step
        measurement = simulate_observation(true_state, station_coords) + np.random.multivariate_normal([0,0], measurement_noise_cov)

        # Belief Update
        new_mean, new_cov = unscented_update(current_belief[0], current_belief[1], measurement, measurement_noise_cov)
        new_belief = (new_mean, new_cov)

        # Create child node
        child = Node(new_belief, node.depth + 1)
        node.add_child(child)
        path.append(child)

        # Simulation / rollout (could be expanded)
        reward = reward_function(current_belief[1], new_cov)

        # Backpropagation
        for ancestor in path:
            ancestor.visits += 1
            ancestor.total_reward += reward

    # Return best action
    best_child = max(root.children, key=lambda n: n.get_average_reward())
    return best_child

# --- MAIN SIMULATION LOOP ---

# Initial belief: satellite initial state + some uncertainty
initial_covariance = np.diag([1e6, 1e6, 1e6, 1e2, 1e2, 1e2])
initial_belief = (np.array(satellite_initial_state), initial_covariance)

# Measurement noise covariance (arcsec-level errors)
measurement_noise_cov = np.deg2rad(np.array([1.0/3600, 1.0/3600]))**2 * np.eye(2)

# Station coordinates in inertial frame (simplified)
station_coords = np.array([6371.0e3, 0, 0])  # Approx Earth radius in meters

# Create MCTS root
root = Node(initial_belief, depth=0)

# Perform MCTS
best_action = mcts_search(root, iterations=100, station_coords=station_coords, measurement_noise_cov=measurement_noise_cov)

# Output
print("Best action belief mean:", best_action.belief[0])
print("Best action belief covariance:", best_action.belief[1])
