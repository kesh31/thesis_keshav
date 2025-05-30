# Load required tudatpy modules
import numpy as np
from matplotlib import pyplot as plt
from tudatpy.data.horizons import HorizonsQuery
from tudatpy.interface import spice
from tudatpy.astro import time_conversion, element_conversion
from tudatpy.math import interpolators
from tudatpy.numerical_simulation import environment_setup, environment, propagation_setup
from tudatpy.numerical_simulation import estimation, estimation_setup
from tudatpy.numerical_simulation.estimation_setup import observation
from datetime import datetime
import matplotlib.dates as mdates
from itertools import zip_longest

spice.load_standard_kernels()



def plot_combined_elevation(
        target, 
        station_names, 
        start_epoch, 
        end_epoch, 
        time_step = '10m',
        time_format = '%Y-%m-%d %H:%M:%S',
        global_frame_origin = 'Earth', 
        global_frame_orientation = 'J2000', 
        geodetic_positions = None,
        custom_ephemeris = None):

    """
    Plots the combined elevation of a target celestial body as seen from multiple ground stations.

    Parameters:
    -----------
    target : str
        Name of the celestial body being observed.
    station_names : list of str
        List of ground station names from which the target's elevation is to be computed.
    start_epoch : str
        Start epoch of the observation in the specified time format.
    end_epoch : str
        End epoch of the observation in the specified time format.
    time_step : str, optional (default='10m')
        Time step for ephemeris queries (e.g., '10m' for 10 minutes).
    time_format : str, optional (default='%Y-%m-%d %H:%M:%S')
        Format of the input time strings.
    global_frame_origin : str, optional (default='Earth')
        Origin of the global reference frame.
    global_frame_orientation : str, optional (default='J2000')
        Orientation of the global reference frame.
    geodetic_positions : dict, optional (default=None)
        Dictionary containing geodetic positions (altitude, latitude, longitude) for ground stations not present in Tudat's default list.
    custom_ephemeris : tudatpy.ephemeris.Ephemeris, optional (default=None)
        Custom ephemeris to be used for the target body instead of querying Horizons.

    Returns:
    --------
    None
        This function generates a plot of the elevation angles for the target as seen from the specified ground stations.
    """

    bodies_to_create = ["Earth"]

    # Convert start and end epochs to Julian Day
    jd_start_epoch = time_conversion.calendar_date_to_julian_day(datetime.strptime(start_epoch, time_format))
    jd_end_epoch = time_conversion.calendar_date_to_julian_day(datetime.strptime(end_epoch, time_format))
    n_day_buffer = 1

        # Add a buffer to the user-defined start and end epochs.
    # This ensures that the simulation interpolator operates without errors nor warnings.
    # While the buffered epochs will differ from the original user-defined range, 
    # the code will later filter out any dates outside the requested range. 
    calendar_date_simulation_start_epoch = time_conversion.julian_day_to_calendar_date(jd_start_epoch - n_day_buffer)
    calendar_date_simulation_end_epoch = time_conversion.julian_day_to_calendar_date(jd_end_epoch + n_day_buffer)

        # Convert the start and end epochs to seconds since the epoch for simulation. 
    # This conversion is needed, as the 'get_default_body_settings_time_limited' function
    # (as well as other Tudat functions) - which we will use later on -
    # only accept floats (in seconds from J2000) as arguments. 
    simulation_seconds_start_epoch = time_conversion.julian_day_to_seconds_since_epoch(jd_start_epoch - n_day_buffer)
    simulation_seconds_end_epoch = time_conversion.julian_day_to_seconds_since_epoch(jd_end_epoch + n_day_buffer)

    # Actual (user-queried) start epoch. Later on, we will filter our results based on this epoch.
    actual_seconds_start_epoch  = time_conversion.julian_day_to_seconds_since_epoch(jd_start_epoch)
    actual_seconds_end_epoch  = time_conversion.julian_day_to_seconds_since_epoch(jd_end_epoch)


    # Create default Earth and target settings. Using time limited settings helps code run faster
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create, simulation_seconds_start_epoch, simulation_seconds_end_epoch, global_frame_origin, global_frame_orientation)
    body_settings.add_empty_settings(target)

    # # Add Earth's shape settings. 
    # # We go for an oblate spherical shape, using a radius of 6378 km and the current accepted flattening value
    # equatorial_radius = 6378*1e3 # in meters
    # flattening = 1/298
    # body_settings.get('Earth').shape_settings = environment_setup.shape.oblate_spherical(
    #     equatorial_radius = equatorial_radius,
    #     flattening = flattening,
    # )


    # Add IERS Earth Rotation Model
    body_settings.get('Earth').rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs(
        environment_setup.rotation_model.iau_2006, global_frame_orientation,
        interpolators.interpolator_generation_settings_float(interpolators.cubic_spline_interpolation(),
                                                             simulation_seconds_start_epoch, simulation_seconds_end_epoch, 60),
        interpolators.interpolator_generation_settings_float(interpolators.cubic_spline_interpolation(),
                                                             simulation_seconds_start_epoch, simulation_seconds_end_epoch, 60),
        interpolators.interpolator_generation_settings_float(interpolators.cubic_spline_interpolation(),
                                                             simulation_seconds_start_epoch, simulation_seconds_end_epoch, 60))
    
    # Retrieving ephemiris data for target from JPL Horizons 
    # Horizons query to retrieve ephemeris.
    query_ephemerides = HorizonsQuery(
        query_id=target,
        location = '@399',
        epoch_start=calendar_date_simulation_start_epoch,
        epoch_end=calendar_date_simulation_end_epoch,
        epoch_step= time_step
    )

    # Allows for ephemeris retrieval
    horizons_ephemeris = query_ephemerides.create_ephemeris_tabulated(
        frame_origin=global_frame_origin,
        frame_orientation= global_frame_orientation
    )

    # Set the fetched ephemeris as the target's ephemeris
    body_settings.get(target).ephemeris_settings = horizons_ephemeris


    # Initialize the (sub)plots.
    fig, axes = plt.subplots(len(station_names), 1, figsize=(17, 15))
    fig.tight_layout(pad=5)  # Adjust padding between subplots


    # ground_stations_settings = list()
    geodetic_stations_dict = geodetic_positions
    horizons_coordinates_dict = dict()

    for current_station in station_names:
        if not geodetic_positions:
            print(f'Station: {current_station} not found in Tudat list of stations, and no custom geodetic position found in the geodetic_position dictionary.\n'
                    f'Please provide it.\n Aborting...')
            exit()
        elif geodetic_positions and not geodetic_positions[current_station]:
            print(f'Geodetic Position for {current_station} not found in the geodetic_position dictionary. Please provide it.\n Aborting...')
            exit()
        if current_station not in geodetic_positions.keys():
            print(f'Station Location for station: {current_station} not found in Tudat list of stations, '
                    f'nor in ground_station_positions. Please provide a valid geodetic location.'
                    )
            exit()

        else:
            geodetic_position = geodetic_positions[current_station]
            ground_stations_settings_list= environment_setup.ground_station.basic_station(
                current_station,
                [geodetic_position[0], # in meters
                    np.deg2rad(geodetic_position[1]), # latitude in radians
                    np.deg2rad(geodetic_position[2]) # longitude in radians
                    ],
                element_conversion.geodetic_position_type
            )


            horizons_coordinates_dict[current_station] = {
                'lon': geodetic_position[2], # in degrees
                'lat': geodetic_position[1], # in degrees
                'elevation': geodetic_position[0]/1000, # in km
            }

    # In order to visualize the plot better, we also decide to record the rise and set times (if any) over the user-specified timespan.
    # These are defined as the points in time where the target rises above or sets below the local horizon, on a given topocentric location.
    rise_set_times_dict = dict()
    for idx, station_name in enumerate(station_names):
        horizons_coord = horizons_coordinates_dict[station_name]

        # The following query to Horizons and subsequent application of the Tudatpy 'interpolated_station_angles' method 
        # allow to retrieve inteprolated azimuths and elevations.
        # This call to horizons differs from the ones we have done before, in that it does not provide Ephemeris, 
        # but rather a list of azimuths and elevations (observables). 
        if not custom_ephemeris:
            query_az_el = HorizonsQuery(
                query_id=target,
                location=horizons_coord,
                epoch_start=calendar_date_simulation_start_epoch,
                epoch_end=calendar_date_simulation_end_epoch,
                epoch_step= time_step
            )
            
            horizons_antenna_az_el = query_az_el.interpolated_observations(degrees = True)
            start_end_condition = np.logical_and(
                horizons_antenna_az_el[:, 0] >= actual_seconds_start_epoch,
                horizons_antenna_az_el[:, 0] <= actual_seconds_end_epoch
            )            

            horizons_azimuth = horizons_antenna_az_el[:,1][start_end_condition]
            horizons_elevation = horizons_antenna_az_el[:,2][start_end_condition]

            # The 0th column of horizons_antenna_az_el represents the queried times (in seconds).
            # Earlier, we added a buffer to the user-defined start and end epochs in order for our simulation to run smoothly.
            # Notice that we could have added a buffer to the start_epoch only, but in that case we would have gotten some 
            # extrapolation warnings related to the interpolator. Although these do not compromise the final results, 
            # it is good to be aware of it.
            # This is the moment where we filter out the unwanted data, only keeping the times >= actual_seconds_start_epoch
            horizons_seconds_ephemeris = horizons_antenna_az_el[:,0][start_end_condition]

            # Apply conversions for convenience
            horizons_calendar_ephemeris = [time_conversion.seconds_since_epoch_to_julian_day(horizons_second_ephemeris) for horizons_second_ephemeris in horizons_seconds_ephemeris]
            horizons_datetime_times = [time_conversion.julian_day_to_calendar_date(horizons_calendar_time) for horizons_calendar_time in horizons_calendar_ephemeris]

            # We would like the times at which the elevation and azimuth is computed in Tudat to be the same as the Horizons ones
            # This is done for ease of comparison. 
            tudat_seconds_ephemeris = horizons_seconds_ephemeris

        else:
            # If the custom_ephemeris flag is on, the times have to be retrieved from the provided ephemeris file. 
            # Moreover, in this case, we cannot make a call to Horizons, hence we will only get the Tudat elevation and azimuth plots.
            tudat_seconds_ephemeris = np.linspace(simulation_seconds_start_epoch, simulation_seconds_end_epoch, int(time_step.strip('m'))*10)
            tudat_seconds_ephemeris = tudat_seconds_ephemeris[start_end_condition]

        tudat_calendar_ephemeris = [time_conversion.seconds_since_epoch_to_julian_day(tudat_second_ephemeris) for tudat_second_ephemeris in tudat_seconds_ephemeris]
        tudat_datetime_times = [time_conversion.julian_day_to_calendar_date(tudat_calendar_time) for tudat_calendar_time in tudat_calendar_ephemeris]

        # We struggled quite a bit to set all the body settings, but now we can finally create our bodies object.
        bodies = environment_setup.create_system_of_bodies(body_settings)
        print(bodies.get("Earth"))
        environment_setup.add_ground_station(
            bodies.get("Earth"),
            ground_stations_settings_list
        )

        # Notice how - so far - we have not run any simulation yet. 
        # As we approach doing that, it is time to create the link ends for the simulation. 
        # The given antenna is set as a receiver through the function: 'body_reference_point_link_end_id' of the observation module. 
        # (See: https://py.api.tudat.space/en/latest/observation.html#tudatpy.numerical_simulation.estimation_setup.observation.body_reference_point_link_end_id)
        link_ends = {
            observation.receiver: observation.body_reference_point_link_end_id('Earth', station_name),
        }
        # Create a single link definition from the link ends
        link_definition = observation.LinkDefinition(link_ends)

        # Finally, the Tudat function 'compute_target_angles_and_range' can be used to simulate the observed azimuth and elevation 
        # from the given receiver, at the given observation times (notice how we filter tudat_seconds_ephemeris by start_end_condition)
        # as we only want to cover the user-defined timespan. 
        station_id = ('Earth', station_name)
        angles_range_dictionary = estimation.compute_target_angles_and_range(
            bodies,
            station_id,
            target,
            observation_times=tudat_seconds_ephemeris,
            is_station_transmitting=False
        )

        # Initialize lists for azimuth and elevation
        tudat_azimuth_list = []
        tudat_elevation_list = []

        # The following code snippet allows to compute the rise and set times, as seen from each station, 
        # and to populate the rise_set_dict dictionary. 
        set_times = []
        rise_times = []
        keys_list = sorted(angles_range_dictionary.keys())

        initial_elevation = np.rad2deg(angles_range_dictionary[keys_list[0]][0])
        flag = initial_elevation <= 0  # True if initially below horizon

        for idx_key, key in enumerate(keys_list):  # Ensure iteration follows sorted order
            azimuth_deg = np.rad2deg(angles_range_dictionary[key][1]) % 360
            elevation_deg = (np.rad2deg(angles_range_dictionary[key][0]) + 90) % 180 - 90  # Normalize elevation
            tudat_azimuth_list.append(azimuth_deg)
            tudat_elevation_list.append(elevation_deg)

            if elevation_deg > 0:  # Object above horizon
                if flag:  # Transition from below horizon to above
                    rise_times.append(tudat_datetime_times[idx_key])
                    flag = False  # Update flag
            else:  # Object below horizon
                if not flag:  # Transition from above horizon to below
                    set_times.append(tudat_datetime_times[idx_key])
                    flag = True  # Update flag

        rise_set_times_dict[station_name] = [rise_times, set_times]
        rise_times = rise_set_times_dict[station_name][0]
        set_times = rise_set_times_dict[station_name][1]

        # Convert to numpy arrays
        azimuth_array = np.array(tudat_azimuth_list)
        elevation_array = np.array(tudat_elevation_list)

        # The last bit of this Jupyter Notebook deals with visualization of the obtained data. 
        # Antennas observing at low elevation quickly lose performance. Moreover, if the elvation is negative, the target is not observable.
        # That is why we plot horizontal lines at 15 degrees and 0 degrees.
        if len(station_names) == 1:
            if not custom_ephemeris:
                # axes.scatter(np.array(horizons_datetime_times)[horizons_elevation > 15], horizons_elevation[horizons_elevation > 15], label=f"Horizons Elevation",
                                #   s=10, c='pink', alpha=1)
                # axes.scatter(np.array(horizons_datetime_times)[horizons_elevation > 15], horizons_azimuth[horizons_elevation > 15], label=f"Horizons Azimuth",
                                #   s=10, c='green', alpha=1)
                axes.scatter(np.array(horizons_datetime_times)[horizons_elevation <= 15], horizons_elevation[horizons_elevation <= 15],
                                  s=7, c='pink', alpha=0.3)
                axes.scatter(np.array(horizons_datetime_times)[horizons_elevation <= 15], horizons_azimuth[horizons_elevation <= 15],
                                  s=7, c='green', alpha=0.3)
                axes.scatter(np.array(horizons_datetime_times)[horizons_elevation <= 0], horizons_elevation[horizons_elevation <= 0],
                                  s=5, c='pink', alpha=0.1)
                axes.scatter(np.array(horizons_datetime_times)[horizons_elevation <= 0], horizons_azimuth[horizons_elevation <= 0],
                                  s=5, c='green', alpha=0.1)

            axes.scatter(np.array(tudat_datetime_times)[elevation_array > 15], elevation_array[elevation_array > 15],
                              s=5, label='Tudat Elevation', c='blue')
            axes.scatter(np.array(tudat_datetime_times)[elevation_array > 15], azimuth_array[elevation_array > 15],
                              s=5, label='Tudat Azimuth', c='red')
            axes.scatter(np.array(tudat_datetime_times)[elevation_array <= 15], elevation_array[elevation_array <= 15],
                              s=3, c='blue', alpha=0.5)
            axes.scatter(np.array(tudat_datetime_times)[elevation_array <= 15], azimuth_array[elevation_array <= 15],
                              s=3, c='red', alpha=0.5)
            axes.scatter(np.array(tudat_datetime_times)[elevation_array <= 0], elevation_array[elevation_array <= 0],
                              s=1, c='blue', alpha=0.1)
            axes.scatter(np.array(tudat_datetime_times)[elevation_array <= 0], azimuth_array[elevation_array <= 0],
                              s=1, c='red', alpha=0.1)

            # Add a horizontal line for elevation = 15
            axes.axhline(y=15, color='k', linestyle='--', alpha=0.5, label = '15 deg')
            axes.axhline(y=0, color='green', linestyle='--', alpha=0.5, label = '0 deg')

            for i in range(len(set_times)):
                if i == 0:
                    axes.axvline(x = set_times[i], color='grey', linestyle='--', alpha=0.5, label = 'Set times')
                else:
                    axes.axvline(x = set_times[i], color='grey', linestyle='--', alpha=0.5)

            for i in range(len(rise_times)):
                if i == 0:
                    axes.axvline(x = rise_times[i], color='purple', linestyle='--', alpha=0.5, label = 'Rise times')
                else:
                    axes.axvline(x = rise_times[i], color='purple', linestyle='--', alpha=0.5)

            # Add labels, title, and legend
            axes.set_xlabel('UTC Time')
            axes.set_ylabel('Elevation/Azimuth (Degrees)')
            axes.set_title(f'{target} from {station_name}')
            axes.legend(loc='upper right')

            # Use zip_longest to pair rise_times and set_times, filling missing values with None
            for rise_time, set_time in zip_longest(rise_times, set_times, fillvalue=None):
                
                if not set_time and rise_time:  # If only rise_time exists
                    print(f'Rise from {station_name}:', rise_time)
                elif not rise_time and set_time:  # If only set_time exists
                    print(f'Set from {station_name}:', set_time)
                elif not rise_time and not set_time:  # If neither exists
                    print(f'No rise nor set for {station_name} in the selected timespan.')
                else:  # If both exist
                    print(f'Rise from {station_name}:', rise_time, f'Set from {station_name}:', set_time)

            # Customize x-axis ticks (datetime formatting and positioning)
            axes.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d\n%H:%M'))
            axes.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=30, maxticks=30))
            axes.grid()
            axes.tick_params(axis='x', labelsize=7)  # Adjust 'labelsize' to the desired font size

        else:
            # Horizons will be able to give an Elevation plot only when  custom_ephemeris == False or None
            if not custom_ephemeris: 
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation > 15], horizons_elevation[horizons_elevation > 15], label=f"Horizons Elevation",
                                  s=20, c='pink', alpha=0.5, marker = 's')
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation > 15], horizons_azimuth[horizons_elevation > 15], label=f"Horizons Azimuth",
                                  s=20, c='green', alpha=0.5, marker = 's')
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation <= 15], horizons_elevation[horizons_elevation <= 15],
                                  s=15, c='pink', alpha=0.3, marker = 's')
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation <= 15], horizons_azimuth[horizons_elevation <= 15],
                              s=15, c='green', alpha=0.3, marker = 's')
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation <= 0], horizons_elevation[horizons_elevation <= 0],
                                  s=10, c='pink', alpha=0.1, marker = 's')
                axes[idx].scatter(np.array(horizons_datetime_times)[horizons_elevation <= 0], horizons_azimuth[horizons_elevation <= 0],
                              s=10, c='green', alpha=0.1, marker = 's')
                
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array > 15], elevation_array[elevation_array > 15],
                              s=5, label='Tudat Elevation', c='blue')
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array > 15], azimuth_array[elevation_array > 15],
                              s=5, label='Tudat Azimuth', c='black')
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array <= 15], elevation_array[elevation_array <= 15],
                              s=3, c='blue', alpha=0.5)
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array <= 15], azimuth_array[elevation_array <= 15],
                              s=3, c='black', alpha=0.5)
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array <= 0], elevation_array[elevation_array <= 0],
                              s=1, c='blue', alpha=0.1)
            axes[idx].scatter(np.array(tudat_datetime_times)[elevation_array <= 0], azimuth_array[elevation_array <= 0],
                              s=1, c='black', alpha=0.1)

            # Add a horizontal line for elevation = 15
            axes[idx].axhline(y=15, color='k', linestyle='--', alpha=0.5, label = '15 deg')
            axes[idx].axhline(y=0, color='green', linestyle='--', alpha=0.5, label = '0 deg')

            for i in range(len(set_times)):
                if i == 0:
                    axes[idx].axvline(x = set_times[i], color='grey', linestyle='--', alpha=0.5, label = 'Set times')
                else:
                    axes[idx].axvline(x = set_times[i], color='grey', linestyle='--', alpha=0.5)

            for i in range(len(rise_times)):
                if i == 0:
                    axes[idx].axvline(x = rise_times[i], color='purple', linestyle='--', alpha=0.5, label = 'Rise times')
                else:
                    axes[idx].axvline(x = rise_times[i], color='purple', linestyle='--', alpha=0.5)

            # Add labels, title, and legend
            axes[idx].set_xlabel('UTC Time')
            axes[idx].set_ylabel('Elevation/Azimuth (Degrees)')
            axes[idx].set_title(f'{target} Az/El from {station_name}')
            axes[idx].legend(loc='upper right')

            # Use zip_longest to pair rise_times and set_times, filling missing values with None
            for rise_time, set_time in zip_longest(rise_times, set_times, fillvalue=None):
                
                if not set_time and rise_time:  # If only rise_time exists
                    print(f'Rise from {station_name}:', rise_time)
                elif not rise_time and set_time:  # If only set_time exists
                    print(f'Set from {station_name}:', set_time)
                elif not rise_time and not set_time:  # If neither exists
                    print(f'No rise nor set for {station_name} in the selected timespan.')
                else:  # If both existx1
                    print(f'Rise from {station_name}:', rise_time, f'Set from {station_name}:', set_time)
        
            # Customize x-axis ticks (datetime formatting and positioning)
            axes[idx].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d\n%H:%M'))
            axes[idx].xaxis.set_major_locator(mdates.AutoDateLocator(minticks=10, maxticks=20))
            axes[idx].grid()
            axes[idx].tick_params(axis='x', labelsize=7)  # Adjust 'labelsize' to the desired font size


    # Show the combined plot
    plt.show()

# Select station names
# station_names =['KATH12M', 'YARRA12M']
station_names =['YARRA12M']

# Here, we trick tudat by changing the (already defined in the tudat list of stations) position for KATH12M.
# Hence, we will expect to have the same exact plot for both Yarragadee and Katherine.
# We do this in order to show how an already-existing Tudat station location (such as the one of KATH12M) can be overwritten by the user.  
geodetic_positions = {'YARRA12M': [250,-29.0464, 115.3456]}

# Select the global frame and orientation
global_frame_origin = 'Earth'
global_frame_orientation = 'J2000'

# Retrieving Horizons + Tudat Az/El plot for given start and end epochs
start_epoch = '2022-12-19 00:00:00'
end_epoch = '2022-12-20 00:00:00'

plot_combined_elevation(
    '2014 HK129',
    station_names,
    start_epoch = start_epoch,
    end_epoch = end_epoch,
    time_step= '20m',
    geodetic_positions=geodetic_positions,
    custom_ephemeris = None
)