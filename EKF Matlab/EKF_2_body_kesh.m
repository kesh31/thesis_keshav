% Satellite Orbit Propagation and EKF Tracking with Rotating Earth
clear; clc; close all;

% Load input data file
params = load("orbit_model_inputs_radec.mat", 'P0', 'X0_true', ...
              'Rk','GM','Q','Re','dtheta', 'stat_ecef', 'theta0');

meas = load("orbit_model_meas_radec", "tvec", "obs_data");

% Constants
mu = params.GM * 1e9 ; % GM [m^3/s^2]
Re = params.Re * 1e3 ;         % Earth radius [m]
omega_earth = params.dtheta ; % Earth rotation rate [rad/s]
dt = 60;              % Time step [s]
T = 60*60*4;             % Total simulation time [s]
N = T / dt;           % Number of time steps
utc_current = [2000 1 1 0 0 0];

observations = meas.obs_data;


% Initial satellite state [position; velocity] in ECI
% r0 = [7000e3; 0; 0];
% v0 = [0; 7.5e3; 1e3];
x_true = params.X0_true * 1000;

% --- Ground Station in ECEF ---
% gs_lat = 0; gs_lon = 0; gs_alt = 0;    % deg deg m
% gs_lla = [gs_lat, gs_lon, gs_alt];
% gs_ecef = lla2ecef(gs_lla);
gs_ecef = params.stat_ecef * 1000;
[gs_lat, gs_lon, gs_alt] = ecef2lla(gs_ecef'); % For ENU transformation

% --- EKF Initialization ---
x_est = x_true + [100; -100; 50; 0.1; -0.1; 0.05]; % Initial guess
P_est = params.P0;
Q = diag([diag(params.Q);0;0;0]);
R = params.Rk;

% --- Storage ---
x_true_store = zeros(6, N);
x_est_store = zeros(6, N);
ra_est_store = zeros(1, N);
dec_est_store = zeros(1, N);
ra_true_store = zeros(1,N);
dec_true_store = zeros(1,N);

P_est_store = zeros(6,6*N);

% --- Simulation Loop ---
for k = 1:N
    t = (k-1)*dt;

    % Ground station in ECI
    theta = params.theta0 + omega_earth * t;
    gs_eci = [cos(theta),-sin(theta), 0;sin(theta),cos(theta),0;0,0,1] * gs_ecef';

    % True propogation
    x_true = rk4(@(t, x) twobody(t, x, mu), t, x_true, dt);
    [ra_true, dec_true] = measure(x_true, gs_eci);

    % Observation
    % z_obs = [ra_true; dec_true] + sqrt(R) * randn(2, 1); % Add noise
    z_obs = observations(:,k);

    % --- EKF Prediction ---
    x_pred = rk4(@(t, x) twobody(t, x, mu), t, x_est, dt);
    F = stm_jacobian(x_est, mu);
    Phi = eye(6) + F*dt;
    P_pred = Phi * P_est * Phi' + Q;

    % --- EKF Update ---
    [z_pred, H] = measurement_model(x_pred, gs_eci);
    y_innovation = z_obs - z_pred; 
    y_innovation(1) = wrapToPi(y_innovation(1));
    S_cov = H * P_pred * H' + R;
    K_gain = P_pred * H' / S_cov;
    x_est = x_pred + K_gain * y_innovation;
    P_est = (eye(6) - K_gain * H) * P_pred;

    % Estimated measurements
    [ra_est,dec_est] = state_to_radec(x_est,gs_eci);

    % --- Store Results ---
    x_true_store(:,k) = x_true;
    x_est_store(:,k) = x_est;
    ra_est_store(k) = ra_est;
    dec_est_store(k) = dec_est;
    ra_true_store(k) = ra_true;
    dec_true_store(k) = dec_true;
    P_est_store(:,(k-1)*6+1:k*6) = P_est;

    % Update UCT for ground station ECI computation
    if utc_current(5) < 59
        utc_current(5) = utc_current(5) + 1;
        
    else
        utc_current(4) = utc_current(4) + 1;
        utc_current(5) = 0;
    end

end

% wrap azimuth
ra_true_store = wrapTo2Pi(ra_true_store);

% --- Plot Errors ---
figure;
subplot(2,1,1);
plot(0:dt:T-dt, vecnorm(x_true_store(1:3,:) - x_est_store(1:3,:)));
ylabel('Position Error [m]'); title('EKF Position Estimation Error');

subplot(2,1,2);
plot(0:dt:T-dt, vecnorm(x_true_store(4:6,:) - x_est_store(4:6,:)));
ylabel('Velocity Error [m/s]'); xlabel('Time [s]');
title('EKF Velocity Estimation Error');

% Plot measurements
figure;
subplot(2,1,1);
plot(0:dt:T-dt,rad2deg(ra_est_store),'_');
hold on
plot(0:dt:T-dt,rad2deg(ra_true_store),'.');
xlabel('Time [s]'); ylabel('Azimuth (degrees)')
legend('Estimated azimuth','True azimuth')
subplot(2,1,2);
plot(0:dt:T-dt,rad2deg(dec_est_store),'_');
hold on
plot(0:dt:T-dt,rad2deg(dec_true_store),'.');
xlabel('Time [s]'); ylabel('Elevation (degrees)')
legend('Estimated Elevation','True Elevation')
sgtitle('Observations vs True Angles')



% % --- 3D Orbit ---
figure;
plot3(x_true_store(1,:), x_true_store(2,:), x_true_store(3,:), 'b');
hold on;
plot3(x_est_store(1,:), x_est_store(2,:), x_est_store(3,:), 'r--');
legend('True Orbit', 'Estimated Orbit');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Satellite Orbit in ECI Frame'); axis equal;

%% --- FUNCTIONS ---

function dx = twobody(~, x, mu)
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    dx = [v; a];
end

function x_next = rk4(f, t, x, dt)
    k1 = f(t, x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt, x + dt * k3);
    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function [ra, dec] = measure(x_sat, gs_eci)
    r_topo = x_sat(1:3) - gs_eci;
    r = norm(r_topo);
    ra = atan2(r_topo(2), r_topo(1));
    dec = asin(r_topo(3)/r);
end

function [h, H] = measurement_model(x, gs_eci)
    pos_topo = x(1:3) - gs_eci;
    range = norm(pos_topo);
    ra = atan2(pos_topo(2), pos_topo(1));
    dec = asin(pos_topo(3)/range);
    h = [ra;dec];
    
    % Numerical Jacobian
    % eps = 1e-5;
    % H = zeros(2,6);

    % Compute observation partial derivatives
    H11 = pos_topo(1)/range;
    H12 = pos_topo(2)/range;
    H13 = pos_topo(3)/range;
    H21 = -pos_topo(2)/(pos_topo(1)^2 * (1 + (pos_topo(2)/pos_topo(1))^2));
    H22 = 1/(pos_topo(1) * (1 + (pos_topo(2)/pos_topo(1))^2));
    H31 = - pos_topo(3) * pos_topo(1)...
          /(range^3 * sqrt(1 - pos_topo(3)^2/range^2));
    H32 = - pos_topo(3) * pos_topo(2)...
          /(range^3 * sqrt(1 - pos_topo(3)^2/range^2));
    H33 = (1/range - pos_topo(3)^2/range^3)...
          / sqrt(1 - (pos_topo(3)/range)^2);

    H_w_range = [H11, H12, H13, 0, 0, 0;
              H21, H22,   0, 0, 0, 0;
              H31, H32, H33, 0, 0, 0];
    H = H_w_range(2:end,:);

    % for i = 1:6
    %     dx = zeros(6,1); dx(i) = eps;
    %     rho_eps = x(1:3) + dx(1:3) - gs_eci;
    %     r_eps = norm(rho_eps);
    %     ra_eps = atan2(rho_eps(2), rho_eps(1));
    %     dec_eps = asin(rho_eps(3)/r_eps);
    %     h_eps = [ra_eps; dec_eps];
    %     H(:,i) = (h_eps - h) / eps;
    % end
end

function F = stm_jacobian(X, GM)
    % r = x(1:3);
    % rnorm = norm(r);
    % I = eye(3);
    % A = -mu * (I/rnorm^3 - 3 * (r*r') / rnorm^5);
    % F = [zeros(3), eye(3); A, zeros(3)];

    r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);

    A41 = -GM * (1/r^3 - 3*X(1)^2/r^5);
    A42 = 3*GM/r^5 * X(1) * X(2);
    A43 = 3*GM/r^5 * X(1) * X(3);
    A52 = -GM * (1/r^3 - 3*X(2)^2/r^5);
    A53 = 3*GM/r^5 * X(2) * X(3);
    A63 = -GM * (1/r^3 - 3*X(3)^2/r^5);

    F = [0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 1, 0;
         0, 0, 0, 0, 0, 1;
         A41, A42, A43, 0, 0, 0;
         A42, A52, A53, 0, 0, 0;
         A43, A53, A63, 0, 0, 0];

end

function enu = eci2enu(rho_eci, lat_deg, lon_deg)
    lat = deg2rad(lat_deg); lon = deg2rad(lon_deg);
    R = [-sin(lon),             cos(lon),              0;
         -sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
          cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat)];
    enu = R * rho_eci;
end

function [ra, dec] = state_to_radec(x_eci, gs_eci)
    % Computes RA and DEC from satellite ECI state and ground station ECI position

    rho = x_eci(1:3) - gs_eci;
    r = norm(rho);

    ra = atan2(rho(2), rho(1));      % RA: angle in xy-plane from x-axis
    dec = asin(rho(3)/r);            % DEC: angle from equatorial plane
end


% function [az, el] = state_to_azel(x_eci, gs_eci, lat_deg, lon_deg)
% % Transforms satellite state in ECI to azimuth and elevation as seen from a ground station
% %
% % Inputs:
% %   x_eci     - 6x1 satellite state vector in ECI [position; velocity]
% %   gs_eci    - 3x1 ground station position vector in ECI [m]
% %   lat_deg   - ground station latitude [deg]
% %   lon_deg   - ground station longitude [deg]
% %
% % Outputs:
% %   az        - azimuth angle [rad], measured clockwise from North
% %   el        - elevation angle [rad], measured from horizon
% 
%     % Relative position from ground station to satellite in ECI
%     rho_eci = x_eci(1:3) - gs_eci;
% 
%     % Convert relative vector to ENU frame
%     enu = eci2enu(rho_eci, lat_deg, lon_deg);
% 
%     % Components of ENU vector
%     east = enu(1);
%     north = enu(2);
%     up = enu(3);
% 
%     % Azimuth and elevation
%     az = atan2(east, north);               % Azimuth: clockwise from North
%     el = asin(up / norm(enu));             % Elevation: angle above horizon
% 
%     % Ensure azimuth is in [0, 2π]
%     if az < 0
%         az = az + 2*pi;
%     end
% end

function [lat, lon, alt] = ecef2lla(ecef)
    a = 6378137; f = 1/298.257223563;
    e2 = f * (2 - f);
    x = ecef(1); y = ecef(2); z = ecef(3);
    lon = atan2(y, x);
    r = hypot(x, y);
    E2 = a^2 - (a*(1-f))^2;
    F = 54 * (a*(1-f))^2 * z^2;
    G = r^2 + (1 - e2) * z^2 - e2 * E2;
    c = (e2^2 * F * r^2) / (G^3);
    s = (1 + c + sqrt(c^2 + 2*c))^(1/3);
    P = F / (3 * (s + 1/s + 1)^2 * G^2);
    Q = sqrt(1 + 2 * e2^2 * P);
    r0 = -(P * e2 * r) / (1 + Q) + sqrt(0.5*a^2*(1 + 1/Q) - P*(1-e2)*z^2/(Q*(1+Q)) - 0.5*P*r^2);
    U = sqrt((r - e2*r0)^2 + z^2);
    V = sqrt((r - e2*r0)^2 + (1 - e2)*z^2);
    z0 = (a*(1-f))^2 * z / (a * V);
    alt = U * (1 - (a*(1-f))^2 / (a * V));
    lat = atan((z + e2*z0) / r);
    lat = rad2deg(lat); lon = rad2deg(lon);
end