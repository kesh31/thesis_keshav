% Satellite Orbit Propagation and EKF Tracking with Rotating Earth
clear; clc; close all;
% Load stuff
params = load("orbit_model_inputs_radec.mat", 'P0', 'X0_true', ...
              'Rk','GM','Q','Re','dtheta', 'stat_ecef', 'theta0');
meas = load("orbit_model_meas_radec", "tvec", "obs_data");
true_data = load("orbit_model_truth.mat",'Xt_mat','time');
state_data_true = true_data.Xt_mat;

% Constants
mu = params.GM ; % earth grav constant [m^3/s^2]
omega_earth = params.dtheta ; % earth rotation rate [rad/s]
dt = 60;              % time step [s]
T = 60*60*4 + 60;             % sim time [s]
N = T / dt;           % number of time steps

observations = meas.obs_data;

% % Initial satellite state [position; velocity] in ECI
x_true = params.X0_true;   % 1.0e3 * [0.7577,5.222607,4.8515,0.00221321,0.00467834,-0.0053713]
x_true = 1.0e+03 *[5.222607,0.7577,4.8515,0.00467834,0.00221321,-0.0053713]';
% x_true = [6.06900255e+06,6.09157542e+06,1.05883186e+06,-1.25648342e+03,3.83516911e+03,-5.98898761e+03]'/1000;
n = size(x_true, 1); 

% Ground Station in ECEF
gs_ecef = params.stat_ecef;

% EKF Initialization 
P_est = params.P0;
R = params.Rk;

% Storing values
x_true_store = zeros(6, N);
x_est_store = zeros(6, N);
ra_est_store = zeros(1, N);
dec_est_store = zeros(1, N);
ra_true_store = zeros(1,N);
dec_true_store = zeros(1,N);
z_obs_store = zeros(2,N);
P_est_store = zeros(6,6*N);
Kk_norm = zeros(1,N);

% Set up settings for ODE45 solver 
RelTol = 1e-12;                                 % relative tolerance
AbsTol = 1e-12;                                 % absolute tolerance
options = odeset('RelTol', RelTol,'AbsTol', AbsTol); 
Phi0 = reshape((eye(n)), 1, [])'; 
Xref_Stm0 = [x_true; Phi0];
xhat_k_prev = zeros(6,1);
t_pre = true_data.time(1);

% sim + EKF loop
for k = 1:N
    tk = (k-1)*dt;

    % Ground station in ECI
    theta = params.theta0 + omega_earth * tk;
    rot_matrix = [cos(theta),-sin(theta), 0;sin(theta),cos(theta),0;0,0,1];
    gs_eci = rot_matrix * gs_ecef';

    % True propogation
    % x_true = rk4(@(t, x) twobody(t, x, mu), t, x_true, dt);
    if k == 1
        x_true;
    else
        [~,x_out] = ode45(@(t,x) twobody(t,x,mu),[t_pre tk],x_true,options);
        x_true = x_out(end,:)';
    end
    % x_true = state_data_true(:,k);
    [ra_true, dec_true] = measure(x_true, gs_eci);

    % Observation
    zk = [ra_true; dec_true] + sqrt(R) * randn(2, 1); % Add noise
    % zk = observations(:,k);

    % Prediction
    if k == 1
        Xref_stm_out = Xref_Stm0';
    else
        Xref_stm_prev(1:3) * 1000;
        [~, Xref_stm_out] = ode45(@(t, X) int_twobody_stm(t, X, mu), [t_pre tk], Xref_stm_prev,options);
    end
    Xref_k = Xref_stm_out(end, 1:6)';
    Phi_k = reshape(Xref_stm_out(end, 7:end), 6, 6)';  % STM
    % norm(Xref_k)
    % norm(Phi_k,"fro")

    Gammak = getGammaMatrix(t_pre, tk);

    xhat_k_kminus1 = Phi_k * xhat_k_prev;  % Nominal state
    P_k_kminus1 = Phi_k * P_est * Phi_k' + Gammak * params.Q * Gammak';

    % EKF Update
    [Gk, H] = measurement_model(Xref_k, gs_eci);
    yk = zk - Gk; 
    % yk(1) = wrapToPi(yk(1));
    S_cov = H * P_k_kminus1 * H' + R;
    Kk = P_k_kminus1 * H' / S_cov;
    Kk_norm(k) = norm(Kk,"fro");

    xhat_k = xhat_k_kminus1 + Kk * (yk - H*xhat_k_kminus1);
    P_est = (eye(6) - Kk*H) * P_k_kminus1 * (eye(6) - Kk*H)' + Kk*R*Kk';
    x_est = Xref_k + xhat_k;

    if k~=1 && vecnorm(xhat_k) > 1
        Xref_k = x_est;
        xhat_k = zeros(n, 1);
    end

    % Xref_stm_prev = [Xref_k;reshape(Phi_k,[],1)];
    Xref_stm_prev = [Xref_k;reshape(Phi0,[],1)]';
    xhat_k_prev = xhat_k;
    t_pre = tk;

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
    z_obs_store(:,k) = zk;
end

% wrap RA
ra_true_store = wrapTo2Pi(ra_true_store);
ra_est_store = wrapTo2Pi(ra_est_store);

%%

% Plot true position
figure;
subplot(3,1,1)
plot(0:dt:T-dt,x_true_store(1,:))
xlabel('X position (m)'), ylabel('Time')

subplot(3,1,2)
plot(0:dt:T-dt,x_true_store(2,:))
xlabel('Y position (m)'), ylabel('Time')

subplot(3,1,3)
plot(0:dt:T-dt,x_true_store(3,:))
xlabel('Z position (m)'), ylabel('Time')
sgtitle('Evolution of true position')


% Plot Errors
figure;
subplot(2,1,1);
plot(12000:dt:T-dt, vecnorm(x_true_store(1:3,201:end) - x_est_store(1:3,201:end)));
ylabel('Position Error [m]'); title('EKF Position Estimation Error');

subplot(2,1,2);
plot(12000:dt:T-dt, vecnorm(x_true_store(4:6,201:end) - x_est_store(4:6,201:end)));
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

% Plot true angles
figure;
plot(0:dt:T-dt,rad2deg(ra_true_store),'.')
hold on
plot(0:dt:T-dt,rad2deg(dec_true_store),'.')
xlabel('Degrees'),ylabel('Time')
title('True Angles')
legend("Right Ascension","Declination")


figure;
plot3(x_true_store(1,:), x_true_store(2,:), x_true_store(3,:), 'b');
hold on;
plot3(x_est_store(1,:), x_est_store(2,:), x_est_store(3,:), 'r--');
legend('True Orbit', 'Estimated Orbit');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Satellite Orbit in ECI Frame'); axis equal;

figure;
subplot(3,1,1);
plot(0:dt:T-dt, (x_true_store(1,:) - x_est_store(1,:)));
ylabel('Position Error x axis [m]'); title('EKF Position Estimation Error');

subplot(3,1,2);
plot(0:dt:T-dt, (x_true_store(2,:) - x_est_store(2,:)));
ylabel('Position Error y axis [m]'); title('EKF Position Estimation Error');

subplot(3,1,3);
plot(0:dt:T-dt, (x_true_store(3,:) - x_est_store(3,:)));
ylabel('Position Error z axis [m]'); title('EKF Position Estimation Error');
sgtitle('Individual position errors')

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
    % ra = atan(pos_topo(2)/pos_topo(1));
    ra = atan2(pos_topo(2), pos_topo(1));
    dec = asin(pos_topo(3)/range);
    h = [ra;dec];
    
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


function dX = int_twobody_stm(t, X, mu)
    r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
    
    % Jacobian of two-body dynamics
    A41 = -mu * (1/r^3 - 3*X(1)^2/r^5);
    A42 = 3*mu/r^5 * X(1) * X(2);
    A43 = 3*mu/r^5 * X(1) * X(3);
    A52 = -mu * (1/r^3 - 3*X(2)^2/r^5);
    A53 = 3*mu/r^5 * X(2) * X(3);
    A63 = -mu * (1/r^3 - 3*X(3)^2/r^5);
    
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         A41 A42 A43 0 0 0;
         A42 A52 A53 0 0 0;
         A43 A53 A63 0 0 0];

    % State + STM derivative
    dX = zeros(42,1);
    dX(1:6) = [X(4); X(5); X(6); -mu*X(1:3)/r^3];
    dX(7:end) = reshape((A * reshape(X(7:end), 6, 6)').', 1, []); 
    % Phi = reshape(X(7:end), 6, 6);
    % dPhi = A * Phi;
    % dX(7:end) = reshape(dPhi, 36, 1);
end

function Gamma = getGammaMatrix(t_pre, t_cur)

    timediff = t_cur - t_pre;
    Gamma = zeros(6, 3);
    for k=1:size(Gamma, 2)
        Gamma(k,k) = 1/2 * timediff^2;
        Gamma(k+3,k) = timediff;
    end
end

% save("ekf_es_true_prop_data.mat","x_true_store","x_est_store","z_obs_store","P_est_store","ra_est_store","ra_true_store","dec_true_store","dec_est_store");