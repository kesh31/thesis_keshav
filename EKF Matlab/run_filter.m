
clear;clc;close all;

% Load input data file
params = load("orbit_model_inputs_radec.mat", 'P0', 'X0_true', ...
              'Rk','GM','Q','Re','dtheta', 'stat_ecef', 'theta0');

truth_data = load('orbit_model_truth.mat','Xt_mat','time');
% Load measurements
meas = load("orbit_model_meas_radec", "tvec", "obs_data");


% Initialize input for batch
% X0_ref = params.X0_true;
X0_ref = [6.06900255e+06,6.09157542e+06,1.05883186e+06,-1.25648342e+03,3.83516911e+03,-5.98898761e+03]'/1000; % Debug

intfcn = @int_twobody_stm;
H_fcn = @gen_H_radec;


[Xref_mat, P_mat, resids] = ekf(X0_ref, meas.tvec, meas.obs_data, ...
                                       intfcn, H_fcn, params,truth_data);

save("orbit_model_ekf_output_radec.mat", ...
       "Xref_mat", "P_mat", "resids");
% [Xref_mat, P_mat, resids] = ckf(X0_ref, meas.tvec, meas.obs_data, ...
%                                         intfcn, H_fcn, params);
% save("orbit_model_ckf_output_radec.mat", ...
%         "Xref_mat", "P_mat", "resids");
