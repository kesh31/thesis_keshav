
clear;clc;close all;

% Load input data file
params = load("orbit_model_inputs_radec.mat", 'P0', 'X0_true', ...
              'Rk','GM','Q','Re','dtheta', 'stat_ecef', 'theta0');

% Load measurements
meas = load("orbit_model_meas_radec", "tvec", "obs_data");

% Load old perturbed initial value
%init = load("perturbed_X0_ref_good", "X0_ref");
%X0_ref = init.X0_ref;

% Initialize input for batch
X0_ref = params.X0_true;
% params.Q = diag([0, 0, 0]);
% Po = params.P0;
% random = randn(6,1);
% pert_vect = sqrt(diag(Po)) .* random;
% X0_ref = params.X0_true + pert_vect;

%intfcn = @int_twobody_sigmaPoints;
%H_fcn = @gen_rgradec_sigmaPoints;
intfcn = @int_twobody_stm;
H_fcn = @gen_H_radec;
[Xref_mat, P_mat, resids] = ekf(X0_ref, meas.tvec, meas.obs_data, ...
                                       intfcn, H_fcn, params);

save("orbit_model_ekf_output_radec.mat", ...
       "Xref_mat", "P_mat", "resids");
% [Xref_mat, P_mat, resids] = ckf(X0_ref, meas.tvec, meas.obs_data, ...
%                                         intfcn, H_fcn, params);
% save("orbit_model_ckf_output_radec.mat", ...
%         "Xref_mat", "P_mat", "resids");
