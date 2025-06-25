clear;clc;close all;

% Checking true state simulated by me vs data from linh's mat file

full_sim_data = load("ekf_es_true_prop_data.mat","x_true_store","x_est_store","z_obs_store","P_est_store","ra_est_store","ra_true_store","dec_true_store","dec_est_store");

meas_data = load("orbit_model_meas_radec", "tvec", "obs_data");

true_data = load("orbit_model_truth.mat",'Xt_mat','time');

time = true_data.time;
n = 6;
sigma_eci = zeros(n, length(true_data.time));
for k=1:6
    for k2 = 1:length(true_data.time)
        sigma_eci(k,k2) = full_sim_data.P_est_store(k,(k2-1)*6+1);
    end
end
Xerr_mat = full_sim_data.x_est_store - full_sim_data.x_true_store;
norm_Xerr_mat = vecnorm(Xerr_mat);

ra_err = full_sim_data.ra_est_store - full_sim_data.ra_true_store;
dec_err = full_sim_data.dec_est_store - full_sim_data.dec_true_store;

% Position errors in ECI (X,Y,Z) with 3-sigma bounds
figure;
subplot(3,1,1)
plot(time, Xerr_mat(1, :), '.')
hold on
plot(time, 3*sigma_eci(1, :), '-')
hold on
plot(time, -3*sigma_eci(1, :), '-')
hold on
% xlim([0 2.8833])
% ylim([-1e3 2.5e3])
ylabel("PosX Error [m]")
title('Position error in ECI')

subplot(3,1,2)
plot(time, Xerr_mat(2,:), '.')
hold on
plot(time, 3*sigma_eci(2,:), '-')
hold on
plot(time, -3*sigma_eci(2,:), '-')
% xlim([0 2.8833])
% ylim([-2.5e3 2e2])
ylabel("PosY Error [m]")

subplot(3,1,3)
plot(time, Xerr_mat(3,:), '.')
hold on
plot(time, 3*sigma_eci(3,:), '-')
hold on
plot(time, -3*sigma_eci(3,:), '-')
%xlim([0 2.8833])
% ylim([-1.5e3 1.5e3])
% xlabel("time [h]")
ylabel("PosZ Error [m]")

figure;
subplot(3,1,1)
plot(time, Xerr_mat(4,:), '.')
hold on
plot(time, 3*sigma_eci(4,:), '-')
hold on
plot(time, -3*sigma_eci(4,:), '-')
% xlim([0 2.8833])
%ylim([-2e0 2e0])
ylabel("VelX Error [m/s]")
title('Velocity error in ECI')

subplot(3,1,2)
plot(time, Xerr_mat(5,:), '.')
hold on
plot(time, 3*sigma_eci(5,:), '-')
hold on
plot(time, -3*sigma_eci(5,:), '-')
% xlim([0 2.8833])
%ylim([-2e0 2e0])
ylabel("VelY Error [m/s]")

subplot(3,1,3)
plot(time, Xerr_mat(6,:), '.')
hold on
plot(time, 3*sigma_eci(6,:), '-')
hold on
plot(time, -3*sigma_eci(6,:), '-')
% xlim([0 2.8833])
% ylim([-2e1 2e1])
xlabel("time [h]")
ylabel("VelZ Error [m/s]")

figure;
subplot(2,1,1)
plot(time, ra_err(:)*206265, '.')
% xlim([0 2.5])
% ylim([-20 20])
ylabel("RA Error [arcsec]")

subplot(2,1,2)
plot(time, dec_err(:)*206265, '.')
% xlim([0 2.5])
% ylim([-20 20])
xlabel("time [h]")
ylabel("Dec Error [arcsec]")

figure;
plot(time, norm_Xerr_mat)
%plot(time, norm_Xerr_mat)
%ylim([0 2e5])
xlabel("time [h]")
ylabel("Normalized position error in ECI [m]")

% figure;
% subplot(3,1,1)
% plot(full_sim_data.x_true_store(1,:));
% hold on
% plot(true_data.Xt_mat(1,:)*1000);
% subplot(3,1,2)
% plot(full_sim_data.x_true_store(2,:));
% hold on
% plot(true_data.Xt_mat(2,:)*1000);
% subplot(3,1,3)
% plot(full_sim_data.x_true_store(3,:));
% hold on
% plot(true_data.Xt_mat(3,:)*1000);
% sgtitle('Comparing true state data')
%
% figure;
% subplot(3,1,1)
% plot(full_sim_data.x_true_store(1,:)-true_data.Xt_mat(1,:)*1000);
% subplot(3,1,2)
% plot(full_sim_data.x_true_store(2,:) - true_data.Xt_mat(2,:)*1000);
% subplot(3,1,3)
% plot(full_sim_data.x_true_store(3,:)-true_data.Xt_mat(3,:)*1000);
%
% sgtitle('Error in true position data')
%
% figure;
% subplot(3,1,1)
% plot(full_sim_data.x_true_store(4,:)-true_data.Xt_mat(4,:)*1000);
% subplot(3,1,2)
% plot(full_sim_data.x_true_store(5,:) - true_data.Xt_mat(5,:)*1000);
% subplot(3,1,3)
% plot(full_sim_data.x_true_store(6,:)-true_data.Xt_mat(6,:)*1000);
%
% sgtitle('Error in true velocity data')
