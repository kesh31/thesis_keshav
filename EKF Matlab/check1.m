% Checking true state simulated by me vs data from linh's mat file

full_sim_data = load("true_prop_sim_meas_data.mat","x_true_store","x_est_store","z_obs_store","P_est_store","ra_est_store","ra_true_store","dec_true_store","dec_est_store");

meas_data = load("orbit_model_meas_radec", "tvec", "obs_data");

true_data = load("orbit_model_truth.mat",'Xt_mat','time');

figure;
subplot(3,1,1)
plot(full_sim_data.x_true_store(1,:));
hold on
plot(true_data.Xt_mat(1,:)*1000);
subplot(3,1,2)
plot(full_sim_data.x_true_store(2,:));
hold on
plot(true_data.Xt_mat(2,:)*1000);
subplot(3,1,3)
plot(full_sim_data.x_true_store(3,:));
hold on
plot(true_data.Xt_mat(3,:)*1000);
sgtitle('Comparing true state data')

figure;
subplot(3,1,1)
plot(full_sim_data.x_true_store(1,:)-true_data.Xt_mat(1,:)*1000);
subplot(3,1,2)
plot(full_sim_data.x_true_store(2,:) - true_data.Xt_mat(2,:)*1000);
subplot(3,1,3)
plot(full_sim_data.x_true_store(3,:)-true_data.Xt_mat(3,:)*1000);

sgtitle('Error in true position data')

figure;
subplot(3,1,1)
plot(full_sim_data.x_true_store(4,:)-true_data.Xt_mat(4,:)*1000);
subplot(3,1,2)
plot(full_sim_data.x_true_store(5,:) - true_data.Xt_mat(5,:)*1000);
subplot(3,1,3)
plot(full_sim_data.x_true_store(6,:)-true_data.Xt_mat(6,:)*1000);

sgtitle('Error in true velocity data')
