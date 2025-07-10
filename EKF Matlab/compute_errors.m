% function compute_errors()

    close all;
    % Load truth data
    X_truth = load("orbit_model_truth.mat", "time", "Xt_mat");
    time = X_truth.time;

    % Load estimated data
    X_estimated = load("orbit_model_ekf_output_radec.mat", ...
                       "Xref_mat", "P_mat", "resids");
    n = size(X_truth.Xt_mat,1);     % size of state vector

    % Standard deviations for all the positions and velocities in ECI
    sigma_eci = zeros(n, length(X_truth.time));
    for k=1:1:size(X_truth.Xt_mat,1)
        sigma_eci(k,:) = squeeze(sqrt(X_estimated.P_mat(k,k,:)));
    end

    % Compute error
    Xerr_mat = X_estimated.Xref_mat - X_truth.Xt_mat;
    norm_Xerr_mat = vecnorm(Xerr_mat);

    % Compute position error and standard deviation in RIC frame
    err_ric = zeros(3, length(time));
    cov_ric = zeros(3, 3, length(time));
    sigma_ric = zeros(3, length(time));

    for k=1:length(time)
        rc_vec = X_truth.Xt_mat(1:3,k);
        vc_vec = X_truth.Xt_mat(4:6,k);

        % Position error in RIC
        Qin = Xerr_mat(1:3,k);
        % err_ric(:,k) = eci2ric(rc_vec, vc_vec, Qin);

        % Position covariance in RIC
        Qin = X_estimated.P_mat(1:3,1:3,k);
        % cov_ric(:,:,k) = eci2ric(rc_vec, vc_vec, Qin);
    end

    % Position standard deviation in RIC
    % for k=1:3
    %     sigma_ric(k,:) = squeeze(sqrt(cov_ric(k,k,:)));
    % end

    % Conversion from seconds to hours
    time = time/3600;

    % Plot 
    % Position errors in ECI (X,Y,Z) with 3-sigma bounds
    figure(1)
    subplot(3,1,1)
    plot(time, Xerr_mat(1, :) * 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(1, :) * 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(1, :) * 1e3, '-')
    hold on
    % xlim([0 2.8833])
    % ylim([-1e3  1e3])
    ylabel("PosX Error [m]")
    title('Position error in ECI')

    subplot(3,1,2)
    plot(time, Xerr_mat(2,:)* 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(2,:)* 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(2,:)* 1e3, '-')
    % xlim([0 2.8833])
    % ylim([-1e3  1e3])
    ylabel("PosY Error [m]")

    subplot(3,1,3)
    plot(time, Xerr_mat(3,:)* 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(3,:)* 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(3,:)* 1e3, '-')
    %xlim([0 2.8833])
    % ylim([-1e3  1e3])
    % xlabel("time [h]")
    ylabel("PosZ Error [m]")

    figure(2)
    subplot(3,1,1)
    plot(time, Xerr_mat(4,:) * 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(4,:) * 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(4,:) * 1e3, '-')
    % xlim([0 2.8833])
    %ylim([-2e0 2e0])
    ylabel("VelX Error [m/s]")
    title('Velocity error in ECI')

    subplot(3,1,2)
    plot(time, Xerr_mat(5,:)* 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(5,:)* 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(5,:)* 1e3, '-')
    % xlim([0 2.8833])
    %ylim([-2e0 2e0])
    ylabel("VelY Error [m/s]")

    subplot(3,1,3)
    plot(time, Xerr_mat(6,:)* 1e3, '.')
    hold on
    plot(time, 3*sigma_eci(6,:)* 1e3, '-')
    hold on
    plot(time, -3*sigma_eci(6,:)* 1e3, '-')
    % xlim([0 2.8833])
    ylim([-2e1 2e1])
    xlabel("time [h]")
    ylabel("VelZ Error [m/s]")

    % figure(3)
    % subplot(3,1,1)
    % plot(time, err_ric(1,:) * 1e3, '.')
    % hold on
    % plot(time, 3*sigma_ric(1,:) * 1e3, '-')
    % hold on
    % plot(time, -3*sigma_ric(1,:) * 1e3, '-')
    % % xlim([0 2.8833])
    % % ylim([-2e3 2e3])
    % ylabel("Radial Error [m]")
    % title('Position error in RIC')

    % subplot(3,1,2)
    % plot(time, err_ric(2,:)* 1e3, '.')
    % hold on
    % plot(time, 3*sigma_ric(2,:)* 1e3, '-')
    % hold on
    % plot(time, -3*sigma_ric(2,:)* 1e3, '-')
    % % xlim([0 2.8833])
    % % ylim([-2e3 2e3])
    % ylabel("In track Error [m]")
    % 
    % subplot(3,1,3)
    % plot(time, err_ric(3,:)* 1e3, '.')
    % hold on
    % plot(time, 3*sigma_ric(3,:)* 1e3, '-')
    % hold on
    % plot(time, -3*sigma_ric(3,:)* 1e3, '-')
    % % xlim([0 2.8833])
    % % ylim([-2.5e3 2.5e3])
    % xlabel("time [h]")
    % ylabel("Cross Track Error [m]")

    figure(4)
    % Convert units of residuals
    % range = X_estimated.resids(1,:) * 1e3;  % in [m]
    % ra = X_estimated.resids(2,:) * 206265;  % in [arcsec]
    % dec = X_estimated.resids(3,:) * 206265; % in [arcsec]
    ra = X_estimated.resids(1,:) * 206265;  % in [arcsec]
    dec = X_estimated.resids(2,:) * 206265; % in [arcsec]

    subplot(3,1,1)
    % plot(time, range, '.')
    % % xlim([0 1])
    % % ylim([-15 30])
    % ylabel("Range error [m]")
    % title('Residuals in Topocentric Frame')

    subplot(3,1,2)
    plot(time, ra(:), '.')
    % xlim([0 2.5])
    % ylim([-20 20])
    ylabel("RA Error [arcsec]")

    subplot(3,1,3)
    plot(time, dec(:), '.')
    % xlim([0 2.5])
    % ylim([-20 20])
    xlabel("time [h]")
    ylabel("Dec Error [arcsec]")

    figure(5)
    plot(time(), norm_Xerr_mat() * 1e3)
    %plot(time, norm_Xerr_mat)
    %ylim([0 2e5])
    xlabel("time [h]")
    ylabel("Normalized position error in ECI [m]")

    figure(6)
    plot(time, X_estimated.Xref_mat(1,:) * 1e3)
    xlabel("time [h]")
    ylabel("X Position in ECI [m]")
    
% Position errors in ECI (X,Y,Z) with 3-sigma bounds
figure;
subplot(3,1,1)
plot(time(50:end), Xerr_mat(1, 50:end), '.')
hold on
plot(time(50:end), 3*sigma_eci(1, 50:end), '-')
hold on
plot(time(50:end), -3*sigma_eci(1, 50:end), '-')
hold on
% xlim([0 2.8833])
% ylim([-1e3 2.5e3])
ylabel("PosX Error [m]")
title('Position error in ECI')

subplot(3,1,2)
plot(time(50:end), Xerr_mat(2,50:end), '.')
hold on
plot(time(50:end), 3*sigma_eci(2,50:end), '-')
hold on
plot(time(50:end), -3*sigma_eci(2,50:end), '-')
% xlim([0 2.8833])
% ylim([-2.5e3 2e2])
ylabel("PosY Error [m]")

subplot(3,1,3)
plot(time(50:end), Xerr_mat(3,50:end), '.')
hold on
plot(time(50:end), 3*sigma_eci(3,50:end), '-')
hold on
plot(time(50:end), -3*sigma_eci(3,50:end), '-')
%xlim([0 2.8833])
% ylim([-1.5e3 1.5e3])
% xlabel("time [h]")
ylabel("PosZ Error [m]")