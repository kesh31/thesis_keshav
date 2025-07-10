function [Xref_mat, P_mat, resids] = ekf(X0_ref, t_obs, obs_data, ...
                                         intfcn, H_fcn, params,true_data)
    
    % Initialization
    n = size(X0_ref, 1); 
    x0_bar = zeros(n, 1);
    xhat_pre = x0_bar;
    P_pre = params.P0;
    t_pre = t_obs(1);
    converged = 1;

    % Output initialization
    Xref_mat = zeros(size(X0_ref, 1), length(t_obs));
    P_mat = zeros(size(X0_ref, 1), size(X0_ref, 1), length(t_obs));
    resids = zeros(size(obs_data));
    %xhatk = zeros(1, length(t_obs));

    % Set up settings for ODE45 solver 
    RelTol = 1e-12;                                 % relative tolerance
    AbsTol = 1e-12;                                 % absolute tolerance
    Phi0 = reshape((eye(n)), 1, [])'; 
    Xref_Stm0 = [X0_ref; Phi0];                     % STM initial condition

    % Loop through all measurements
    for k=1:1:length(t_obs)

        if k== 174
            t_obs(k)
        end

        % Update
            tk = t_obs(k);
            % Yk = obs_data(:,k);
            true_state = true_data.Xt_mat(:,k);            % Debug
            [ra_true, dec_true] = measure_debug(true_state,params,tk);% Debug
            Yk = [ra_true; dec_true] + sqrt(params.Rk(2,2)) * randn(2, 1); % Add noise
            Rk = params.Rk;

            % Compute STM and Xref
            options = odeset('RelTol', RelTol,'AbsTol', AbsTol); 

            % ODE45 fails for first measurements
            if k == 1
                Xref_stm = Xref_Stm0';
            else
                [~,Xref_stm] = ode45(@(t, X) intfcn(t, X, params), ...
                                        [t_pre tk], Xref_stm_prev, ...
                                        options);
            end
            Xref_k = Xref_stm(end, 1:n)';
            
            % Predict state deviation vector and covariance
            Phik = reshape(Xref_stm(end, n+1:end), n, n)';
            Gammak = getGammaMatrix(t_pre, tk);

            xk_bar = Phik * xhat_pre;
            Pk_bar = Phik * P_pre * Phik' + Gammak * params.Q * Gammak';

            % Compute expected measurement 
            [Hk_til, Gk] = feval(H_fcn, tk, Xref_k, params);
            yk = Yk - Gk;                           % Define as residuals

            % Corrector     
            Kk = Pk_bar * Hk_til' * inv(Hk_til * Pk_bar * Hk_til' + Rk); 
            xhat = xk_bar + Kk * (yk - Hk_til * xk_bar);
   
            Xref_mat(:,k)  = Xref_k + xhat;
            Pk = (eye(n) - Kk * Hk_til) * Pk_bar ...
                * (eye(n) - Kk * Hk_til)' + Kk * Rk * Kk';
            P_mat(:,:,k) = Pk;
            resids(:,k) = yk - Hk_til * xhat;

            % Extended KF update
            %xhatk(k) = vecnorm(xhat);
            if k~=1 && vecnorm(xhat) > converged
             
                Xref_k = Xref_mat(:,k);
                xhat = zeros(n, 1);
            end

            % Update
            Xref_stm_prev = [Xref_k; Phi0]';
            t_pre = tk;
            P_pre = Pk;
            xhat_pre = xhat;
    end
end

function Gamma = getGammaMatrix(t_pre, t_cur)

    timediff = t_cur - t_pre;
    Gamma = zeros(6, 3);
    for k=1:size(Gamma, 2)
        Gamma(k,k) = 1/2 * timediff^2;
        Gamma(k+3,k) = timediff;
    end

end