function [ra, dec] = measure_debug(x_sat, params,tk)

    % Compute Earth rotaion angle in [rad]
    theta = params.theta0 + params.dtheta * tk;

    % Transform stations position from ECEF to ECI
    if size(params.stat_ecef) == [1, 3]
        stat_ecef = params.stat_ecef';      % ensure to have column vector 
    end
        trans = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
    % trans = ecefToEciSimplified(theta);
    gs_eci =  trans * stat_ecef;
    %stat_eci = [-2958476.164700; 5610449.069874; 669294.973672];

    r_topo = x_sat(1:3) - gs_eci;
    r = norm(r_topo);
    ra = atan2(r_topo(2), r_topo(1));
    dec = asin(r_topo(3)/r);
end