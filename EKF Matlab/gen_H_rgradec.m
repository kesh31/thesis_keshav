function [Hk_til, Gk] = gen_H_rgradec(t, X, params)
 
    % Compute Earth rotaion angle in [rad]
    theta = params.theta0 + params.dtheta * t;

    % Retrieve position of satellite
    pos = X(1:3);

    % Transform stations position from ECEF to ECI
    if size(params.stat_ecef) == [1, 3]
        stat_ecef = params.stat_ecef';      % ensure to have column vector 
    end
    trans = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
    % trans = ecefToEciSimplified(theta);
    stat_eci =  trans * stat_ecef;
    %stat_eci = [-2958476.164700; 5610449.069874; 669294.973672];

    % Compute expected measurements
    pos_topo = pos-stat_eci;
    range = norm(pos_topo);
    ra = atan(pos_topo(2)/pos_topo(1));
    dec = asin(pos_topo(3)/range);
    Gk = [range; ra; dec];

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

    Hk_til = [H11, H12, H13, 0, 0, 0;
              H21, H22,   0, 0, 0, 0;
              H31, H32, H33, 0, 0, 0];
end