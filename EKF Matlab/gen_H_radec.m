function [Hk_til, Gk] = gen_H_radec(t, X, params)

    % This function does exactly the same as gen_H_rgradec() except from
    % returning range as a measurement. Therefore, all related range 
    % component need to be removed from Gk and Hk_til

    [Hk_til_withRange, Gk_withRange] = gen_H_rgradec(t, X, params);
    %[Hk_til_withRange, Gk_withRange] = gen_H_rgradec_orekit(t, X, params);
    Hk_til = Hk_til_withRange(2:end, :);
    Gk = Gk_withRange(2:end);
end