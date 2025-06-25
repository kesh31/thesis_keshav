function dX = int_twobody_stm(t, X, params)

    % Initialization
    GM = params.GM;
    r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);

    A41 = -GM * (1/r^3 - 3*X(1)^2/r^5);
    A42 = 3*GM/r^5 * X(1) * X(2);
    A43 = 3*GM/r^5 * X(1) * X(3);
    A52 = -GM * (1/r^3 - 3*X(2)^2/r^5);
    A53 = 3*GM/r^5 * X(2) * X(3);
    A63 = -GM * (1/r^3 - 3*X(3)^2/r^5);
    
    A = [0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 1, 0;
         0, 0, 0, 0, 0, 1;
         A41, A42, A43, 0, 0, 0;
         A42, A52, A53, 0, 0, 0;
         A43, A53, A63, 0, 0, 0];

    % First six rows store state derivative and last rows derivative of STM
    dX = zeros(42, 1); 
    dX(1:6) = twobody(t, X(1:6), params.GM);
    dX(7:end) = reshape((A * reshape(X(7:end), 6, 6)').', 1, []);   
end


function dx = twobody(~, x, mu)
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    dx = [v; a];
end