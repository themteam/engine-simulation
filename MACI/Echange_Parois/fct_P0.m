function [ theta, M ] = fct_P0( cond_init, ths )
%FCT_P0 Summary of this function goes here
%   Detailed explanation goes here
    
    options = odeset('Mass', 'M(t,y)', 'RelTol', 1e-3, ...
                     'AbsTol', [1e2, 1e-1, 1e-7, 1e-7, 1e-7, 1e-2, 1e-7]);
    
    % Initialisation
    y0 = [cond_init.P0, cond_init.T0, cond_init.m0, ...
          cond_init.mu0, cond_init.mb0, cond_init.f0, ...
          cond_init.mcapa0];
    [theta, M] = ode45('systemeFunction1', 0:0.5:720, y0, options);
    while true
        % Sauvegarde
        M_prev = M;
        % Itération
        y0 = M(end, :);
        [theta, M] = ode45('systemeFunction1', 0:0.5:720, y0, options);
        % Vérification de satisfaction du critère d'arrêt
        if ths < norm(M - M_prev, 'fro')
            break;
        end
    end
    
end

