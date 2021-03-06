t_values = linspace(0,2,500);

y0 = @(t) [sin(t),cos(t)];


n_max = 10;
n_values = 1:n_max;

trapezium_rule_cauchy_sequence = zeros(size(n_values));
simpsons_rule_cauchy_sequence = zeros(size(n_values));
simpsons_38_rule_cauchy_sequence = zeros(size(n_values));
gauss_legendre_cauchy_sequence = zeros(size(n_values));

trapezium_current_trajectory = zeros(2,size(t_values,2));
simpsons_current_trajectory = zeros(2,size(t_values,2));
simpsons_38_current_trajectory = zeros(2,size(t_values,2));
gauss_legendre_current_trajectory = zeros(2,size(t_values,2));
    
    
for n = n_values
    
    spmd(4)
        if labindex == 1
            [trapezium_rule_zeros, trapezium_rule_weights] = computeTrapeziumRuleWeights(0, 1, 6*n+1);
            [~, trapezium_soly] = compute_trajectory_simulation(trapezium_rule_zeros, trapezium_rule_weights, t_values, y0);
            if(n>1)
                trapezium_rule_cauchy_sequence(n-1) = max(vecnorm(trapezium_soly - trapezium_current_trajectory));
            end
            trapezium_current_trajectory = trapezium_soly;
        elseif labindex == 2
            [simpsons_rule_zeros, simpsons_rule_weights] = computeCompositeSimpsonRuleWeights(0, 1, 6*n+1);
            [~, simpsons_soly] = compute_trajectory_simulation(simpsons_rule_zeros, simpsons_rule_weights, t_values, y0);
            if(n>1)
                simpsons_rule_cauchy_sequence(n-1) = max(vecnorm(simpsons_soly - simpsons_current_trajectory));
            end
            simpsons_current_trajectory = simpsons_soly;
        elseif labindex == 3
            [simpsons_38_rule_zeros, simpsons_38_rule_weights] = computeCompositeSimpson38RuleWeights(0, 1, 6*n+1);
            [~, simpsons_38_soly] = compute_trajectory_simulation(simpsons_38_rule_zeros, simpsons_38_rule_weights, t_values, y0);
            if(n>1)
                simpsons_38_rule_cauchy_sequence(n-1) = max(vecnorm(simpsons_38_soly - simpsons_38_current_trajectory));
            end
            simpsons_38_current_trajectory = simpsons_38_soly;
        elseif labindex == 4
            [gauss_legendre_zeros, gauss_legendre_weights] = computeGaussLegendreWeights(0, 1, 6*n+1);
            [~, legendre_soly] = compute_trajectory_simulation(gauss_legendre_zeros, gauss_legendre_weights, t_values, y0);
            if(n>1)
                gauss_legendre_cauchy_sequence(n-1) = max(vecnorm(legendre_soly - gauss_legendre_current_trajectory));
            end
            gauss_legendre_current_trajectory = legendre_soly;
        end
    end
end
figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
hold on;
box on;
xlim([min(6.*n_values+1), max(6.*(n_max-1)+1)]);
loglog(6.*n_values+1, trapezium_rule_cauchy_sequence{1}, '-o', 'DisplayName', 'Trapezium Rule');
loglog(6.*n_values+1, simpsons_rule_cauchy_sequence{2}, '-o', 'DisplayName', "Simpson's Rule");
loglog(6.*n_values+1, simpsons_38_rule_cauchy_sequence{3}, '-o', 'DisplayName', "Simpson's 3/8 Rule");
loglog(6.*n_values+1, gauss_legendre_cauchy_sequence{4}, '-o', 'DisplayName', "Gauss Legendre");
xlabel('{\it N}', 'FontSize', 20);
%legend
set(gca,'XScale', 'log', 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
filename = "nonlinear_example2_cauchy_curves.eps";
print -depsc -tiff -r300 -painters filename;


function [x,y] = compute_trajectory_simulation(delay_times,weights, t_values, y0)
    delays = delay_times;
    if(delay_times(1) == 0)
        delays = delay_times(2:end);
    end
    dydt = @(t,y,Z) linearDerivativeExample(delay_times,weights, t,y,Z);
    options = odeset('RelTol',10e-10,'AbsTol',10e-10);
    sol = dde23(dydt, delays, y0, t_values,options);
    
    x = sol.x;
    y = deval(sol, t_values);
    
end


function derivative = linearDerivativeExample(delay_times, weights, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    
    integrand = cos(delay_times(:)') .* sqrt(1-((function_values(1,:)).^2)).* sign(cos(t-delay_times(:)'));
    integral_approximation = dot(weights(:), integrand);
    
    derivative_x = y(2);
    derivative_y = -y(1) - y(2) * (sin(2)/4 +1/2) - y(1)*(1-cos(2))/4 + integral_approximation;
    derivative = [derivative_x, derivative_y]';
end