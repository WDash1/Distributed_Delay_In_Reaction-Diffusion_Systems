t_values = linspace(10,25,500);

y0 = @(t) sin(sqrt(2)*t)+cos(t);



alpha = -2.3;
beta =3.45;


N = 11;    

[gauss_legendre_zeros, gauss_legendre_weights] = computeGaussLegendreWeights(0, 1, N);
[~, legendre_soly] = compute_trajectory_simulation(gauss_legendre_zeros, ...
                                                    gauss_legendre_weights, ...
                                                    alpha, ...
                                                    beta, ...
                                                    t_values, ...
                                                    y0);
    
figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
hold on;
box on;
xlim([min(t_values), max(t_values)]);

plot(t_values, legendre_soly, '-', 'DisplayName', 'Gauss Legendre');
xlabel('{\it t}', 'FontSize', 20);
ylabel('{\it y(t)}', 'FontSize', 20);

%legend
ax = gca;
ax.FontSize = 20; 

filename = "Distributed_Delay_Example_Linear_Trajectory_alpha="+string(alpha)+"_beta="+string(beta)+".eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);


function [x,y] = compute_trajectory_simulation(delay_times,weights, alpha, beta, t_values, y0)
    delays = delay_times;
    if(delay_times(1) == 0)
        delays = delay_times(2:end);
    end
    dydt = @(t,y,Z) linearDerivativeExample(weights,alpha, beta, t,y,Z);
    options = odeset('RelTol',10e-10,'AbsTol',10e-10);
    sol = dde23(dydt, delays, y0, [0, max(t_values)],options);
    
    x = sol.x;
    y = deval(sol, t_values);
    
end



function derivative = linearDerivativeExample(weights, alpha, beta, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    integrand = function_values;
    integral_approximation = dot(weights, integrand);
    derivative = alpha*y + beta*integral_approximation;
end