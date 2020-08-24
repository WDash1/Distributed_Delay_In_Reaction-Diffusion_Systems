t_values = linspace(0,0.5,500);

y0 = @(t) exp(t);


mu = 5;
%sigma = mu*0.1;

n_amt = 10;
n_values = 2.*(1:n_amt)+1;

sigma_amt = 3;
sigma_values = [mu*0.01, mu*0.1, mu*0.05];

gauss_hermite_error_sequence = zeros(sigma_amt, n_amt);


for sigma_index = 1:sigma_amt
    for n_index = 1:n_amt
        n = n_values(n_index);
        sigma = sigma_values(sigma_index);

        sol = computeGaussHermiteExampleTrajectory(mu,sigma,n,t_values);
                                              
       

        hermite_soly = deval(sol, t_values);
        gauss_hermite_error_sequence(sigma_index,n_index) = max(abs(hermite_soly - y0(t_values)));
        
        
    end
end


figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
hold on;
box on;
xlim([min(n_values), max(n_values)]);
for sigma_index = 1:sigma_amt
    loglog(n_values, gauss_hermite_error_sequence(sigma_index,:), '-o', 'DisplayName', '\sigma='+string(sigma_values(sigma_index)));
end
xlabel('{\it N}', 'FontSize', 20);
%legend
set(gca,'XScale', 'log', 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
filename = "LI_Distributed_Delay_Cauchy_Curve23_a="+string(a)+"_b="+string(b)+"_tau="+string(mu)+"_sigma="+string(sigma);
print('-depsc', '-tiff', '-r300', '-painters', filename+".eps");

function [x,y] = compute_trajectory_simulation(delay_times,weights, a,b,mu,sigma,t_values, u0,v0)
    dist = makedist('Normal', mu,sigma);
    trunc = truncate(dist, 0, 2*mu);
    truncated_gaussian_coefficients = pdf(trunc, delay_times);
    weights2 = weights.*truncated_gaussian_coefficients;

    sol = computeDistributedDelayLISchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                delay_times, ...
                                                weights2, ...
                                                max(t_values), ...
                                                u0, ...
                                                v0, ...
                                                10e-10, ...
                                                10e10);

    x = sol.x;
    y = deval(sol, t_values);
    
end