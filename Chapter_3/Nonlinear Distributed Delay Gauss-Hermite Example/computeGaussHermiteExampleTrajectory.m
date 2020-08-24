function sol = computeGaussHermiteExampleTrajectory(mu,sigma,N,t_mesh, t,y, Z)
    y0 = @(t) exp(t);

    %Compute the Gauss Hermite quadrature weights for our given normal
    %distribution, which lie inside the integration limits [a,b].
    [evaluation_points, integrand_weights] = computeWeights(0,2*mu,mu,sigma,N);

    %The derivative function for our example distributed delay equation.
    dydt = @(t,y,Z) exampleDerivative(sigma*sqrt(2*pi).*integrand_weights(:)', 0,2*mu,mu,sigma,t,y,Z);
    
    %Generate an event function to stop the simulation, if it exceeds the
    %specified threshold.
    eventFunction = @(t,y,Z) terminalEventFcn(10e1,t,y,Z);
    options = odeset('RelTol',10e-10,'AbsTol',10e-10);


        %If zero is an evaluation point in our delay_times vector, then exclude
    %it from the dde23 function call and re-insert it in the underlying
    %distributedDelaySchnakenberg function.
    if(evaluation_points(1) == 0)  
        solution = dde23(dydt, double(evaluation_points(2:end)), y0, [0,max(t_mesh)], options);
    else
        solution = dde23(dydt, double(evaluation_points), y0, [0,max(t_mesh)], options);
    end

    sol = solution;
end


%This function is used to compute the Gauss Hermite weights and evaluation
%points and then filter our the values which lie outside of the domain 
%[a,b]
function [evaluation_points, integrand_weights] = computeWeights(a,b,mu,sigma,N)
    [points, weights] = computeGaussHermiteWeights(mu,sigma,N);
    
    
    valid_indicies = (points > a) .* (points < b);
    
    evaluation_points = points .* valid_indicies;
    integrand_weights = weights .* valid_indicies;
    
    evaluation_points = nonzeros(evaluation_points);
    integrand_weights = nonzeros(integrand_weights);
    
end


%This function is used to compute the derivative of the example distributed
%delay equation.
function derivative = exampleDerivative(weights, a, b, mu, sigma, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    
    
    
    integrand = function_values(1,:).^(2.*mu/(sigma.^2));
    integral_approximation = dot(weights, integrand);
    
    
    
    mu1 = (2*t*sigma - mu) / (2*sigma.^2 - 1);
    sigma1 = sigma/(sqrt(1-2*sigma^2));
    
    
    derivative_x = y(1) + integral_approximation ...
                    -sqrt(2*pi)*(normcdf((b+mu)/sigma)-normcdf((a+mu)/sigma));
    
    derivative=double(derivative_x);
    
                
     %derivative_y = y(1) + integral_approximation - y(2).^2 .* ...
     %               exp((((2*(sigma^2) - mu).^2)-mu.^2)/(sigma.^2)) .* ...
     %               (normcdf(b+2*sigma.^2, mu, sigma) - normcdf(a+2*sigma.^2, mu, sigma));
    
     %derivative = [double(derivative_x), double(derivative_y)]';
end


function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end