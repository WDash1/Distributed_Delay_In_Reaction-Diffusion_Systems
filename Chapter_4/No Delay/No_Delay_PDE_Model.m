% The parameter settings to be used for simulations of the Schnakenberg
% reaction-diffusion system.
a = 0.2;
b = 1.3;

%The diffusion coefficients for simulations of the PDE model.
Du = 1;
Dv = 100;

%The fixed point for our model.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The number of points which we wish to divide the space interval into.
m = 50;


%The time mesh for our model.
t_mesh = linspace(0,300,500);


%The maximum and minimum spatial values that we wish to use for
%simulations.
x_max = 50;
x_min = 0;

%The discretised space values that will be used in simulations of our PDE.
x_values = (x_min + ((0:m)/m) .*(x_max-x_min))';

%The space step size for our numerical simulations.
h = (x_max - x_min)/m;

%The initial data at t=0 for our simulation.
u0 = @(x) u_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
v0 = @(x) v_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
%initial_data = [u0(x_values) ; v0(x_values)];


u0 = u_fixed+0.001.*rand(1,m+1);
v0 = v_fixed+0.001.*rand(1,m+1);
initial_data = [u0,v0];
initial_data =  [1.5001, ...
    1.5003, ...
    1.5005, ...
    1.5007, ...
    1.5004, ...
    1.5008, ...
    1.5007, ...
    1.5010, ...
    1.5005, ...
    1.5003, ...
    1.5001, ...
    1.5006, ...
    1.5008, ...
    1.5004, ...
    1.5001, ...
    1.5003, ...
    1.5002, ...
    1.5003, ...
    1.5004, ...
    1.5005, ...
    1.5005, ...
    1.5009, ...
    1.5005, ...
    1.5009, ...
    1.5006, ...
    1.5010, ...
    1.5002, ...
    1.5007, ...
    1.5003, ...
    1.5007, ...
    1.5007, ...
    1.5001, ...
    1.5003, ...
    1.5002, ...
    1.5007, ...
    1.5008, ...
    1.5003, ...
    1.5008, ...
    1.5007, ...
    1.5000, ...
    1.5006, ...
    1.5004, ...
    1.5009, ...
    1.5000, ...
    1.5005, ...
    1.5004, ...
    1.5005, ...
    1.5008, ...
    1.5003, ...
    1.5008, ...
    1.5005, ...
    0.5778, ...
    0.5780, ...
    0.5785, ...
    0.5783, ...
    0.5779, ...
    0.5781, ...
    0.5784, ...
    0.5780, ...
    0.5785, ...
    0.5780, ...
    0.5787, ...
    0.5780, ...
    0.5785, ...
    0.5780, ...
    0.5781, ...
    0.5779, ...
    0.5784, ...
    0.5785, ...
    0.5783, ...
    0.5782, ...
    0.5784, ...
    0.5784, ...
    0.5785, ...
    0.5784, ...
    0.5787, ...
    0.5780, ...
    0.5785, ...
    0.5780, ...
    0.5779, ...
    0.5784, ...
    0.5782, ...
    0.5782, ...
    0.5784, ...
    0.5785, ...
    0.5781, ...
    0.5784, ...
    0.5782, ...
    0.5786, ...
    0.5786, ...
    0.5780, ...
    0.5784, ...
    0.5784, ...
    0.5783, ...
    0.5786, ...
    0.5780, ...
    0.5781, ...
    0.5779, ...
    0.5787, ...
    0.5784, ...
    0.5783, ...
    0.5784];
%The derivative function for our model.
dydt = @(t,y) computeNoDelaySchnakenbergDerivative(a,b,Du,Dv,m,h,y);



%Simulate our PDE and extract the u and v values.
[time_mesh,solution_mesh] = ode15s(dydt,t_mesh,initial_data);
u_values = solution_mesh(:,1:m+1);
v_values = solution_mesh(:,m+1+(1:(m+1)));
[XX,YY] = meshgrid(x_values,t_mesh);

%Plot the u values from the simulation.
figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
h1=surf(XX,YY, u_values);
view([0  90])
axis square
grid off
box on
set(h1, 'EdgeColor','none')
xlabel('{\it x}');
ylabel('{\it t}');
set(gca,'FontSize',20)
title('u values');
xlim([min(x_values), max(x_values)])
ylim([min(t_mesh), max(t_mesh)])

set(gca,'XTick',min(x_values):10:max(x_values));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    
set(gca,'YTick',min(t_mesh):60:max(t_mesh));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))
colorbar
filename1 = 'Output_Images/No_Delay_PDE_Model_u_values_a='+string(a)+'_b='+string(b)+'_m='+string(m)+'_t='+max(t_mesh);
print(filename1+'.png', '-dpng', '-r300');
savefig(filename1+'.fig');

%Plot the u values from the simulation.
figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
h1=surf(XX,YY, v_values);
view([0  90])
axis square
grid off
box on
set(h1, 'EdgeColor','none')
xlabel('{\it x}');
ylabel('{\it t}');
set(gca,'FontSize',20)
title('v values');
xlim([min(x_values), max(x_values)])
ylim([min(t_mesh), max(t_mesh)])

set(gca,'XTick',min(x_values):10:max(x_values));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    
set(gca,'YTick',min(t_mesh):60:max(t_mesh));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))
colorbar
filename2 = 'Output_Images/No_Delay_PDE_Model_v_values_a='+string(a)+'_b='+string(b)+'_m='+string(m)+'_t='+max(t_mesh);
print(filename2+'.png', '-dpng', '-r300');
savefig(filename2+'.fig');

