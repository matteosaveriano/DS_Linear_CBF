% Constrained DS motion using control barrier functions
% 
% Author: Matteo Saveriano
% Date: 27.06.24

close all
clear

addpath(genpath('./Data'));
addpath(genpath('./Lib'));

light_blue = [0.3010,0.7450,0.9330];

%% Define linear constraints as CBFs
xLim.max = [0.8 0.8];
xLim.min = [-0.8 -.8];

xGoal = [0;0];
dt = 0.01;
x =  [0.7; 0.7];
xOrig = x;
lambda = diag([5,30]);

constrCoef = -[ 1  0 -xLim.max(1);...
               -1  0  xLim.min(1);...
                0  1 -xLim.max(2);...
                0 -1  xLim.min(2);...
               -1  0.3 -0.4 ];           
                      
%% Define stable nonlinear dynamics (SEDS)
SEDS_model = load('CShape.mat');

dsHandle = @(x)GMR(SEDS_model.Priors, SEDS_model.Mu, SEDS_model.Sigma, x, 1:2, 3:4); % @(x)globally_stable_DS(x); %@(x)limit_cycle_DS(x); % % %@(x)2*(xGoal - x); %% 
D = [xLim.min(1) xLim.max(1) xLim.min(2) xLim.max(2)];

%% Define CBF function handle
dsHandleZeroBarrier = @(x)constrained_ds_2D(dsHandle, x, xGoal, constrCoef, 'zcbf');

%% Simulation loop
i = 1;
L = 2000;
u = [];
runningTimes = [];
while i<=L 
    % Constrained DS motion
    xd = dsHandleZeroBarrier(x(:,i));
    x(:,i+1) = x(:,i) + xd*dt;
    
    % Original position (for comparison)
    xdOrig = dsHandle(xOrig(:,i)-xGoal);
    xOrig(:,i+1) = xOrig(:,i) + xdOrig*dt;
       
    % Check for convergence
    i = i + 1;
    if(norm(x(:,i)-xGoal)<1e-3)
        i = L+1;
    end
end

%% Plot Results
figure(1)
plot(x(1,:),x(2,:),'LineWidth',2,'Color',light_blue)
hold on
plot(x(1,1),x(2,1),'o','Color',light_blue)
plot(x(1,end),x(2,end),'*','Color',light_blue)
plot(xOrig(1,:),xOrig(2,:),'r')
plot_linear_constraints(constrCoef, xLim, 'k')
grid on;
axis([-1.2,1.2,-1.2,1.2])