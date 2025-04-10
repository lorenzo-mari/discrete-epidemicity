%% Epidemicity analysis for viral respiratory infections (MATLAB(R) code)
%
%%

%% Model parameters
% Values for SARS-CoV-2

% clean-up
clearvars
close all
clc

AV=5.2; % average generation time (days)
SD=1.7; % standard deviation of generation time distribution (days)
gamma=1/11.7; % inverse of infection duration (1/days)
R_0=3.6; % basic reproduction number
R_c=0.95; % control reproduction number

% parameter checks
flag_SD=false;
if SD>AV
    warning(['The standard deviation of the generation time ', ...
        'distribution should not be longer than the average ', ...
        'generation time'])
    flag_SD=true;
end

flag_gamma=false;
if 1/gamma<AV
    warning(['The average duration of infection should not be shorter ',...
        'than the average generation time'])
    flag_gamma=true;
end

flag_R=false;
if R_c>R_0
   warning(['The control reproduction number should not be higher ',...
       'than the basic reproduction number'])
   flag_R=true;
end

if flag_SD || flag_gamma || flag_R
    return
end

%% Numerical discretization
% Daily time step

% finite infectivity
[f_j,q]=eval_f(AV,SD,gamma); 

% probability of still being infected 0.5 days post-exposure
sigma_0=exp(-gamma/2); 

% probability of still being infected from one day to another
sigma=exp(-gamma); 

%% Analytical epidemicity analysis 
% L1-norm

% first-timestep metrics (epidemicity index, reproduction number threshold, 
% herd epidemicity threshold
[E,R_star,HET]=analytics(R_c,max(f_j),sigma_0,sigma,R_0);
display(E)
display(R_star)
display(HET) 

% herd immunity threshold
HIT=max(0,1-1/R_0);
display(HIT) 

%% Numerical epidemicity analysis
% L1-norm

% Leslie matrix
L=build_leslie(q,R_c,f_j,sigma_0,sigma); 

% amplification envelope
timespan=150; % time window to evaluate transient dynamics (days)
A=eval_ampenv(L,timespan); 

% maximum amplification overall and time to peak amplification
[A_max,k_max]=max(A(2:end)); 
display(A_max)
display(k_max)

% maximum duration of the outbreak
k_end=find(A<1,1,'first')-1; 
display(k_end)

% amplification envelope check
if isempty(k_end)
   warning('Increase the simulation timespan')
   return
end

%% Numerical simulations

% most-amplified perturbation at time 1
u_0=zeros(q,1);
[~,i]=max(sum(L));
u_0(i)=1; 

% most-amplified perturbation overall
v_0=zeros(q,1);
[~,i]=max(sum(L^k_max));
v_0(i)=1; 

% initializations
I_u=zeros(q,timespan+1); 
I_u(:,1)=u_0; % for most-amplified perturbation at time 1
I_v=zeros(q,timespan+1);
I_v(:,1)=v_0; % for most-amplified perturbation overall

% simulations
for i=1:timespan
    I_u(:,i+1)=L*I_u(:,i);
    I_v(:,i+1)=L*I_v(:,i);
end

% normalized prevalence (most-amplified perturbation at time 1)
P_u=sum(I_u); 

% normalized prevalence (most-amplified perturbation overall)
P_v=sum(I_v); 

% plotting the simulation results
figure
hold on
box on
cols=lines(2);

plot(0:timespan,P_u,'Color',cols(1,:),...
    'DisplayName','Most amplified perturbation after 1 day')
plot(0:timespan,P_v,'Color',cols(1,:),'LineStyle','--',...
    'DisplayName','Most amplified perturbation overall')
plot(0:timespan,A,'Color',[cols(2,:) 0.5],'LineWidth',2,...
    'DisplayName','Amplification envelope')
xlabel('Time k (days)')
ylabel('Normalized prevalence P(k)')
legend('show','Location','best')

%% Local functions

function [f_j,q]=eval_f(AV,SD,gamma)
    % EVAL_F - Evaluate finite infectivity (daily discretization)
    
    % Gamma-distribution paramenters
    kappa=AV^2/SD^2; % shape parameter
    theta=SD^2/AV; % scale parameter

    % age cutoff for discretization 
    cutoff=1e-4; % numerical threshold
    q1=-1/gamma*log(cutoff);
    q2=fzero(@(x) gamcdf(x,kappa,theta)-(1-cutoff),100);
    q=ceil(max([q1,q2])); % maximum age of infection

    % finite (daily) infectivity
    f_j=NaN(1,q); 
    for j=1:q
        f_j(j)=integral(@(tau) gampdf(tau,kappa,theta)./exp(-gamma*tau),...
            j-1,j); % trapezoid formula
    end
end

function [E,R_star,HET]=analytics(R,f_max,sigma_0,sigma,R_0)
    % ANALYTICS - Analytical epidmicity analysis (L1-norm)
    
    % epidemicity index
    E=sigma_0*R*f_max+sigma;

    % threshold reproduction number to avoid epidemicity
    R_star=1/sigma_0*(1-sigma)/f_max;

    % herd epidemicity threshold
    HET=max(0,1-R_star/R_0);
end

function L=build_leslie(q,R,f,sigma0,sigma)
    % BUILD_LESLIE - Build the epidemiological Leslie projection matrix
    
    L=zeros(q,q);
    L(1,:)=R*sigma0*f;
    L(2:q,1:q-1)=eye(q-1)*sigma;
end

function A=eval_ampenv(L,timespan)
    % EVAL_AMPENV - Evaluate the amplification envelope (L1-norm)
    
    A=NaN(1,timespan+1);
    A(1)=1;

    for j=1:timespan
        A(j+1)=norm(L^j,1);
    end
end