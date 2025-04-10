%% Reactivity analysis for invasive species (MATLAB(R) code)
%
%%

%% Model parameters
% Values for _Ailanthus altissima_ 

% clean-up
clearvars
close all
clc

b3=3110.44; % per-capita yearly seed production for individuals in 
            % development stage 3
v1=1; v2=0.751; % [vX] viability probability for individuals in development 
                % stage X = 1, 2
e1=1.34e-3; e2=1.33e-3;  % [eX] emergence probability for individuals in 
                % development stage X = 1, 2
s2=1; s3=0.76; s4=0.75; % [sX] survival probability for individuals in 
                % development stage X = 2, 3, 4
g2=0.20; g3=1.32e-2; % [gX] growth probability for individuals in 
                % development stage X = 2, 3
r4=0; % retrogression probability for individuals in development stage 4


%% Projection matrix

A=[v2*(1-e2) 0 0 b3*v1*(1-e1);
    v2*e2 0 0 b3*v1*e1; 
    0 s2*g2 s3*(1-g3) s4*r4; 
    0 0 s3*g3 s4*(1-r4)];

display(A)

%% Asymptotic stability and reactivity analysis

% dominant eigenvalue of the projection matrix
lambda=eigs(A,1,'largestabs');
display(lambda)

% net reproduction value
N=b3*v1*(e1-e1*v2+e2*v2)*s2*g2*s3*g3/(1-v2+e2*v2)/...
    (1-s3*(1-g3)-s4*(1-r4)+s3*s4*(1-g3-r4));
display(N)

% first-timestep reactivity
G=b3*v1+s4;
display(G)

%% Controls

% net reproduction value threshold to prevent transient growth
Nstar=(1-s4)/b3/v1*N;
display(Nstar)

% relative reduction in reproduction value to prevent transient growth
RRR=1-(1-s4)/b3/v1;
display(RRR)
