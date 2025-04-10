%% Epidemicity analysis for tick-borne zoonoses (MATLAB(R) code)
%
%%

%% Model parameters
% Values for _Borrelia burgdorferi_ transmitted by _Ixodes ricinus_

% clean-up
clearvars
close all
clc

% pathogen-related parameters
v=615; % duration of systemic infection (days) 
zLL=0.5; zLN=0.5; zLA=0.5; % [zLX] transmission efficiency from ticks in 
                           % larval stage to ticks in stage X = L, N, A
zNL=0.5; zNN=0.5; zNA=0.5; % [zNX] transmission efficiency from ticks in 
                           % nymph stage to ticks in stage X = L, N, A
zAL=0.5; zAN=0.5; zAA=0.5; % [zAX] transmission efficiency from ticks in 
                           % adult stage to ticks in stage X = L, N, A
pL=0.5; pN=0.5; pA=0.385; % [pX] transmission efficiency from competent 
                          % hosts to ticks in stage X = L, N, A
qL=0.8; qN=0.8; qA=0.8; % [qX] transmission efficiency from ticks in stage 
                        % X = L, N, A to competent hosts
r=0.0005; % transmission efficiency from adults to eggs

% tick-related parameters
e=2000; % mean number of eggs per adult 
sL=0.105; % survival probability from egg to feeding larva
sN=0.105; % survival probability from larva to feeding nymph
sA=0.135; % survival probability from nymph to feeding adult
cLL=30; cLN=20; cLA=0.0505; % [cLX] mean number of ticks in larval stage 
                            % co-feeding with a tick in stage X = L, N, A
cNL=2.55; cNN=1; cNA=0.0505; % [cNX] mean number of ticks in nymph stage 
                             % co-feeding with a tick in stage X = L, N, A
cAL=0.0505; cAN=0.0505; cAA=0.0505; % [cAX] mean number of ticks in adult 
                                    % stage co-feeding with a tick in stage 
                                    % X = L, N, A
g=1; % mean fraction of ticks on a host feeding close enough for 
     % non-systemic transmission
nL=15.5; nN=1.505; nA=0.0501; % [nX] mean number of ticks in stage 
                              % X = L, N, A on competent hosts
dL=3; dN=4.5; dA=10; % [dX] duration of attachment to the host for ticks in 
                     % stage X = L, N, A (days)
hNS=0.5; % fraction of blood meals on hosts competent for non-systemic 
         % transmission
hSYS=0.5; % fraction of blood meals on hosts competent for systemic 
          % transmission

%% Next-generation matrix

k11=sL*sN*sA*e*r;
k12=sN*sA*e*r;
k13=sA*e*r;
k14=e*r;
k15=0;
k21=(sL*zLL*cLL+sL*sN*zNL*cLN+sL*sN*sA*zAL*cLA)*g*hNS;
k22=(sN*zNL*cLN+sN*sA*zAL*cLA)*g*hNS;
k23=sA*zAL*cLA*g*hNS;
k24=0;
k25=pL*v*nL/dL;
k31=(sL*zLN*cNL+sL*sN*zNN*cNN+sL*sN*sA*zAN*cNA)*g*hNS;
k32=(sN*zNN*cNN+sN*sA*zAN*cNA)*g*hNS;
k33=sA*zAN*cNA*g*hNS;
k34=0;
k35=pN*v*nN/dN;
k41=(sL*zLA*cAL+sL*sN*zNA*cAN+sL*sN*sA*zAA*cAA)*g*hNS;
k42=(sN*zNA*cAN+sN*sA*zAA*cAA)*g*hNS;
k43=sA*zAA*cAA*g*hNS;
k44=0;
k45=pA*v*nA/dA;
k51=(sL*qL+sL*sN*qN+sL*sN*sA*qA)*hSYS;
k52=(sN*qN+sN*sA*qA)*hSYS;
k53=sA*qA*hSYS;
k54=0;
k55=0;

K=[k11 k12 k13 k14 k15;
    k21 k22 k23 k24 k25;
    k31 k32 k33 k34 k35;
    k41 k42 k43 k44 k45;
    k51 k52 k53 k54 k55];

display(K)

%% Asymptotic stability and epidemicity analysis

% basic reproduction number
R0=eigs(K,1,'largestabs'); 
display(R0)

% epidemicity index
E0=max(sum(K));
display(E0)

%% Controls
% Assuming simultaneous reductions in all transmission routes

systemic=true; % reduction in systemic transmission
nonsystemic=true; % reduction in non-systemic (co-feeding) tranmission 
transovarial=true; % reduction in trans-ovarial (vertical) transmission

% transmission reductions to avoid endemicity and first-timestep 
% epidemicity
[ctrl_R,ctrl_E]=numerical_controls(K,systemic,nonsystemic,transovarial);

% reproduction number threshold to avoid first-timestep epidemicity
Rstar=eigs(NGM_control(K,systemic,nonsystemic,transovarial,ctrl_E),1,...
    'largestabs'); 
display(Rstar)

% herd immunity threshold
HIT=1-ctrl_R; 
display(HIT)

% herd epidemicity threshold
HET=1-ctrl_E; 
display(HET)

%% Local functions

function Kc=NGM_control(K,systemic,nonsystemic,transovarial,ctrl)
    % NGM_CONTROL - Evaluate the next-generation matrix in the presence
    % of controls

    Kc=K;
    if transovarial
        Kc(1,1:4)=Kc(1,1:4)*ctrl;
    end
    if systemic
        Kc(end,:)=Kc(end,:)*ctrl;
        Kc(:,end)=Kc(:,end)*ctrl;
    end
    if nonsystemic
        Kc(2:4,1:3)=Kc(2:4,1:3)*ctrl;
    end
end

function [ctrl_R,ctrl_E]=numerical_controls(K,systemic,nonsystemic,...
    transovarial)
    % NUMERICAL_CONTROLS - Find the transmission reductions to avoid 
    % endemicity and first-timestep epidemicity

    ctrl_R=fzero(@(ctrl) eigs(NGM_control(K,systemic,nonsystemic,...
        transovarial,ctrl),1,'largestabs')-1,0.5);
    ctrl_E=fzero(@(ctrl) max(sum(NGM_control(K,systemic,nonsystemic,...
        transovarial,ctrl)))-1,0.5);
end
