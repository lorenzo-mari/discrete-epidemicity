clearvars
% close all
% clc

rng default

pathogen=cell(7,1);
pathogen{1}='\textit{B. burgdorferi}'; % Borrelia burgdorferi - Ixodes ricinus
pathogen{2}='TBEV'; % Tick-borne encephalitis virus - Ixodes ricinus
pathogen{3}='\textit{A. phagocytophilum}'; % Anaplasma phagocytophilum - Ixodes ricinus
pathogen{4}='CCHFV'; % Crimean–Congo hemorrhagic fever virus - Hyalomma
pathogen{5}='\textit{R. rickettsii}'; % Rickettsia rickettsii - Dermacentor
pathogen{6}='KFDV'; % Kyasanur Forest disease virus - Haemaphysalis spinigera
pathogen{7}='THOV'; % Thogoto virus - Rhipicephalus appendiculatus

isVirus=[0 1 0 1 0 1 1];
tickGenus={'Ixodes','Ixodes','Ixodes','Hyalomma','Dermacentor','Haemaphysalis','Rhipicephalus'};

ord=[4 6 2 7 3 1 5];
pathogen=pathogen(ord);
isVirus=isVirus(ord);
tickGenus=tickGenus(ord);

%% pathogen-related parameters
pathogen_pars_min=[...
    30      1       6       1       1       2       3; ...
    0.2     0.2     0.01    0.01    0.096   0       0.45; ...
    0.1     0.6     0.2     0.02    0.12    0.4     0.59; ...
    0.1     0.6     0.015   0.05    0.12    0.4     0.59; ...
    0.31    0.6     0.01    0.11    0.12    0.4     0.59; ...
    0.6     0.8     0.041   0.13    0.7     0.13    0.6; ...
    0.6     0.8     0.041   0.13    0.7     0.13    0.6; ...
    0.6     0.8     0.041   0.13    0.7     0.13    0.6; ...
    0       0.01    0       0.14    0.8     0       0];

pathogen_pars_max=[...
    1200    4       2000    7       25      7       4; ...
    0.8     0.8     0.21    0.04    0.144   1       0.81; ...
    0.9     1       0.83    0.67    0.4     0.75    0.82; ...
    0.9     1       0.2     0.67    0.4     0.75    0.82; ...
    0.46    1       0.6     0.33    0.4     0.75    0.82; ...
    1       1       0.2     0.71    1       0.15    0.9; ...
    1       1       0.2     0.71    1       0.15    0.9; ...
    1       1       0.2     0.71    1       0.15    0.9; ...
    0.001   0.03    0.01    0.2     1       0.1     0.01];

np=size(pathogen_pars_min,2);

%% tick-related parameters
tick_pars_min=[...
    1500    4258    4000    250     3000; ...
    0.04    0.02    0.03    0.05    0.03; ...
    0.04    0.02    0.03    0.05    0.03; ...
    0.1     0.07    0.08    0.15    0.08; ...
    24      48      6.24    3.84    107; ...
    0.1     32.8    0       3.84    129; ...
    0.001   0       0       3.76    58; ...
    16      61.6    0       3.84    123; ...
    0       102.4   5.84    3.84    77; ...
    0.001   0       2.9     3.76    55; ...
    0.001   0       0       3.84    103; ...
    0.001   0       2.9     3.84    67; ...
    0.001   30.2    2.9     3.76    54; ...
    1       0.1     0.1     0.1     0.1; ...
    1       60      0.73    0.05    16; ...
    0.01    14      0.25    0.02    6.8; ...
    0.0001  5.5     0.01    0.003   4.56; ...
    2       4       3       2       3; ...
    3       6       4       7       4; ...
    6       6       7       5       7; ...
    0.1     0.1     0.1     0.1     0.1; ...
    0.1     0.1     0.1     0.1     0.1];

tick_pars_max=[...
    2500    9476    6000    1000    5000; ...
    0.17    0.11    0.13    0.25    0.14; ...
    0.17    0.11    0.13    0.25    0.14; ...
    0.17    0.11    0.13    0.25    0.17; ...
    36      72      9.36    5.76    161; ...
    5       49.2    0       5.76    193; ...
    0.1     0       0       5.64    88; ...
    24      92.4    0       5.76    185; ...
    2       153.6   8.76    5.76    115; ...
    0.1     0       4.3     5.64    83; ...
    0.1     0       0       5.76    155; ...
    0.1     0       4.3     5.76    101; ...
    0.1     45.4    4.3     5.64    80; ...
    1       0.5     0.5     0.5     0.5; ...
    30      90      23.7    1.7     24; ...
    3       21      7.7     1.7     10.2; ...
    0.1     8.5     13.9    1.7     6.84; ...
    4       9       7       5       7; ...
    6       11      8       8       8; ...
    14      18      12      7       12; ...
    0.9     0.9     0.9     0.9     0.9; ...
    0.9     0.9     0.9     0.9     0.9];

tick_pars_min=[repmat(tick_pars_min(:,1),1,3) tick_pars_min(:,2:end)];
tick_pars_max=[repmat(tick_pars_max(:,1),1,3) tick_pars_max(:,2:end)];

%% main loop
nn=1e6;

R0=NaN(np,nn);
E0=NaN(np,nn);

for z=1:nn

    % rand pars
    pathogen_pars_rand=pathogen_pars_min+(pathogen_pars_max-pathogen_pars_min).*rand(size(pathogen_pars_min));
    pathogen_pars_rand=pathogen_pars_rand(:,ord);
    tick_pars_rand=tick_pars_min+(tick_pars_max-tick_pars_min).*rand(size(tick_pars_min));
    tick_pars_rand=tick_pars_rand(:,ord);
    
    I=pathogen_pars_rand(1,:);
    theta=pathogen_pars_rand(2,:);
    pL=pathogen_pars_rand(3,:);
    pN=pathogen_pars_rand(4,:);
    pA=pathogen_pars_rand(5,:);
    qL=pathogen_pars_rand(6,:);
    qN=pathogen_pars_rand(7,:);
    qA=pathogen_pars_rand(8,:);
    rA=pathogen_pars_rand(9,:);
    
    thetaLL=theta; 
    thetaLN=theta; 
    thetaLA=theta; 
    thetaNL=theta; 
    thetaNN=theta; 
    thetaNA=theta; 
    thetaAL=theta; 
    thetaAN=theta; 
    thetaAA=theta; 
    
    E=tick_pars_rand(1,:);
    sL=tick_pars_rand(2,:);
    sN=tick_pars_rand(3,:);
    sA=tick_pars_rand(4,:);
    cLL=tick_pars_rand(5,:);
    cNL=tick_pars_rand(6,:);
    cAL=tick_pars_rand(7,:);
    cLN=tick_pars_rand(8,:);
    cNN=tick_pars_rand(9,:);
    cAN=tick_pars_rand(10,:);
    cLA=tick_pars_rand(11,:);
    cNA=tick_pars_rand(12,:);
    cAA=tick_pars_rand(13,:);
    cS=tick_pars_rand(14,:);
    nLH=tick_pars_rand(15,:);
    nNH=tick_pars_rand(16,:);
    nAH=tick_pars_rand(17,:);
    dL=tick_pars_rand(18,:);
    dN=tick_pars_rand(19,:);
    dA=tick_pars_rand(20,:);
    hCN=tick_pars_rand(21,:);
    hCS=tick_pars_rand(22,:);

    for p=1:np

        k11=sL(p)*sN(p)*sA(p)*E(p)*rA(p);
        k12=sN(p)*sA(p)*E(p)*rA(p);
        k13=sA(p)*E(p)*rA(p);
        k14=E(p)*rA(p);
        k15=0;
        k21=(sL(p)*thetaLL(p)*cLL(p)+sL(p)*sN(p)*thetaNL(p)*cLN(p)+sL(p)*sN(p)*sA(p)*thetaAL(p)*cLA(p))*cS(p)*hCN(p);
        k22=(sN(p)*thetaNL(p)*cLN(p)+sN(p)*sA(p)*thetaAL(p)*cLA(p))*cS(p)*hCN(p);
        k23=sA(p)*thetaAL(p)*cLA(p)*cS(p)*hCN(p);
        k24=0;
        k25=pL(p)*I(p)*nLH(p)/dL(p);
        k31=(sL(p)*thetaLN(p)*cNL(p)+sL(p)*sN(p)*thetaNN(p)*cNN(p)+sL(p)*sN(p)*sA(p)*thetaAN(p)*cNA(p))*cS(p)*hCN(p);
        k32=(sN(p)*thetaNN(p)*cNN(p)+sN(p)*sA(p)*thetaAN(p)*cNA(p))*cS(p)*hCN(p);
        k33=sA(p)*thetaAN(p)*cNA(p)*cS(p)*hCN(p);
        k34=0;
        k35=pN(p)*I(p)*nNH(p)/dN(p);
        k41=(sL(p)*thetaLA(p)*cAL(p)+sL(p)*sN(p)*thetaNA(p)*cAN(p)+sL(p)*sN(p)*sA(p)*thetaAA(p)*cAA(p))*cS(p)*hCN(p);
        k42=(sN(p)*thetaNA(p)*cAN(p)+sN(p)*sA(p)*thetaAA(p)*cAA(p))*cS(p)*hCN(p);
        k43=sA(p)*thetaAA(p)*cAA(p)*cS(p)*hCN(p);
        k44=0;
        k45=pA(p)*I(p)*nAH(p)/dA(p);
        k51=(sL(p)*qL(p)+sL(p)*sN(p)*qN(p)+sL(p)*sN(p)*sA(p)*qA(p))*hCS(p);
        k52=(sN(p)*qN(p)+sN(p)*sA(p)*qA(p))*hCS(p);
        k53=sA(p)*qA(p)*hCS(p);
        k54=0;
        k55=0;
    
        K=[k11 k12 k13 k14 k15;
            k21 k22 k23 k24 k25;
            k31 k32 k33 k34 k35;
            k41 k42 k43 k44 k45;
            k51 k52 k53 k54 k55];

        R0(p,z)=eigs(K,1,'largestabs');
        E0(p,z)=max(sum(K));

    end

    if mod(z,1e3)==0
        progress=z/nn*100;
        display(progress)
    end
end

fR=min(1,1./R0);
HIT=1-fR;
fE=min(1,1./E0);
HET=1-fE;
Rstar=R0.*fE;

%%
figure('Renderer','painters','Units','centimeters','Position',[40 10 17.8 14])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

nexttile(1)
boxplot(log10(R0'),'Symbol','')
set(gca,'XLim',[-0.25 np+1.25],'XTick',1:np,'XTickLabel',[],'YLim',[-1 2],'YTick',-1:2,'YTickLabel',{'$10^{-1}$';'$10^0$';'$10^1$';'$10^2$'},'TickDir','out','TickLabelInterpreter','latex')
ylabel('Basic reproduction number, $\mathcal{R}_0$','Interpreter','latex')
title('(a)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(2)
boxplot(log10(E0'),'Symbol','')
set(gca,'XLim',[-0.25 np+1.25],'XTick',1:np,'XTickLabel',[],'YLim',[0 5],'YTick',0:5,'YTickLabel',{'$10^0$';'$10^1$';'$10^2$';'$10^3$';'$10^4$';'$10^5$'},'TickDir','out','TickLabelInterpreter','latex')
ylabel('First-timestep epidemicity, $\mathcal{E}$','Interpreter','latex')
title('(b)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(3)
boxplot(HIT','Symbol','')
set(gca,'XLim',[-0.25 np+1.25],'XTick',1:np,'XTickLabel',pathogen,'YLim',[0 1],'YTick',0:0.2:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel('Herd immunity threshold, HIT','Interpreter','latex')
title('(c)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(4)
boxplot(HET','Symbol','')
set(gca,'XLim',[-0.25 np+1.25],'XTick',1:np,'XTickLabel',pathogen,'YLim',[0.8 1],'YTick',0.8:0.05:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel('Herd epidemicity threshold, HET','Interpreter','latex')
title('(d)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

set(gcf,'NumberTitle','off','Name','Figure S12')
