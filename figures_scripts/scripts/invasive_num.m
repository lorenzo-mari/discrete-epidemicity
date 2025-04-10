clearvars
% close all
% clc

F=cell(8,1);
T=cell(8,1);
species=cell(8,1);

% Ailanthus altissima (invasive)
i=1;
F_i=zeros(4); 
F_i(1,4)=3106.28; 
F_i(2,4)=4.16;
F{i}=F_i;
T{i}=[0.75 0 0 0; 
    0.001 0 0 0; 
    0 0.20 0.75 0; 
    0 0 0.01 0.75];
species{i}='\textit{Ailanthus altissima}';

% Alliaria petiolata (invasive)
i=2;
F_i=zeros(3);
F_i(1,3)=658.59;
F_i(2,3)=105.58;
F{i}=F_i;
T{i}=[0.65 0 0; 
    0.04 0 0; 
    0 0.39 0];
species{i}='\textit{Alliaria petiolata}';

% Allium vineale (invasive)
i=3;
F_i=zeros(3);
F_i(1,3)=1.96;
F_i(2,3)=8.38;
F{i}=F_i;
T{i}=[0.24 0 0; 0.001 0.01 0.02; 0 0.14 0.14];
species{i}='\textit{Allium vineale}';

% Cerastium fontanum (invasive)
i=4;
F_i=zeros(2);
F_i(1,2)=87.64; 
F_i(2,2)=4.44;
F{i}=F_i;
T{i}=[0.30 0; 0.001 0];
species{i}='\textit{Cerastium fontanum}';

% Lonicera maackii (invasive)
i=5;
F_i=zeros(4);
F_i(1,4)=141.72;
F_i(2,4)=16.10;
F{i}=F_i;
T{i}=[0.44 0 0 0; 0.004 0.21 0 0; 0 0.44 0.83 0.18; 0 0 0.15 0.82];
species{i}='\textit{Lonicera maackii}';

% Rosa multiflora (invasive)
i=6;
F_i=zeros(4);
F_i(1,4)=24.46;
F_i(2,4)=6.46;
F{i}=F_i;
T{i}=[0.12 0 0 0; 0.001 0.14 0 0; 0 0.86 0.57 0.28; 0 0 0.33 0.72];
species{i}='\textit{Rosa multiflora}';

% Taraxacum officinale (invasive)
i=7;
F_i=zeros(3);
F_i(1,3)=3039.24;
F_i(2,3)=706.80;
F{i}=F_i;
T{i}=[0.18 0 0; 3.33e-3 0.39 0; 0 0.13 0.55];
species{i}='\textit{Taraxacum officinale}';

% Veronica arvensis (invasive)
i=8;
F_i=zeros(2);
F_i(1,2)=89.69;
F{i}=F_i;
T{i}=[0.38 0; 0.001 0.59];
species{i}='\textit{Veronica arvensis}';

n=i;
lambda=NaN(1,n);
R0=NaN(1,n);
E0=NaN(1,n);
Amax=NaN(1,n);
alpha=NaN(1,n);

for i=1:n

    A=F{i}+T{i};
    Q=F{i}/(diag(ones(size(T{i},1),1))-T{i});

    lambda(i)=eigs(A,1,'largestabs');
    R0(i)=eigs(Q,1,'largestabs');
    E0(i)=max(sum(A));
    Amax(i)=max(ampenv(A));
    alpha(i)=findalpha(F{i},T{i});

end
Amax(R0>1)=NaN;

%%
figure('Renderer','painters','Units','centimeters','Position',[40 10 17.8 14])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

nexttile(1)
bar(1:n,R0)
set(gca,'XLim',[-0.25 n+1.25],'XTick',1:n,'XTickLabel',[],'YScale','log','YLim',[0.1 1e3],'YTick',10.^(-1:3),'TickDir','out','TickLabelInterpreter','latex')
ylabel('Net reproduction value, $\mathcal{N}$','Interpreter','latex')
title('(a)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(2)
bar(1:n,E0)
set(gca,'XLim',[-0.25 n+1.25],'XTick',1:n,'XTickLabel',[],'YScale','log','YLim',[1 1e4],'YTick',10.^(0:4),'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Maximum first-timestep';'reactivity, $\mathcal{G}$'},'Interpreter','latex')
title('(b)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(3)
bar(1:n,alpha.*R0)
set(gca,'XLim',[-0.25 n+1.25],'XTick',1:n,'XTickLabel',species,'YLim',[0 0.12],'YTick',0:0.03:0.12,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Net reproduction value threshold';'to prevent transient growth, $\mathcal{N}_*$'},'Interpreter','latex')
title('(c)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile(4)
bar(1:n,1-alpha)
set(gca,'XLim',[-0.25 n+1.25],'XTick',1:n,'XTickLabel',species,'YLim',[0.9 1],'YTick',0.9:0.02:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Relative reduction in reproduction value';'to prevent transient growth, RRR'},'Interpreter','latex')
title('(d)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

set(gcf,'NumberTitle','off','Name','Figure S13')

%%
function ae=ampenv(M)

    n=100;
    t=0:n;
    ae=NaN(1,n+1);
    ae(1)=1;

    for i=2:n+1
        ae(i)=max(sum(M^t(i)));
    end

end

function alpha=findalpha(F,T)

    Ec=@(alpha) max(sum(alpha*F+T))-1;
    alpha=fzero(Ec,[0 1]);

end
