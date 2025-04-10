clearvars
% close all
% clc

% parameters
gtd_ave=5.2;
gtd_std=1.7;
gamma=1/11.7; 
R0=3.6;

% HIT
HIT=1-1/R0;

% Rc example
Rc=0.95;

% eval Leslie matrix
[f,q]=eval_f(gtd_ave,gtd_std,gamma);
sigma0=exp(-gamma/2);
sigma=exp(-gamma);
L=build_leslie(Rc,f,sigma0,sigma);

[E_2_an,R_2_an,E_1_an,R_1_an]=analytics(Rc,f,sigma0,sigma);

ae_1=eval_ampenv(L,300,1);
ae_2=eval_ampenv(L,300,2);
dur_1=find(ae_1<1,1,'first')-1;
dur_2=find(ae_2<1,1,'first')-1;

wind=ceil(dur_1/50)*50;
ae_1(wind+2:end)=[];
ae_2(wind+2:end)=[];
E_1_ae=ae_1(2);
E_2_ae=ae_2(2);

% optimal perturbation at time 0 (l2 norm)
[~,~,w]=svd(L); 
I0_Z_2=abs(w(:,1));

% optimal perturbation overall (l2 norm)
[ae_2_max,T2]=max(ae_2(2:end)); 
[~,~,w]=svd(L^T2);
I0_T_2=abs(w(:,1));

% optimal perturbation at time 0 (l1 norm)
I0_Z_1=zeros(q,1);
[~,i]=max(sum(L));
I0_Z_1(i)=1;

% optimal perturbation overall (l1 norm)
[ae_1_max,T1]=max(ae_1(2:end)); 
I0_T_1=zeros(q,1);
[~,i]=max(sum(L^T1));
I0_T_1(i)=1;

% HET
HET_1=1-R_1_an/R0;
HET2HIT_1=HET_1/HIT;

HET_2=1-R_2_an/R0;
HET2HIT_2=HET_2/HIT;

% simulations
I_Z_2=zeros(q,wind+1);
I_Z_2(:,1)=I0_Z_2;
nI_Z=NaN(1,wind+1);
nI_Z(1)=norm(I0_Z_2);

I_T_2=zeros(q,wind+1);
I_T_2(:,1)=I0_T_2;
nI_T=NaN(1,wind+1);
nI_T(1)=norm(I0_T_2);

I_Z_1=zeros(q,wind+1);
I_Z_1(:,1)=I0_Z_1;
pI_Z=NaN(1,wind+1);
pI_Z(1)=sum(I0_Z_1);

I_T_1=zeros(q,wind+1);
I_T_1(:,1)=I0_T_1;
pI_T=NaN(1,wind+1);
pI_T(1)=sum(I0_T_1);

for i=1:wind
    I_Z_2(:,i+1)=L*I_Z_2(:,i);
    nI_Z(i+1)=norm(I_Z_2(:,i+1));

    I_T_2(:,i+1)=L*I_T_2(:,i);
    nI_T(i+1)=norm(I_T_2(:,i+1));

    I_Z_1(:,i+1)=L*I_Z_1(:,i);
    pI_Z(i+1)=sum(I_Z_1(:,i+1));

    I_T_1(:,i+1)=L*I_T_1(:,i);
    pI_T(i+1)=sum(I_T_1(:,i+1));
end

% R threshold
n=101;
R_vect=linspace(0,0.999,n);
E_1_vect=NaN(1,n);
A_1_vect=NaN(1,n);
T_1_vect=NaN(1,n);
W_1_vect=NaN(1,n);
E_2_vect=NaN(1,n);
A_2_vect=NaN(1,n);
T_2_vect=NaN(1,n);
W_2_vect=NaN(1,n);

for k=1:n
    R_k=R_vect(k);

    L_k=L;
    L_k(1,:)=L(1,:)/Rc*R_k;

    ae_1_k=eval_ampenv(L_k,1,1);
    E_1_vect(k)=ae_1_k(2);

    if ae_1_k(2)>1
        wind_1=find_duration(L_k,1);
        ae_1_k=eval_ampenv(L_k,wind_1,1);      
        [a,t]=max(ae_1_k(2:end));
        A_1_vect(k)=a;
        T_1_vect(k)=t;
        W_1_vect(k)=wind_1-1;
    else
        A_1_vect(k)=ae_1_k(2);
        T_1_vect(k)=0;
        W_1_vect(k)=0;
    end

    ae_2_k=eval_ampenv(L_k,1,2);
    E_2_vect(k)=ae_2_k(2);

    if ae_2_k(2)>1
        wind_2=find_duration(L_k,2);
        ae_2_k=eval_ampenv(L_k,wind_2,2); 
        [a,t]=max(ae_2_k(2:end));
        A_2_vect(k)=a;
        T_2_vect(k)=t;
        W_2_vect(k)=wind_2-1;
    else
        A_2_vect(k)=ae_2_k(2);
        T_2_vect(k)=0;
        W_2_vect(k)=0;
    end        
    
end

%%
cols=lines(2);

figure('Renderer','painters','Units','centimeters','Position',[40 5 16 16])
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')

nexttile
hold on
plot(0:wind,pI_Z/pI_Z(1),'LineWidth',1,'Color',cols(1,:))
plot(0:wind,pI_T/pI_T(1),'LineWidth',1,'Color',cols(1,:),'LineStyle','--')
plot(0:wind,ae_1,'LineWidth',3,'Color',[cols(2,:) 0.5])
set(gca,'XLim',[0 wind],'XTick',0:30:wind,'YLim',[0.75 2.25],'YTick',0.5:0.5:2.5,'TickDir','out','TickLabelInterpreter','latex')
xlabel('Time $k$ (days)','Interpreter','latex')
ylabel({'Normalized prevalence $\mathcal{P}(k)$'},'Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10; 
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
leg=legend(...
    'Most crit. pert. at $k = 0$',...
    'Most crit. pert. overall',...
    'Amplification envelope $\mathcal{A}(k)$',...
    'Location','northeast','Box','off','Interpreter','latex','FontSize',10,'AutoUpdate','off');
leg.ItemTokenSize=[20 20];
title('(a)','Interpreter','latex','Units','normalized','Position', [0.04, 0.86, 0],'FontSize',10)

scatter(1,E_1_an,'o','LineWidth',1,'MarkerEdgeColor',cols(2,:),'MarkerFaceColor',cols(2,:))
text(1+1,E_1_an-0.05,'$(1,\mathcal{E})$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
scatter(T1,ae_1_max,'o','LineWidth',1,'MarkerEdgeColor',cols(2,:),'MarkerFaceColor','w')
text(T1+1,ae_1_max,'$(k_{\max},\mathcal{A}_{\max})$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom')
scatter(dur_1,ae_1(dur_1+1),'o','LineWidth',1,'MarkerEdgeColor',cols(2,:),'MarkerFaceColor','w')
text(dur_1+2,ae_1(dur_1+1),'$(k_{\mathrm{end}},1)$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom')

nexttile
hold on
plot(R_vect,E_1_vect,'LineWidth',1,'Color',cols(1,:))
plot(R_vect,A_1_vect,'LineWidth',1,'Color',cols(2,:))
set(gca,'XLim',[0 1],'XTick',0:0.2:1,'YLim',[0.5 2.5],'YTick',0.5:0.5:2.5,'TickDir','out','TickLabelInterpreter','latex')
xlabel('Control reproduction number $\mathcal{R}_{\mathrm{c}}$','Interpreter','latex')
ylabel('Amplification factor','Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10; 
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
leg=legend(...
    'First-timestep epidemicity $\mathcal{E}$',...
    'Maximum epidemicity $\mathcal{A}_{\max}$',...
    'Location','north','Box','off','Interpreter','latex','FontSize',10,'AutoUpdate','off');
leg.ItemTokenSize=[20 20];
title('(b)','Interpreter','latex','Units','normalized','Position', [0.04, 0.86, 0],'FontSize',10)

scatter(R_1_an,1,'LineWidth',1,'MarkerEdgeColor',cols(1,:),'MarkerFaceColor',cols(1,:))
text(R_1_an,1-0.08,'$(\mathcal{R}_*,1)$','Interpreter','latex','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','top')

nexttile
hold on
plot(R_vect,T_1_vect,'LineWidth',1,'Color',cols(1,:))
plot(R_vect,W_1_vect,'LineWidth',1,'Color',cols(2,:))
set(gca,'XLim',[0 1],'XTick',0:0.2:1,'YLim',[1 1e3],'YTick',10.^(0:3),'YScale','log','TickDir','out','TickDir','out','TickLabelInterpreter','latex')
xlabel('Control reproduction number $\mathcal{R}_{\mathrm{c}}$','Interpreter','latex')
ylabel('Time to\ldots (days)','Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10; 
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
leg=legend(...
    '\ldots~largest outbreak''s peak $k_{\max}$',...
    '\ldots~longest outbreak''s end $\mathcal{A}_{\max}$',...
    'Location','west','Box','off','Interpreter','latex','FontSize',10,'AutoUpdate','off');
leg.ItemTokenSize=[20 20];
title('(c)','Interpreter','latex','Units','normalized','Position', [0.04, 0.86, 0],'FontSize',10)

set(gcf,'NumberTitle','off','Name','Figure S1')

%%
function [f,q]=eval_f(gtd_ave,gtd_std,gamma)

    kappa=gtd_ave^2/gtd_std^2;
    theta=gtd_std^2/gtd_ave;

    thresh=1e-4;
    q1=-1/gamma*log(thresh);
    q2=fzero(@(x) gamcdf(x,kappa,theta)-(1-thresh),100);
    q=ceil(max([q1,q2]));

    f=NaN(1,q);
    for i=1:q
        f(i)=integral(@(tau) gampdf(tau,kappa,theta)./exp(-gamma*tau),i-1,i);
    end

end

function L=build_leslie(R,f,sigma0,sigma)

    q=length(f);    
    
    L=zeros(q,q);
    L(1,:)=R*sigma0*f;
    L(2:q,1:q-1)=eye(q-1)*sigma;

end

function ae=eval_ampenv(L,wind,nrm)

    ae=NaN(1,wind+1);
    ae(1)=1;

    for j=1:wind
        ae(j+1)=norm(L^j,nrm);
    end

end

function [E_2_an,R_2_an,E_1_an,R_1_an]=analytics(R,f,sigma0,sigma)

    E_2_an=sqrt(1/2*(sigma^2+R^2*sigma0^2*sum(f.^2)+sqrt((sigma^2+R^2*sigma0^2*sum(f.^2))^2-4*sigma^2*R^2*sigma0^2*f(end)^2)));
    R_2_an=1/sigma0*sqrt((1-sigma^2)/sum(f.^2));

    E_1_an=sigma0*R*max(f)+sigma;
    R_1_an=1/sigma0*(1-sigma)/max(f);

end

function wind=find_duration(L,nrm)

    wind=ceil(fzero(@(x) norm(L^ceil(x),nrm)-1,[1 1e8]));

end
