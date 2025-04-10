clearvars
% close all
% clc

comp_flag=0;

load gtd_data
n=numel(data_names);
gamma=1./duration;

if comp_flag

    E_1=NaN(n,1);
    E_2=NaN(n,1);
    R_thresh_1=NaN(n,1);
    R_thresh_2=NaN(n,1);
    Rtheor_thresh_1=NaN(n,1);
    Rtheor_thresh_2=NaN(n,1);
    HET_1=NaN(n,1);
    HET_2=NaN(n,1);
    HIT=NaN(n,1);
    HETtoHIT_1=NaN(n,1);
    HETtoHIT_2=NaN(n,1);
    maxP_1=NaN(n,1);
    maxP_2=NaN(n,1);
    maxT_1=NaN(n,1);
    maxT_2=NaN(n,1);
    dur_1=NaN(n,1);
    dur_2=NaN(n,1);
    q_vect=NaN(n,1);
    
    R=0.95; 
    
    %%
    for z=1:n
        [f,q]=eval_f(data_ave(z),data_std(z),gamma(z));
        q_vect(z)=q;
        sigma0=exp(-gamma(z)/2);
        sigma=exp(-gamma(z));
        L=build_leslie(R,f,sigma0,sigma);
    
        [R_an_1,R_an_2]=threshold(f,sigma0,sigma);
        R_thresh_1(z)=R_an_1;
        R_thresh_2(z)=R_an_2;
    
        [R_theor_1,R_theor_2]=threshold_theor(data_ave(z),data_std(z),gamma(z),q);
        Rtheor_thresh_1(z)=R_theor_1;
        Rtheor_thresh_2(z)=R_theor_2;
    
        wind_1=find_duration(L,1);
        ae_1=eval_ampenv(L,wind_1,1);
        E_1(z)=ae_1(2);
        [m,j]=max(ae_1(2:end)); 
        maxP_1(z)=m;
        maxT_1(z)=j;
        dur_1(z)=wind_1-1;
    
        wind_2=find_duration(L,2);
        ae_2=eval_ampenv(L,wind_2,2);
        E_2(z)=ae_2(2);
        [m,j]=max(ae_2(2:end)); 
        maxP_2(z)=m;
        maxT_2(z)=j;
        dur_2(z)=wind_2-1;
    
        HET_1(z)=max(0,1-R_thresh_1(z)/R0(z));
        HET_2(z)=max(0,1-R_thresh_2(z)/R0(z));
    
        HIT(z)=max(0,1-1/R0(z));
    
        HETtoHIT_1(z)=HET_1(z)/HIT(z);
        HETtoHIT_2(z)=HET_2(z)/HIT(z);
    
        progress=round(z/n*100);
        display(progress)

        save all_diseases
    end
else
    load all_diseases
end

%% 
figure('Renderer','painters','Units','centimeters','Position',[40 10 17.8 14])
tiledlayout(2,2,'TileSpacing','tight','Padding','compact')

nexttile
bar(1:n,[E_1,maxP_1],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',[],'YLim',[1 2.5],'YTick',1:0.5:2.5,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Maximum amplification';'in one step/overall'},'Interpreter','latex')
leg=legend('$\mathcal{E}$','$\mathcal{A}_{\max}$','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(a)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,[maxT_1,dur_1],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',[],'YLim',[1 1e3],'YTick',[1 10 100 1000],'YScale','log','TickDir','out','TickLabelInterpreter','latex')
ylabel({'Time to peak/','outbreak end (days)'},'Interpreter','latex')
leg=legend('$k_{\max}$','$k_{\mathrm{end}}$','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(b)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,R_thresh_1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',data_names,'XTickLabelRotation',45,'YLim',[0 0.6],'YTick',0:0.2:0.6,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Control reproduction number $\mathcal{R}_*$';'to prevent epidemicity'},'Interpreter','latex')
title('(c)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,[HIT,HET_1],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',data_names,'XTickLabelRotation',45,'YLim',[0 1.2],'YTick',0:0.2:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Herd immunity/';'epidemicity threshold'},'Interpreter','latex')
leg=legend('HIT','HET','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(d)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

set(gcf,'NumberTitle','off','Name','Figure 2')

%%
figure('Renderer','painters','Units','centimeters','Position',[40 10 17.8 14])
tiledlayout(2,2,'TileSpacing','tight','Padding','compact')

nexttile
bar(1:n,[E_2,maxP_2],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',[],'YLim',[1 1.8],'YTick',1:0.2:1.8,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Maximum amplification';'in one step/overall'},'Interpreter','latex')
leg=legend('$\mathcal{E}^{\ell^2}$','$\mathcal{A}_{\max}^{\ell^2}$','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(a)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,[maxT_2,dur_2],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',[],'YLim',[0.8 3e2],'YTick',[1 10 100],'YScale','log','TickDir','out','TickLabelInterpreter','latex')
ylabel({'Time to peak/','outbreak end (days)'},'Interpreter','latex')
leg=legend('$k_{\max}^{\ell^2}$','$k_{\mathrm{end}}^{\ell^2}$','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(b)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,R_thresh_2)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',data_names,'XTickLabelRotation',45,'YLim',[0 1],'YTick',0:0.2:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Control reproduction number $\mathcal{R}_*^{\ell^2}$';'to prevent epidemicity'},'Interpreter','latex')
title('(c)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

nexttile
bar(1:n,[HIT,HET_2],1)
set(gca,'XLim',[-0.2 n+1.2],'XTick',1:n,'XTickLabel',data_names,'XTickLabelRotation',45,'YLim',[0 1.2],'YTick',0:0.2:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel({'Herd immunity/';'epidemicity threshold'},'Interpreter','latex')
leg=legend('HIT','HET$^{\ell^2}$','Location','north','Box','off','NumColumns',2,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
title('(d)','Interpreter','latex','Units','normalized','Position', [0.05, 0.9, 0],'FontSize',10)
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
box off

set(gcf,'NumberTitle','off','Name','Figure S10')

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

function [R_an_1,R_an_2]=threshold(f,sigma0,sigma)

    R_an_1=1/sigma0*(1-sigma)/max(f);
    R_an_2=1/sigma0*sqrt((1-sigma^2)/sum(f.^2));

end

function [R_theor_1,R_theor_2]=threshold_theor(gtd_ave,gtd_std,gamma,q)

    kappa=gtd_ave^2/gtd_std^2;
    theta=gtd_std^2/gtd_ave;

    [~,maxphi]=fminsearch(@(tau) -gampdf(tau,kappa,theta)./exp(-gamma*tau),10);
    maxphi=-1*maxphi;
    intphi2=integral(@(tau) (gampdf(tau,kappa,theta)./exp(-gamma*tau)).^2,0,q);

    R_theor_1=gamma/maxphi;
    R_theor_2=sqrt(2*gamma/intphi2);

end

function ae=eval_ampenv(L,wind,nrm)

    ae=NaN(1,wind+1);
    ae(1)=1;
    for j=1:wind
        ae(j+1)=norm(L^j,nrm);
    end

end

function wind=find_duration(L,nrm) 

    wind=ceil(fzero(@(x) norm(L^ceil(x),nrm)-1,[1 1e5]));

end
