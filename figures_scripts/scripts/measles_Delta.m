clearvars
% close all
% clc

% parameters
gtd_ave=11.7;
gtd_std=2.5;
gamma=1/14;
R0=14.5;

Rc=0.95;

Delta_ex=[1/4 1 3.5]; % days
nd=numel(Delta_ex);

figure('Renderer','painters','Units','centimeters','Position',[40 5 12 15])
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')

cols=lines(nd);
hl=NaN(2,nd);

%%
for j=1:nd
    [f,q]=eval_f_Delta(gtd_ave,gtd_std,gamma,Delta_ex(j));
    sigma0=exp(-gamma*Delta_ex(j)/2);
    sigma=exp(-gamma*Delta_ex(j));
    L=build_leslie(Rc,f,sigma0,sigma);
    
    [E_2_an,~,E_1_an,~]=analytics(Rc,f,sigma0,sigma);
    
    wind=30;
    ae_1=eval_ampenv(L,ceil(wind/Delta_ex(j)),1);
    ae_2=eval_ampenv(L,ceil(wind/Delta_ex(j)),2);

    nexttile(1)
    hold on
    hl(1,j)=plot(linspace(0,wind,ceil(wind/Delta_ex(j))+1),ae_1,'LineWidth',1,'Color',cols(j,:));
    scatter(Delta_ex(j),E_1_an,'o','MarkerEdgeColor',cols(j,:),'MarkerFaceColor',cols(j,:))
    set(gca,'XLim',[0 wind],'XTick',0:10:wind,'YLim',[0.5 2],'YTick',0.5:0.5:2,'TickDir','out','TickLabelInterpreter','latex')
    xlabel('Time $t$ (days)','Interpreter','latex')
    ylabel('$\mathcal{A}(t)$','Interpreter','latex')
    ax=gca;
    ax.XAxis.FontSize=8;
    ax.XAxis.Label.FontSize=10; 
    ax.YAxis.FontSize=8;
    ax.YAxis.Label.FontSize=10;
    if j==2
        text(Delta_ex(j)+0.5,E_1_an,'$(\Delta,\mathcal{E}(\Delta))$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
    end
    if j==3
        leg=legend(hl(1,:),['$\Delta$ = ',num2str(Delta_ex(1)),' days'],...
            ['$\Delta$ = ',num2str(Delta_ex(2)),' day'],...
            ['$\Delta$ = ',num2str(Delta_ex(3)),' days'],...
            'Box','off','Location','south','Interpreter','latex','FontSize',10,'NumColumns',3);
        leg.ItemTokenSize=[20 15];
        title('(a)','Interpreter','latex','Units','normalized','Position',[0.04, 0.86, 0],'FontSize',10)
    end

    progress=j/nd*100;
    display(progress)
end

%% discretization step
nd=41;
Delta_vect=logspace(-3,1,nd); % days
E_1_vect=NaN(1,nd);
E_2_vect=NaN(1,nd);
R_1_vect=NaN(1,nd);
R_2_vect=NaN(1,nd);
D_1_vect=NaN(1,nd);
D_2_vect=NaN(1,nd);

[maxphi,intphi2]=phi_funs(gtd_ave,gtd_std,gamma,q);
R_1_lim=gamma/maxphi;
R_2_lim=sqrt(2*gamma/intphi2);
D_1_lim=maxphi*Rc-gamma;
D_2_lim=Rc^2/2*intphi2-gamma;
HET_1_lim=1-R_1_lim/R0;
HET_2_lim=1-R_2_lim/R0;

for j=1:nd
    [f,q]=eval_f_Delta(gtd_ave,gtd_std,gamma,Delta_vect(j));
    sigma0=exp(-gamma*Delta_vect(j)/2);
    sigma=exp(-gamma*Delta_vect(j));
    
    [E_2_an,R_2_an,E_1_an,R_1_an]=analytics(Rc,f,sigma0,sigma);

    E_1_vect(j)=E_1_an;
    E_2_vect(j)=E_2_an;

    R_1_vect(j)=R_1_an;
    R_2_vect(j)=R_2_an;

    D_1_vect(j)=(E_1_vect(j)-1)/Delta_vect(j);
    D_2_vect(j)=(E_2_vect(j)-1)/Delta_vect(j);

    progress=j/nd*100;
    display(progress)

end

HET_1_vect=1-R_1_vect/R0;
HET_2_vect=1-R_2_vect/R0;

%%
nexttile(2)
hold on
box on
yyaxis left
plot(Delta_vect,E_1_vect,'LineWidth',1)
scatter(Delta_vect(1),1,'o','MarkerEdgeColor',cols(1,:),'MarkerFaceColor',cols(1,:))
text(Delta_vect(1)+0.0002,1-0.05,'$(\delta,\mathcal{E}(\delta))$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
set(gca,'XScale','log','YLim',[0.5 2.5],'YTick',0.5:0.5:2.5,'TickDir','out','TickLabelInterpreter','latex')
xlabel('Discretization timestep $\Delta$ (days)','Interpreter','latex')
ylabel('$\mathcal{E}$','Interpreter','latex')
yyaxis right
plot(Delta_vect,D_1_vect,'LineWidth',1)
scatter(Delta_vect(1),D_1_lim,'o','MarkerEdgeColor',cols(2,:),'MarkerFaceColor',cols(2,:))
text(Delta_vect(1)+0.0002,D_1_lim-0.01,'$(\delta,(\mathcal{E}(\delta)-1)/\delta)$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
set(gca,'YLim',[0 0.4],'YTick',0:0.1:0.4,'TickDir','out','TickLabelInterpreter','latex')
ylabel('$(\mathcal{E} - 1) / \Delta$ (days)$^{-1}$','Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10; 
ax.YAxis(1).FontSize=8;
ax.YAxis(2).FontSize=8;
ax.YAxis(1).Label.FontSize=10;
ax.YAxis(2).Label.FontSize=10;
title('(b)','Interpreter','latex','Units','normalized','Position',[0.04, 0.86, 0],'FontSize',10)

nexttile(3)
hold on
box on
yyaxis left
plot(Delta_vect,R_1_vect,'LineWidth',1)
scatter(Delta_vect(1),R_1_lim,'o','MarkerEdgeColor',cols(1,:),'MarkerFaceColor',cols(1,:))
text(Delta_vect(1)+0.0002,R_1_lim-0.015,'$(\delta,\mathcal{R}_*(\delta))$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
set(gca,'XScale','log','YLim',[0.1 0.5],'YTick',0.1:0.1:0.5,'TickDir','out','TickLabelInterpreter','latex')
xlabel('Discretization timestep $\Delta$ (days)','Interpreter','latex')
ylabel('$\mathcal{R}_*$','Interpreter','latex')
yyaxis right
plot(Delta_vect,HET_1_vect,'LineWidth',1)
scatter(Delta_vect(1),HET_1_lim,'o','MarkerEdgeColor',cols(2,:),'MarkerFaceColor',cols(2,:))
text(Delta_vect(1)+0.0002,HET_1_lim-0.001,'$(\delta,\mathrm{HET}(\delta))$','Interpreter','latex','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top')
set(gca,'YLim',[0.96 1],'YTick',0.96:0.01:1,'TickDir','out','TickLabelInterpreter','latex')
ylabel('HET','Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10; 
ax.YAxis(1).FontSize=8;
ax.YAxis(2).FontSize=8;
ax.YAxis(1).Label.FontSize=10;
ax.YAxis(2).Label.FontSize=10;
title('(c)','Interpreter','latex','Units','normalized','Position',[0.04, 0.86, 0],'FontSize',10)

set(gcf,'NumberTitle','off','Name','Figure S3')

%%
function [f,q]=eval_f_Delta(gtd_ave,gtd_std,gamma,Delta)

    kappa=gtd_ave^2/gtd_std^2;
    theta=gtd_std^2/gtd_ave;

    thresh=1e-4;
    q1=-1/gamma*log(thresh);
    q2=fzero(@(x) gamcdf(x,kappa,theta)-(1-thresh),100);
    q=ceil(ceil(max([q1,q2]))/Delta);

    f=NaN(1,q);
    for i=1:q
        f(i)=integral(@(tau) gampdf(tau,kappa,theta)./exp(-gamma*tau),Delta*(i-1),Delta*i);
    end
    
end

function [maxphi,intphi2]=phi_funs(gtd_ave,gtd_std,gamma_par,q)

    kappa=gtd_ave^2/gtd_std^2;
    theta=gtd_std^2/gtd_ave;

    [~,maxphi]=fminsearch(@(tau) -gampdf(tau,kappa,theta)./exp(-gamma_par*tau),10);
    maxphi=-1*maxphi;
    intphi2=integral(@(tau) (gampdf(tau,kappa,theta)./exp(-gamma_par*tau)).^2,0,q);

    % maxphi_an=((kappa-1)/(1/theta-gamma_par))^(kappa-1)*exp(1-kappa)/theta^kappa/gamma(kappa)

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
