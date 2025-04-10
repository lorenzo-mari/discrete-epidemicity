clearvars
% close all
% clc

% parameters
gtd_ave=5.2; gtd_std=1.7; lat=4.05;
dur=11.7-lat;
gamma=1/(lat+dur);

% sim pars
wind=90;
tpertave=3.5;
psizeave=1;

Rc_vect=[0.20 0.95];
nR=numel(Rc_vect);
prev=NaN(nR,wind+1);
inc=NaN(nR,wind+1);

% perturbations
rng('default')

delta_t=(exprnd(tpertave*ones(1,wind)));
delta_t(find(cumsum(delta_t)>wind,1,'first'):end)=[];
t_in=round(cumsum(delta_t))+1;

np=numel(t_in);
n_in=exprnd(psizeave*ones(1,np));

rlat=round(lat);
in=zeros(rlat,wind+1);
for i=1:np
        k=ceil(rand*rlat);
        in(k,t_in(i))=in(k,t_in(i))+n_in(i);
end

% Leslie
[f,q]=eval_f(gtd_ave,gtd_std,gamma);
sigma0=exp(-gamma/2);
sigma=exp(-gamma);

for rr=1:nR
    Rc=Rc_vect(rr);

    L=build_leslie(Rc,f,sigma0,sigma);

    % simulation
    input=[in;zeros(q-rlat,wind+1)];

    I=zeros(q,wind+1);
    I(:,1)=input(:,1);
    for i=1:wind
        I(:,i+1)=L*I(:,i)+input(:,i+1);
    end

    prev(rr,:)=sum(I);
    inc(rr,:)=I(1,:)-in(1,:);
end

%%
figure('Renderer','painters','Units','centimeters','Position',[40 5 8.7 8.7])
tiledlayout(4,1,'TileSpacing','tight','Padding','compact')

nexttile(1)
bar(0:wind,in','stacked','LineStyle','n')
set(gca,'XAxisLocation','top','XLim',[0 wind],'XTick',0:30:wind,'XTickLabel','','YLim',[0 ceil(max(n_in)/3)*3+3],'YTick',0:3:ceil(max(n_in)/3)*3+3,'YDir','reverse','TickDir','out','TickLabelInterpreter','latex')
ylabel({'Imports';'(per $10^4$ pop.)'},'Interpreter','latex')
box off
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
leg=legend('1','2','3','4','Location','southeast','Box','off','NumColumns',4,'Interpreter','latex','FontSize',10);
leg.ItemTokenSize=[12 12];
leg.Position=leg.Position-[0 0.01 0 0];
title(leg,'Age of infection (days)')
newcolors=parula(8);
colororder(gca,newcolors(3:6,:))
title('(a)','Interpreter','latex','Units','normalized','Position', [0.05, 0.7, 0],'FontSize',10)

nexttile([3,1])
hold on
for rr=1:nR
    plot(0:wind,prev(rr,:),'LineWidth',1,'DisplayName',['$\mathcal{R}_{\mathrm{c}} =$ ',num2str(Rc_vect(rr))])
end
set(gca,'XLim',[0 wind],'XTick',0:30:wind,'YLim',[0 60],'YTick',0:20:60,'TickDir','out','TickLabelInterpreter','latex') 
xlabel('Time (days)','Interpreter','latex')
ylabel({'Active cases';'(per $10^4$ population)'},'Interpreter','latex')
ax=gca;
ax.XAxis.FontSize=8;
ax.XAxis.Label.FontSize=10;
ax.YAxis.FontSize=8;
ax.YAxis.Label.FontSize=10;
leg=legend('show');
leg.Interpreter='latex';
leg.ItemTokenSize=[9 9];
leg.FontSize=10;
leg.Location='east';
leg.Box='off';
newcolors=lines(2);
colororder(gca,newcolors)
title('(b)','Interpreter','latex','Units','normalized','Position', [0.05, 0.91, 0],'FontSize',10)

set(gcf,'NumberTitle','off','Name','Figure 1')

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
