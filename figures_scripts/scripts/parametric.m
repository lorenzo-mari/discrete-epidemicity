clearvars
% close all
% clc

flag_comp=0;

if flag_comp==1

    n=101;
        
    rec_ave_vect=[5 10 20 30];
    gamma_vect=1./rec_ave_vect;
    ng=numel(gamma_vect);

    %%
    R_1_mat=NaN(n,n,ng);
    R_2_mat=NaN(n,n,ng);
    
    for j=1:ng
        gtd_ave_vect=linspace(0.1,min(30,1/gamma_vect(j)),n);
        gtd_std_vect=linspace(0.1,min(10,1/gamma_vect(j)),n);
        [gtd_ave_mat,gtd_std_mat]=meshgrid(gtd_ave_vect,gtd_std_vect);
    
        R_1_temp=NaN(n,n);
        R_2_temp=NaN(n,n);
    
        for i=1:n^2
            [R_1,R_2]=threshold(gtd_ave_mat(i),gtd_std_mat(i),gamma_vect(j));
            R_1_temp(i)=R_1;
            R_2_temp(i)=R_2;
        end
        R_1_mat(:,:,j)=min(R_1_temp,1,'includenan');
        R_2_mat(:,:,j)=min(R_2_temp,1,'includenan');
    
        progress=round(j/ng*100);
        display(progress)
        
    end

    save parametric
else
    load parametric
end

%%
R0_vect=[0.5 2.5 5 10];
nr=numel(R0_vect);

figure('Renderer','painters','Units','centimeters','Position',[40 10 17.8 16])
htl3=tiledlayout(ng,nr,'TileSpacing','tight','Padding','compact');

z=0;
for j=1:ng
    gtd_ave_vect=linspace(0.1,min(30,1/gamma_vect(j)),n);
    gtd_std_vect=linspace(0.1,min(10,1/gamma_vect(j)),n);
    [gtd_ave_mat,gtd_std_mat]=meshgrid(gtd_ave_vect,gtd_std_vect);
    xt=linspace(0,max(gtd_ave_vect),6);
    xt(1)=min(gtd_ave_vect);
    yt=linspace(0,max(gtd_std_vect),6);
    yt(1)=min(gtd_std_vect);

    for k=1:nr
        z=z+1;

        nexttile(htl3)
        hold on
        [~,hc]=contourf(gtd_ave_mat,gtd_std_mat,max(0,1-R_1_mat(:,:,j)/R0_vect(k),'includenan'),[0 0.5 0.75 0.9 0.95]);
        contour(gtd_ave_mat,gtd_std_mat,max(0,1-R_1_mat(:,:,j)/R0_vect(k),'includenan'),[eps eps],'LineWidth',1,'LineColor','w','LineStyle','-')
        colormap('parula')
        clim([0 1])
        plot([0 gtd_ave_vect],[0 gtd_ave_vect],'w','LineWidth',1)
        set(gca,'XLim',[min(gtd_ave_vect) max(gtd_ave_vect)],'YLim',[min(gtd_std_vect) max(gtd_std_vect)],'Layer','top','TickDir','out','TickLabelInterpreter','latex','XTick',xt,'YTick',yt)
        box on
        ax=gca;
        ax.XAxis.FontSize=8;
        ax.YAxis.FontSize=8;

        if j==1
            title(['$\mathcal{R}_0$ = ',num2str(R0_vect(k))],'Interpreter','latex','FontSize',10)
        end

        if k==4
            yyaxis right
            set(gca,'YTick',[])
            ax=gca;
            ax.YAxis(2).Color='k';
            ylabel(['$1/\gamma$ = ',num2str(rec_ave_vect(j)),' (days)'],'Interpreter','latex','FontSize',10,'Color','k')
        end
    
    end
end

xlabel(htl3,'Average generation time AV (days)','Interpreter','latex','FontSize',10)
ylabel(htl3,'Standard deviation of generation time distribution SD (days)','Interpreter','latex','FontSize',10)

set(gcf,'NumberTitle','off','Name','Figure 4')

%%
function [R_1,R_2]=threshold(gtd_ave,gtd_std,gamma)

    kappa=gtd_ave^2/gtd_std^2;
    
    if kappa<=1
        R_1=NaN;
        R_2=NaN;
    else
        theta=gtd_std^2/gtd_ave;

        thresh=1e-4;
        q1=-1/gamma*log(thresh);
        q2=fzero(@(x) gamcdf(x,kappa,theta)-(1-thresh),100);
        q=ceil(max([q1,q2]));

        f=NaN(1,q);
        for i=1:q
            f(i)=integral(@(tau) gampdf(tau,kappa,theta)./exp(-gamma*tau),i-1,i);
        end

        sigma0=exp(-gamma/2);
        sigma=exp(-gamma);

        R_1=1/sigma0*(1-sigma)/max(f);
        R_2=1/sigma0*sqrt((1-sigma^2)/sum(f.^2));
    end

end
