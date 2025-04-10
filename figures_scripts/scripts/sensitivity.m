clearvars
% close all
% clc

flag_comp=0;

if flag_comp==1
    load gtd_data
    gamma=1./duration;
    n=numel(gamma);
    
    Rc=0.95;
    
    ns=11;
    sens_vect=1+linspace(-1,1,ns)*0.2;
    
    ave_mat=NaN(n,ns);
    E1_ave=NaN(n,ns);
    E2_ave=NaN(n,ns);
    R1_ave=NaN(n,ns);
    R2_ave=NaN(n,ns);
    HET1_ave=NaN(n,ns);
    HET2_ave=NaN(n,ns);
    maxP1_ave=NaN(n,ns);
    maxP2_ave=NaN(n,ns);
    maxT1_ave=NaN(n,ns);
    maxT2_ave=NaN(n,ns);
    dur1_ave=NaN(n,ns);
    dur2_ave=NaN(n,ns);
    
    std_mat=NaN(n,ns);
    E1_std=NaN(n,ns);
    E2_std=NaN(n,ns);
    R1_std=NaN(n,ns);
    R2_std=NaN(n,ns);
    HET1_std=NaN(n,ns);
    HET2_std=NaN(n,ns);
    maxP1_std=NaN(n,ns);
    maxP2_std=NaN(n,ns);
    maxT1_std=NaN(n,ns);
    maxT2_std=NaN(n,ns);
    dur1_std=NaN(n,ns);
    dur2_std=NaN(n,ns);
    
    infper_mat=NaN(n,ns);
    E1_infper=NaN(n,ns);
    E2_infper=NaN(n,ns);
    R1_infper=NaN(n,ns);
    R2_infper=NaN(n,ns);
    HET1_infper=NaN(n,ns);
    HET2_infper=NaN(n,ns);
    maxP1_infper=NaN(n,ns);
    maxP2_infper=NaN(n,ns);
    maxT1_infper=NaN(n,ns);
    maxT2_infper=NaN(n,ns);
    dur1_infper=NaN(n,ns);
    dur2_infper=NaN(n,ns);
    
    R0_mat=NaN(n,ns);
    E1_R0=NaN(n,ns);
    E2_R0=NaN(n,ns);
    R1_R0=NaN(n,ns);
    R2_R0=NaN(n,ns);
    HET1_R0=NaN(n,ns);
    HET2_R0=NaN(n,ns);
    maxP1_R0=NaN(n,ns);
    maxP2_R0=NaN(n,ns);
    maxT1_R0=NaN(n,ns);
    maxT2_R0=NaN(n,ns);
    dur1_R0=NaN(n,ns);
    dur2_R0=NaN(n,ns);
    
    k=0;
    for z=1:n
    
        ave_mat(z,:)=data_ave(z)*sens_vect;
        std_mat(z,:)=data_std(z)*sens_vect;
        infper_mat(z,:)=duration(z)*sens_vect;
        R0_mat(z,:)=R0(z)*sens_vect;
    
        for w=1:ns
            % gtd_ave
            f=eval_f(ave_mat(z,w),data_std(z),gamma(z));
            sigma0=exp(-gamma(z)/2);
            sigma=exp(-gamma(z));
            L=build_leslie(Rc,f,sigma0,sigma);
    
            wind1=find_duration(L,1);
            dur1_ave(z,w)=wind1-1;
            wind2=find_duration(L,2);
            dur2_ave(z,w)=wind2-1;
    
            [R1_an,R2_an]=threshold(f,sigma0,sigma);
            R1_ave(z,w)=R1_an;
            HET1_ave(z,w)=max([0 1-R1_an/R0(z)]);
            R2_ave(z,w)=R2_an;
            HET2_ave(z,w)=max([0 1-R2_an/R0(z)]);
    
            ae_1=eval_ampenv(L,wind1,1);
            [m,j]=max(ae_1(2:end)); 
            E1_ave(z,w)=ae_1(2);
            maxP1_ave(z,w)=m;
            maxT1_ave(z,w)=j;
    
            ae_2=eval_ampenv(L,wind2,2);
            [m,j]=max(ae_2(2:end)); 
            E2_ave(z,w)=ae_2(2);
            maxP2_ave(z,w)=m;
            maxT2_ave(z,w)=j;
    
            % gtd_std
            f=eval_f(data_ave(z),std_mat(z,w),gamma(z));
            sigma0=exp(-gamma(z)/2);
            sigma=exp(-gamma(z));
            L=build_leslie(Rc,f,sigma0,sigma);
    
            wind1=find_duration(L,1);
            dur1_std(z,w)=wind1-1;
            wind2=find_duration(L,2);
            dur2_std(z,w)=wind2-1;
    
            [R1_an,R2_an]=threshold(f,sigma0,sigma);
            R1_std(z,w)=R1_an;
            HET1_std(z,w)=max([0 1-R1_an/R0(z)]);
            R2_std(z,w)=R2_an;
            HET2_std(z,w)=max([0 1-R2_an/R0(z)]);
    
            ae_1=eval_ampenv(L,wind1,1);
            [m,j]=max(ae_1(2:end)); 
            E1_std(z,w)=ae_1(2);
            maxP1_std(z,w)=m;
            maxT1_std(z,w)=j;
    
            ae_2=eval_ampenv(L,wind2,2);
            [m,j]=max(ae_2(2:end)); 
            E2_std(z,w)=ae_2(2);
            maxP2_std(z,w)=m;
            maxT2_std(z,w)=j;
    
            % duration
            f=eval_f(data_ave(z),data_std(z),1/infper_mat(z,w));
            sigma0=exp(-1/infper_mat(z,w)/2);
            sigma=exp(-1/infper_mat(z,w));
            L=build_leslie(Rc,f,sigma0,sigma);
    
            wind1=find_duration(L,1);
            dur1_infper(z,w)=wind1-1;
            wind2=find_duration(L,2);
            dur2_infper(z,w)=wind2-1;
    
            [R1_an,R2_an]=threshold(f,sigma0,sigma);
            R1_infper(z,w)=R1_an;
            HET1_infper(z,w)=max([0 1-R1_an/R0(z)]);
            R2_infper(z,w)=R2_an;
            HET2_infper(z,w)=max([0 1-R2_an/R0(z)]);
    
            ae_1=eval_ampenv(L,wind1,1);
            [m,j]=max(ae_1(2:end)); 
            E1_infper(z,w)=ae_1(2);
            maxP1_infper(z,w)=m;
            maxT1_infper(z,w)=j;
    
            ae_2=eval_ampenv(L,wind2,2);
            [m,j]=max(ae_2(2:end)); 
            E2_infper(z,w)=ae_2(2);
            maxP2_infper(z,w)=m;
            maxT2_infper(z,w)=j;
    
            % R0
            f=eval_f(data_ave(z),data_std(z),gamma(z));
            sigma0=exp(-gamma(z)/2);
            sigma=exp(-gamma(z));
            L=build_leslie(Rc,f,sigma0,sigma);
    
            wind1=find_duration(L,1);
            dur1_R0(z,w)=wind1-1;
            wind2=find_duration(L,2);
            dur2_R0(z,w)=wind2-1;
    
            [R1_an,R2_an]=threshold(f,sigma0,sigma);
            R1_R0(z,w)=R1_an;
            HET1_R0(z,w)=max([0 1-R1_an/R0_mat(z,w)]);
            R2_R0(z,w)=R2_an;
            HET2_R0(z,w)=max([0 1-R2_an/R0_mat(z,w)]);
    
            ae_1=eval_ampenv(L,wind1,1);
            [m,j]=max(ae_1(2:end)); 
            E1_R0(z,w)=ae_1(2);
            maxP1_R0(z,w)=m;
            maxT1_R0(z,w)=j;
    
            ae_2=eval_ampenv(L,wind2,2);
            [m,j]=max(ae_2(2:end)); 
            E2_R0(z,w)=ae_2(2);
            maxP2_R0(z,w)=m;
            maxT2_R0(z,w)=j;
    
            % progress
            k=k+1;
            progress=k/n/ns*100;
            display(progress)
        end
    
    end

    save sensitivity
else
    load sensitivity
end

ref=ceil(ns/2);

%%
sens_value=20;

i1=find(sens_vect==1-sens_value/100);
i2=find(sens_vect==1+sens_value/100);

for j=1:6
    switch j
        case 1
            sens_ave=E1_ave;
            sens_std=E1_std;
            sens_infper=E1_infper;
            sens_R0=E1_R0;
        case 2
            sens_ave=maxP1_ave;
            sens_std=maxP1_std;
            sens_infper=maxP1_infper;
            sens_R0=maxP1_R0;
        case 3
            sens_ave=maxT1_ave;
            sens_std=maxT1_std;
            sens_infper=maxT1_infper;
            sens_r0=maxT1_R0;
        case 4
            sens_ave=dur1_ave;
            sens_std=dur1_std;
            sens_infper=dur1_infper;
            sens_R0=dur1_R0;
        case 5
            sens_ave=R1_ave;
            sens_std=R1_std;
            sens_infper=R1_infper;
            sens_R0=R1_R0;
        case 6
            sens_ave=HET1_ave;
            sens_std=HET1_std;
            sens_infper=HET1_infper;
            sens_R0=HET1_R0;
    end
    
    if j<6
        figure('Renderer','painters','Units','centimeters','Position',[40 5 8 18])
        tiledlayout(3,1,'TileSpacing','tight','Padding','compact')
    else
        figure('Renderer','painters','Units','centimeters','Position',[40 5 16 12])
        tiledlayout(2,2,'TileSpacing','tight','Padding','compact')
    end

    tornado=(sens_ave(:,[i1 i2])-sens_ave(:,ref))./sens_ave(:,ref)*100;
    [~,s]=sort(sum(abs(tornado),2),'ascend');
    xl=max(abs(tornado(:)));
    if xl<100
        xl=max([10 ceil(xl/10)*10]);
    else
        xl=ceil(xl/100)*100;
    end

    nexttile
    hold on
    barh(1:n,tornado(s,1))
    barh(1:n,tornado(s,2))
    set(gca,'YTick',1:n,'YTickLabel',data_names(s),'XLim',[-xl xl],'TickLabelInterpreter','latex','XTick',-xl:xl/2:xl,'TickDir','out')
    title('AV ($\pm$ 20\%)','Interpreter','latex')
    ax=gca;
    ax.XAxis.FontSize=8;
    ax.XAxis.Label.FontSize=10; 
    ax.YAxis.FontSize=8;
    text(xl*0.7,1,'(a)','Interpreter','latex','FontSize',10)
    
    tornado=(sens_std(:,[i1 i2])-sens_std(:,ref))./sens_std(:,ref)*100;
    [~,s]=sort(sum(abs(tornado),2),'ascend');
    xl=max(abs(tornado(:)));
    if xl<100
        xl=max([10 ceil(xl/10)*10]);
    else
        xl=ceil(xl/100)*100;
    end
    
    nexttile
    hold on
    barh(1:n,tornado(s,1))
    barh(1:n,tornado(s,2))
    set(gca,'YTick',1:n,'YTickLabel',data_names(s),'XLim',[-xl xl],'TickLabelInterpreter','latex','XTick',-xl:xl/2:xl,'TickDir','out')
    title('SD ($\pm$ 20\%)','Interpreter','latex')
    ax=gca;
    ax.XAxis.FontSize=8;
    ax.XAxis.Label.FontSize=10; 
    ax.YAxis.FontSize=8;
    text(xl*0.7,1,'(b)','Interpreter','latex','FontSize',10)

    tornado=(sens_infper(:,[i1 i2])-sens_infper(:,ref))./sens_infper(:,ref)*100;
    [~,s]=sort(sum(abs(tornado),2),'ascend');
    xl=max(abs(tornado(:)));
    if xl<100
        xl=max([10 ceil(xl/10)*10]);
    else
        xl=ceil(xl/100)*100;
    end

    nexttile
    hold on
    barh(1:n,tornado(s,1))
    barh(1:n,tornado(s,2))
    set(gca,'YTick',1:n,'YTickLabel',data_names(s),'XLim',[-xl xl],'TickLabelInterpreter','latex','XTick',-xl:xl/2:xl,'TickDir','out')
    if j<6
        switch j
            case 1
                xlabel('$\mathcal{E}$ sensitivity (\%)','Interpreter','latex')
            case 2
                xlabel('$\mathcal{A}_{\max}$ sensitivity (\%)','Interpreter','latex')
            case 3
                xlabel('$k_{\max}$ sensitivity (\%)','Interpreter','latex')
            case 4
                xlabel('$k_{\mathrm{end}}$ sensitivity (\%)','Interpreter','latex')
            case 5
                xlabel('$\mathcal{R}_*$ sensitivity (\%)','Interpreter','latex')
        end
    end
    title('$1 / \gamma$ ($\pm$ 20\%)','Interpreter','latex')
    ax=gca;
    ax.XAxis.FontSize=8;
    ax.XAxis.Label.FontSize=10; 
    ax.YAxis.FontSize=8;
    text(xl*0.7,1,'(c)','Interpreter','latex','FontSize',10)

    if j==6
        tornado=(sens_R0(:,[i1 i2])-sens_R0(:,ref))./sens_R0(:,ref)*100;
        [~,s]=sort(sum(abs(tornado),2),'ascend');
        xl=max(abs(tornado(:)));
        if xl<100
            xl=max([10 ceil(xl/10)*10]);
        else
            xl=ceil(xl/100)*100;
        end
            
        nexttile
        hold on
        barh(1:n,tornado(s,1))
        barh(1:n,tornado(s,2))
        set(gca,'YTick',1:n,'YTickLabel',data_names(s),'XLim',[-xl xl],'TickLabelInterpreter','latex','XTick',-xl:xl/2:xl,'TickDir','out')
        xlabel('HET sensitivity (\%)','Interpreter','latex')
        title('$\mathcal{R}_0$ ($\pm$ 20\%)','Interpreter','latex')   
        ax=gca;
        ax.XAxis.FontSize=8;
        ax.XAxis.Label.FontSize=10; 
        ax.YAxis.FontSize=8;
        text(xl*0.7,1,'(d)','Interpreter','latex','FontSize',10)

        nexttile(3)
        xlabel('HET sensitivity (\%)','Interpreter','latex')
        ax=gca;
        ax.XAxis.Label.FontSize=10; 
    end

    switch j
        case 1
            set(gcf,'NumberTitle','off','Name','Figure S4')
        case 2
            set(gcf,'NumberTitle','off','Name','Figure S5')
        case 3
            set(gcf,'NumberTitle','off','Name','Figure S6')
        case 4
            set(gcf,'NumberTitle','off','Name','Figure S7')
        case 5
            set(gcf,'NumberTitle','off','Name','Figure S8')
        case 6
            set(gcf,'NumberTitle','off','Name','Figure S9')
    end
    
end

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

function [R1_an,R2_an]=threshold(f,sigma0,sigma)

    R1_an=1/sigma0*(1-sigma)/max(f);
    R2_an=1/sigma0*sqrt((1-sigma^2)/sum(f.^2));

end

function ae=eval_ampenv(L,wind,nrm)

    ae=NaN(1,wind+1);
    ae(1)=1;
    for j=1:wind
        ae(j+1)=norm(L^j,nrm);
    end

end

function wind=find_duration(L,nrm) 

    wind=ceil(fzero(@(x) norm(L^ceil(x),nrm)-1,[1 1e4]));

end

