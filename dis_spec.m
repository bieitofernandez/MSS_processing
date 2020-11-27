function [eps,Kc,epsN,KcN,MAD,MADf, MADc,SPEC]=dis_spec(pres,x0,px0,K1,K2,Kn,visco, sL, sOV, plt)
    eps=NaN;Kc=NaN;epsN=NaN;KcN=NaN;
    
    if nargin<6
        plt=0;
    end
    
    %deletes nan and detrend
    I = isfinite(x0);
    x0 = x0(I);
    px0 = px0(I);
    pres = pres(I);
    x = detrend(pres,x0);
    px = detrend(pres,px0);


    I = isfinite(x);
    x = x(I);
    px = px(I);
    pres = pres(I);
    
    
    if isempty(x) |  sum(x==0)==length(pres)
        return;
    end
    
    Fs=length(pres)./ (max(pres)-min(pres));

    %if sL = 0, sL = length(x)/2
    if sL == 0
        sL = round(length(x)/2)-1;
    end
    if sOV == 0
        sOV = round(sL/2);
    end
    
    %degrees of freedom acording to me
    %Nseg = floor((length(x)-sOV)/(sL-sOV));
    %dof = 2*Nseg;
    %more complex (methods in oceanography book)
    NsegM = floor(length(x)/(sL/2));
    dof = 0.92*8/3*NsegM;
    MADc = sqrt(2/dof);
    
    [PSD0,K] = pwelch(x-nanmean(x),hann(sL),sOV,sL,Fs,'onesided');
    if sum(abs(px))>0
        [PSDps, ~] = pwelch(px-nanmean(px),hann(sL),sOV,sL,Fs,'onesided');
        [CPS, ~] = cpsd(x-nanmean(x),px-nanmean(px),hann(sL),sOV,sL,Fs,'onesided');
        H = CPS./PSDps;
        PSDcont = abs(H).^2.*PSDps;
        PSD = PSD0 - PSDcont;
    else
        PSD = PSD0;
    end

    %output for spectra
    SPEC.K = K;
    SPEC.S = PSD;
    SPEC.Snc = PSD0;
    SPEC.dof = dof;

    
    dK = K(2) - K(1);
    mpres=mean(pres);
         
    iK1=find(K>=K1,1,'first');
    iK2=find(K>=K2,1,'first');
    iKn=find(K>=Kn,1,'first');
    
    %EPSILON ITERATIVE CALCULATION
    iK3=iK2;
    K3=K(iK3);
    Kc0=K2;
    flag=0;
    while flag==0

        eps=7.5*visco*sum( (PSD(iK1+1:iK3)+PSD(iK1:iK3-1)).* (K(iK1+1:iK3)-K(iK1:iK3-1)) )/2;
        Kc=1/(2*pi())*(eps./visco^3).^(1/4);
        if abs(Kc-Kc0)<=2*dK
            flag = 1;
            iK3 = find(K>=Kc,1,'first');
        elseif Kc<=Kn
            inc = (Kc-K3)/abs(Kc-K3);
            iK3= iK3+inc; %find(K >= Kc,1,'first');
            K3=K(iK3);
            Kc0 = Kc;
        else
            iK3=iKn;
            flag = 2;
        end
              
    end

    eps=7.5*visco*sum( (PSD(iK1+1:iK3)+PSD(iK1:iK3-1)).* (K(iK1+1:iK3)-K(iK1:iK3-1)) )/2;
    Kc=1/(2*pi())*(eps./visco^3).^(1/4); 

    
     % THE ITERATIVE CALCULATION HAS TO BE CORRECTED FOR LOST VARIANCE
     % (either assuming NasmYth or with defined polynoms)
     
     %(OPTION 1) correction for lost  of variance (assuming NasmYsth)
     %first corrects for the time-response
     p1 = [0.00245744,0.06590227,1.627074,2.099292];
     log_eps=log10(eps);
     log_eps=polyval( p1 , log_eps);
     eps=10.^log_eps;
     %now for lost variance
     NAS=nasmyth(eps,visco,K);
     varianceN=sum( (NAS(iK1+1:iK3)+NAS(iK1:iK3-1)).* (K(iK1+1:iK3)-K(iK1:iK3-1)) )/2;
     epsUN=7.5*visco*varianceN;
     if abs(eps-epsUN)/epsUN>0.05 
        eps=eps*eps/epsUN;
        Kc=1/(2*pi())*(eps./visco^3).^(1/4);       
     end

     
     
    %(OPTION 2)correction for lost variance (with polynom)
     %        %Correction polynom new (Malaspina)
    % $$$         p1 = [0.00245744,0.06590227,1.627074,2.099292];
    %  $$$         p2 = [-0.006687726,-0.1359909,0.08578571,-2.017015];
    % $$$         p3 = [0.005100397,0.1503782,2.462239,4.690549];
    % $$$        %Correction polynom old (Trynitrop)
    % $$$        p1 = fliplr([1.621873,1.506306,0.05198044,0.00175014]);
    % $$$        p2 = fliplr([4.690549,2.462239,0.1503782,0.005100397]);
    % $$$        p3 = fliplr([-2.017015,0.08578571,-0.1359909,-0.006687726]);
       %PNS06: desde malaspina en adiante. Non sei por que os p2 e p3
        %estaban cambiados
%         p1 = [0.00245744,0.06590227,1.627074,2.099292];
%         p2 = [0.005100397,0.1503782,2.462239,4.690549];
%         p3 = [-0.006687726,-0.1359909,0.08578571,-2.017015];
% 
%        log_eps=log10(eps);
%        log_eps=polyval( p1 , log_eps); %corrects spatial response for PNS06
%        log_eps=polyval( p2 , log_eps); %variance below 2 cpm
%        log_eps=polyval( p3 , log_eps); %corrects variance above 30 cpm
%        eps=10.^log_eps;
%         
       
       
       %EPSILON CALCULATION BY FITTING TO NASMYTH
       epsN=nlinfit(K(iK1:iKn),log(PSD(iK1:iKn)),@(e,kk)log(nasmyth(e,visco,kk)),10^-9);
       KcN=1/(2*pi())*(eps./visco^3).^(1/4);
       
       K3=K(iK3);
       
       NAS2=nasmyth(epsN,visco,K);
       NAS=nasmyth(eps,visco,K);

       MAD =  mean(abs(log10(PSD(iK1:iK3)./NAS(iK1:iK3))));
       MADf = mean(abs(log10(PSD(iK1:iKn)./NAS(iK1:iKn))));
       
       if plt~=0
            figure(1)
            clf
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [14 9]); 
            set(gcf, 'PaperPositionMode', 'manual') ;
            set(gcf, 'PaperPosition', [0 0 14 9]);
            axes('position',[0.1,0.15,0.2,0.8])
            hold on
            plot(x,pres,'k')
            axis ij
            grid('on')
            xlabel('shear (s^{-1})')
            ylabel('p (db)')
            
            axes('position',[0.4,0.15,0.55,0.8])
            loglog(logspace(0,3,100),nasmyth(logspace(-11,-1,11),visco,logspace(0,3,100)),'color',[0.60 0.60 0.60])
            hold on
            for i=1:11
               text(1.5,min(nasmyth(10^(-12+i),visco,[1.5 2])),['10^{',num2str(-12+i),'} W/kg'],'color',[0.60 0.60 0.60])
            end
            loglog(K,PSD0,'color',[0.5,0.5,0.5])
            %loglog(Kr,PSDr,'b')
            loglog(K,PSD,'k')    
            loglog(K(1:iKn),PSD(1:iKn),'ko', 'markersize',4)
            loglog(K,NAS,'k','linewidth',2)
            loglog(K,NAS2,'--k','linewidth',2)
            ylim([10^-9,10])
            xlim([1 10^3])
            line([Kc,Kc],ylim,'color','k')
            text(1.1*Kc,0.1,'K_c')
            line([Kn,Kn],ylim,'color','k','linestyle','--')
            text(1.2,10^-0.5,['\epsilon _{iter} = ', num2str(eps,'%1.3e'),' W/kg'],'horizontalalignment','left')
            text(1.2,10^-1,['\epsilon _{Nas}= ', num2str(epsN,'%1.3e'),' W/kg'],'horizontalalignment','left')
            text(2, 10^-8, ['MAD = ', num2str(MAD, '%1.2f'),'(', num2str(MADc*2, '%1.2f'),')'],'horizontalalignment','left')
            text(2, 10^-8.5, ['MADf = ', num2str(MADf, '%1.2f')],'horizontalalignment','left')
            title([num2str(mpres), 'm'])
            xlabel('k (cpm)')
            ylabel('S_{sh} (s^{-2} cpm^{-1})')
            pause()
       end
       

