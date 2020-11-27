function [Xiv,XiC,Xif,KBT,fit_flag, sXif,sKBT,MAD, MADf,MADc,MLKH,LKHratio,maxK,SPEC]  = Xi_spec(pres,x0,K1,fmax,KB,W,Tdis, sL, sOV,fit_noise,plt)

    Xiv = NaN;XiC=NaN;Xif=NaN;KBT=NaN;sXif = nan; sKBT = nan; MAD = NaN; MADf = NaN;MLKH=nan; LKHratio = nan; maxK = nan; fit_flag = 0;
    Pr = 7.56;
    DT = 1.44e-7;
    if nargin<8
        plt=0;
    end
    if nargin<7
        fit_noise = 1;
    end
    %reduces number of points
    %x0 = x0(1:4:end);
    %pres = pres(1:4:end);
    
    %detrend
    x = detrend(pres,x0,1);
    I=find(isfinite(x));
    x=x(I);
    pres=pres(I);

    mpres=mean(pres);

    Fs=length(pres)./ (max(pres)-min(pres));
    
    if isempty(x) | sum(x==0)==length(pres)
        return;
    end
    
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

    %calculates spectrum
    [PSDT,K]=pwelch(x-nanmean(x),hann(sL),sOV,sL,Fs,'onesided');
    PSD = (2*pi()*K).^2.*PSDT; %gradient spectra
    fr = W*K;
    Kn = fmax/W;
     
    iK1=find(K>K1,1,'first');
    iKn0=find(K>=Kn,1,'first');
    iKB=find(K>=KB,1,'first'); if isempty(iKB), iKB=length(K); end
    
    %correccion for time response
    PSDnc = PSD;
    %Goto 2016
    %tau = 0.005*W.^-0.32;
    %(Gregg and Meagher, 1995)/ Peterson & Fer2014
    tau0 = 0.012;
    tau = tau0*W.^-0.32;
    %simple
    %tau = 0.012;
    H = 1./(1+(2*pi()*tau*fr).^2).^2;
    PSD = PSDnc./H;
    
    %output for spectra
    SPEC.K = K;
    SPEC.S = PSD;
    SPEC.Snc = PSDnc;
    SPEC.dof = dof;
     
    if fit_noise == 1
        %fits noise function
         noiseC = false;
         Kni = Kn;
         Knf = min([Kni*5,max(K)]);
         fac = 0;
         noise0 = [-9.5,-0.1];
         while ~noiseC
             try
                 pnoise=nlinfit(fr(K>Kni & K<Knf),log(PSDT(K>Kni & K<Knf)),@(params0,fr0)log(FP07noise(params0,fr0)),noise0+[fac,0]);
                 noiseC = true;
             end
             fac = fac + 0.25;
         end
         %pnoise
    else
        %takes a fixed value (with the new spectral calculation works well)
        pnoise = [-9.5,-0.1];
    end
    Sn = FP07noise(pnoise,fr);
    Sn = Sn.*(2*pi()*K).^2;
    Sn = Sn./H;
    Sn(1) = 0;

%     figure(33)
%     clf
%     loglog(K,PSD,'k')
%     hold on
%     loglog(K, Sn)
%     pause()
    %Deletes undesired part of the spectrum
    iKn = min([iKn0,find(H<0.01,1,'first')]); 
    maxK = K(iKn);
    iKnM0 = find(PSD(iK1:end)<2*Sn(iK1:end),1,'first')+iK1-1;

    %first guess
    cont = false;
    ikcor = iK1:iKnM0;
    Xiv =6*DT*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    Xin = 6*DT*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    if Xiv>Xin
        cont = true;
    else
        Xiv = 0;
    end
    Xiv = Xiv -Xin;

    %continues only if Xiv is detectable
    if cont
         %calculates Xi using theoretical spectrum to correct for unresolved
         BAT=Tspec(Tdis,Xiv,KB,K);
         XiT=6*DT*sum( (BAT(ikcor(2:end))+BAT(ikcor(1:end-1))).*(K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
         XiC=Xiv.*Xiv./XiT; 
         BAT=Tspec(Tdis,XiC,KB,K);
         MAD = meanabsdev( PSD(ikcor), BAT(ikcor), Sn(ikcor) );
 
         %fits parameters in the noise free region
         KF = 2*K(iKnM0);
         Ktest = linspace(max([5*K1,KF/5]),KF*5,40);
         cost = nan(size(Ktest));
         Xif0 = nan(size(Ktest));
         dK = Ktest(2)-Ktest(1);
         for i = 1:length(Ktest)
             ks = 0.04*Pr^(-0.5)*Ktest(i);
             Kf1 = max([ks,K1]);
             iKf1 = find(K>=Kf1,1,'first');
             ikfit = iKf1:iKn;
             [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, dof, Sn);
         end

         Ktest = Ktest(isfinite(cost));
         cost = cost(isfinite(cost));
         LKHtest = - cost;

         ML = max(LKHtest);
         iML = find(LKHtest == ML);
         KB00 = Ktest(iML);
         ks = 0.04*Pr^(-0.5)*KB00;
         Kf1 = max([ks,K1]);    
         iKf1 = find(K>=Kf1,1,'first');
         ikfit = iKf1:iKn;
         MLp = -cost_T_fit(K, PSD, KB00+dK, ikfit, Tdis, dof, Sn);
         MLm = -cost_T_fit(K, PSD, KB00-dK, ikfit, Tdis, dof, Sn);
         deltak0=abs((2*dK)/sqrt(2*ML-MLm-MLp));
         deltak = max([deltak0,dK]);

         kmin2 = max([KB00-deltak,K(iK1+1)]);
         %kmax2 = min([KB00+deltak, (K(iKn)-3*dK)/0.04*Pr^0.5 ])
         kmax2 = min([KB00+deltak,5*K(iKn)]);
         Ktest = linspace(kmin2,kmax2,40);
         clear cost Xif0
         cost = nan(size(Ktest));
         Xif0 = nan(size(Ktest));
         for i = 1:length(Ktest)
             ks = 0.04*Pr^(-0.5)*Ktest(i);
             Kf1 = max([ks,K1]);
             iKf1 = find(K>=Kf1,1,'first');
             ikfit = iKf1:iKn;
             [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, dof, Sn);
         end

         Ktest = Ktest(isfinite(cost));
         Xif0 = Xif0(isfinite(cost));
         cost = cost(isfinite(cost));
         LKHtest = - cost;
         MLKH = max(LKHtest);
         iML = find(LKHtest == MLKH);
         Xif = Xif0(iML);
         KBT = Ktest(iML);
         BATf = Tspec(Tdis,Xif,KBT,K);
         ks = 0.04*Pr^(-0.5)*KBT;
         Kf1 = max([ks,K1]);
         iKf1 = find(K>=Kf1,1,'first');
         ikfit = iKf1:iKn;
         MADf = meanabsdev( PSD(ikfit), BATf(ikfit), Sn(ikfit) );

         %calculates uncertainties in the fitting parameters
         [MLp,Xip] = cost_T_fit(K, PSD, KBT+deltak, ikfit, Tdis, dof, Sn);
         [MLm,Xim] = cost_T_fit(K, PSD, KBT-deltak, ikfit, Tdis, dof, Sn);
         MLp = - MLp;
         MLm = - MLm;
         sKBT=abs((2*deltak)/sqrt(2*MLKH-MLm-MLp));
         sXif = abs((Xip-Xim)/sqrt(2*MLKH-MLm-MLp));

         %fits to polynom (avoiding the noisy part)
         ikfitA = ikfit(1):min([ikfit(end),find(BATf<2*Sn,1,'first')]);
         logK = log(K(ikfitA));
         logS = log(PSD(ikfitA));
         pp=polyfit(logK, logS,1);
         Sm = exp(polyval(pp, log(K)));
         LKHpol = -cost_MLE(PSD(ikfit), Sm(ikfit), dof, Sn(ikfit));

         LKHratio = (MLKH - LKHpol);
         LKHratio = log10(exp(1))*LKHratio;
         %likelyhood ratio is log10(Pteo/Pexp) and C = log(Pteo), to converto


         if LKHratio>0 && MADf<MADc*2 && KBT<3*maxK && abs(sKBT)<0.5*abs(KBT)
             fit_flag =1;
         end
     end

        %plots
        iKBT=find(K>=KBT,1,'first');if isempty(iKBT), iKBT=length(K); end

        if plt~=0
            figure(1)
            clf
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [14 9]); 
            set(gcf, 'PaperPositionMode', 'manual') ;
            set(gcf, 'PaperPosition', [0 0 14 9]);
            axes('position',[0.1,0.15,0.2,0.8])
            plot(x0-mean(x0),pres,'linestyle','--','color',[0.5,0.5,0.5])
            hold on
            plot(x,pres,'k')
            axis ij
            grid('on')
            xlabel('T'' (^oC)')
            ylabel('p (db)')
 
            axes('position',[0.4,0.15,0.55,0.8])
            loglog(K,PSD,'k')
            hold on
            %loglog(K2,PSD2,'b')
            loglog(K,PSDnc, 'color' ,[0.5,0.5,0.5])
            if cont
                loglog(K(ikfit),PSD(ikfit),'ko', 'markersize',4)%,'markerfacecolor','k')
                loglog(K(ikfit),Sm(ikfit)+Sn(ikfit),'b') 
                loglog(K,Sn+BATf, 'color',[0.5,0.5,0.5], 'linewidth',2);
                loglog(K,BATf, 'color',[0.5,0.5,0.5], 'linewidth',1,'linestyle','--');
                %loglog(K, Sn, 'r')
                loglog(K,Sn+BAT,'color','k', 'linewidth',2);
                loglog(K,BAT,'color','k', 'linewidth',1, 'linestyle','--');
                line([K(iKBT) K(iKBT)], ylim,'color','k')
            end
            loglog(K,Sn, 'color','r', 'linewidth',1,'linestyle','--');
            ylim([10^-9,1e3])
            line([K(iKn) K(iKn)], ylim,'color',[0.5,0.5,0.5],'linestyle','-')
            line([K(iK1) K(iK1)], ylim,'color',[0.5,0.5,0.5],'linestyle','-')
            
            text(1.1*KBT,0.1,'K_B')
            if Tdis=='B'
                text(1.2,10^2.5,'Batchelor Spectrum','horizontalalignment','left')
            elseif Tdis=='K'
                 text(1.2,10^2.5,'Kraichnan Spectrum','horizontalalignment','left')
            end           

            text(1.2,10^1.5,['\chi_{var. cor.} = ', num2str(XiC,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^0.75,['\chi_{fit} = ', num2str(Xif,'%1.3e'),'\pm',num2str(sXif,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^-0,['\chi_{var} = ', num2str(Xiv,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^-0.75,['K_B^{sh}= ', num2str(KB,'%1.0f'),' cpm'],'horizontalalignment','left')
            text(1.2,10^-1.5,['K_{B}= ', num2str(KBT,'%1.0f'),'\pm',num2str(sKBT,'%1.0f'),' cpm'],'horizontalalignment','left')
            text(100, 1e-7, ['MADf = ', num2str(MADf, '%1.2f'),'(', num2str(MADc*2, '%1.2f'),')'],'horizontalalignment','left')
            text(100, 1e-8, ['MAD = ', num2str(MAD, '%1.2f')],'horizontalalignment','left')
            text(100, 1e-6, ['LKHratio = ', num2str(LKHratio, '%1.1f')],'horizontalalignment','left')
            title([num2str(mpres), 'm'],'Fontsize',12)
            xlabel('K (cpm)')
            ylabel('PSD (K^2 m^{-2} cpm^{-1})')
            xlim([1,1000])
            fit_flag

            %text(100,10^-5,['K_B^{Batchelor} = ', num2str(KBT,'%1.0f'),' cpm'],'horizontalalignment','left')
            %saveas(gcf,'Tspec.png')
            pause()
            

        end

end
