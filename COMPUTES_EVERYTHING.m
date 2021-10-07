function [OUTPUT,SPEC]=COMPUTES_EVERYTHING(filename,nL,p0,pint,pend,pintD,pintGR,sL,sOV,Tdis,shsens,fit_noise,plt)

    if nargin<13
        plt=0;
        if nargin<12
           fit_noise = 1;
           if nargin<11
                shsens=12;
                if nargin<10
                    Tdis='n';
                    if nargin<3
                        p0=11;
                        pint=1;
                        pintD=1;
                        pintGR = 1;
                        pend=300;
                    end
                end
           end
        end
    end

    if nargout >1
        save_spec = 1;
        K = logspace(0,3,1024);
        nK = length(K);
    else
        save_spec = 0;
    end
   

    K1=2;
    K2=14;%14;
    Kmax = 30;
    %%%%%%%%%%%%%%%%% OS POLINOMIOS CAMBIANSE EN DIS_SPEC %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=1.44e-7;
    fmaxT=60;%estaba en 40 cpm
    kcut = K1; %wavenumber cutoff for highpass filter
    fs = 1024; %sampling frequency 

    g=9.80665;
    
    if ischar(filename)
        data=importdata(filename,' ',nL); %nL = 50 Malaspina
    else
        data=filename;
    end
    date = data.textdata{3};
    
    pres0=data.data(:,3);
    T0=data.data(:,4);
    ntc0=data.data(:,2);
    COND0=data.data(:,5);
    pshear=data.data(:,12);
    fluo0 = data.data(:,6);
    turb0 = data.data(:,7);   
    
    %smooth pressure
    %pres0=move_av(pres0,1000);

    %defines the grid
    pres=transpose([p0:pint:pend]);


    %calculo salinidade
    S0=sw_salt(COND0./sw_c3515,T0,pres0);
    pden0=sw_pden(S0,T0,pres0,0);

    %Corrects temperature time response
    T0=response(T0,19,160,1);
    
    
    %thermodynamic variables
    alpha0=sw_alpha(S0,T0,pres0);
    beta0=sw_beta(S0,T0,pres0);

    
    %NTC (According to hatmut)
    %ntc0=move_av(ntc0,40);
    %grNTC0=grad_II(-pres0,ntc0);
    %grNTC0=move_av(grNTC0,20);

    %save FP07 for export
    OUTPUT.pres_raw = pres0;
    OUTPUT.T_FP07_raw = ntc0;
    
    if plt
       figure(2)
       clf
       subplot(2,3,1)
       plot(T0,-pres0)
       hold on
       plot(ntc0,-pres0)
       legend('CTD','FP07')
       xlabel('T (oC)')
       ylabel('pres (db)')
       
       subplot(2,3,2)
       plot(S0,-pres0)
       xlabel('S')
       yticklabels([])
       
       subplot(2,3,3)
       plot(move_av(grad_II(-pres0,ntc0),20),-pres0)
       xlabel('dTdz (K/m)')
       yticklabels([])
       
       
       subplot(2,3,4)
       plot(data.data(:,10), -pres0)
       xlabel('Sh1 (1/s)')
       ylabel('pres (db)')
       
       subplot(2,3,5)
       plot(data.data(:,11), -pres0)
       xlabel('Sh2 (1/s)')
       yticklabels([])
       
       subplot(2,3,6)
       plot(pshear, -pres0)
       xlabel('pSh (1/s)')
       yticklabels([])
       pause()
    end

     %sorting
    [T_sort0,isT]=sort(T0,1,'descend');
    [ntc_sort0,isNTC]=sort(ntc0,1,'descend');
    [pden_sort0,isd]=sort(pden0);
   
    disT=pres0-pres0(isT);
    disNTC=pres0-pres0(isNTC);
    disden=pres0-pres0(isd);

    %filters T 
    mW = (max(pres0)-min(pres0))/(length(pres0))*fs;
    fcut = kcut*mW;
    [b,a] = butter(1,fcut/(fs/2),'high');
    ntcf=filtfilt(b,a,ntc0);
    %pshear=filtfilt(b,a,pshear);

    
    T=pres_av(pres0,T0,pres,pint,2.7);
    S=pres_av(pres0,S0,pres,pint,2.7);
    pden=pres_av(pres0,pden_sort0,pres,pint,2.7);
    alpha=pres_av(pres0,alpha0,pres,pint,2.7);
    beta=pres_av(pres0,beta0,pres,pint,2.7);
    fluo = pres_av(pres0,fluo0,pres,pint,2.7);
    turb = pres_av(pres0,turb0,pres,pint,2.7);

    grT=mean_grad(pres0,T0,pres,pintGR);
    grS=mean_grad(pres0,S0,pres,pintGR);
    grpden=mean_grad(pres0,pden_sort0,pres,pintGR);
    
    W = nan(size(T));
    
    epsilon = nan(size(T));
    epsN = nan(size(T));
    peps = nan(size(T));
    varpsh = nan(size(T));
    if shsens ==12
        eps1 = nan(size(T));
        epsN1 = nan(size(T));
        MADs1 = nan(size(T));
        MADfs1 = nan(size(T)); 
        MADcs1 = nan(size(T));
        
        eps2 = nan(size(T));
        epsN2 = nan(size(T));
        MADs2 = nan(size(T));
        MADfs2 = nan(size(T)); 
        MADcs2 = nan(size(T));
    else
        MADs = nan(size(T));
        MADfs = nan(size(T)); 
        MADcs = nan(size(T));
    end
    
    Xi  = nan(size(T));
    Xif = nan(size(T));
    KBT = nan(size(T));
    Xiv = nan(size(T));
    epsT = nan(size(T));
    
    MADx = nan(size(T));
    MADxf = nan(size(T));
    X_fit_flag = nan(size(T));
    MADop=nan(size(pres));
    LKHR =nan(size(pres));
    maxKT =nan(size(pres));
    
    if save_spec
        Ssh1 = nan(length(pres), nK);
        Ssh10 = nan(length(pres), nK);
        Ssh2 = nan(length(pres), nK);
        Ssh20 = nan(length(pres), nK);
        SgrT = nan(length(pres), nK);
        SgrT0 = nan(length(pres), nK);
    end
    
    
    
   
    %reads shear and calculates epsilon for 1 or 2 sensors

    if shsens==12
        %if there are two sensors
        shear1=data.data(:,10);
        shear2=data.data(:,11);
        %high-passes
        %shear1=filtfilt(b,a,shear1);
        %shear2=filtfilt(b,a,shear2);
        
        l= find(pres-0.5*pintD>pres0(1),1,'first');

        while l<=length(pres)
            if pres(l)+0.5*pintD>max(pres0)
                break
            end
            j=find(pres0>=(pres(l)-0.5*pintD) & pres0<=(pres(l)+0.5*pintD) );

            W(l) = (pintD)/(length(j))*fs;
            if j(end)>length(pres0)
                break;
            end
            if pres0(j(end))>pres(end)+0.5*pintD
                break;
            end
            
            visco=viscosity(nanmean(T0(j)));
            
            %shear microstructure
            if plt ==1
                fprintf('\n Shear Sensor 1');
            end
            [eps1(l),~,epsN1(l),~,MADs1(l),MADfs1(l),MADcs1(l),SPECsh1]=dis_spec(pres0(j),shear1(j),pshear(j),K1,K2,Kmax,visco,sL, sOV, plt);
            if plt ==1
                fprintf('\n Shear Sensor 2');
            end
            [eps2(l),~,epsN2(l),~,MADs2(l),MADfs2(l),MADcs2(l),SPECsh2]=dis_spec(pres0(j), shear2(j),pshear(j),K1,K2,Kmax,visco,sL, sOV, plt);
            if plt ==1
            fprintf('\n Pseudoshear');
            end
            [peps(l),Kcp]=dis_spec(pres0(j),pshear(j),zeros(length(j)),K1,K2,Kmax,visco, sL, sOV, plt);
            varpsh(l) = var(pshear(j)); %pseudoshear variance
            
            %Test fits and merge sensors
            epsilon(l)=merge_sensors(eps1(l),eps2(l)); %Merge epsilon sensors
            epsN(l)=merge_sensors(epsN1(l),epsN2(l)); %Merge epsilon sensors
            

            %Temperature micro-structure    
            if Tdis~= 'n'
                  if plt ==1
                    fprintf('\n FP07 Temperature Sensor');
                   end
                  if sum(isfinite(ntcf(j)))==length(j) 
                      KB=1/(2*pi())*(epsilon(l)/(visco*k^2))^(1/4);   
                      [Xiv(l),Xi(l),Xif(l),KBT(l),X_fit_flag(l),~,~, MADx(l), MADxf(l), MADop(l),~,LKHR(l),maxKT(l),SPECT]=Xi_spec(pres0(j),ntcf(j),K1,fmaxT,KB,W(l),Tdis,sL, sOV, fit_noise,plt);
                      epsT(l) = visco*k^2*(2*pi()*KBT(l))^4;
                  end
            end
            l=l+1;

            %interpolates and saves spectra    
            if save_spec
                Ssh1(l,:) = interp1(SPECsh1.K, SPECsh1.S,K);
                Ssh10(l,:) = interp1(SPECsh1.K, SPECsh1.Snc,K);
                Ssh2(l,:) = interp1(SPECsh2.K, SPECsh2.S,K);
                Ssh20(l,:) = interp1(SPECsh2.K, SPECsh2.Snc,K);
                SgrT(l,:) = interp1(SPECT.K,SPECT.S,K);
                SgrT0(l,:) = interp1(SPECT.K,SPECT.Snc,K);
            end
            
        end
    else
        %if there only is one sensor
        if shsens==1
            shear0=data.data(:,10);
        elseif shsens==2
            shear0=data.data(:,11);
        end    
        %highpasses
        %shear0=filtfilt(b,a,shear0);
        
        l= find(pres-0.5*pintD>pres0(1),1,'first');
        while l<=length(pres)
            if pres(l)+0.5*pintD>max(pres0)
                break
            end
            
            j=find(pres0>=(pres(l)-0.5*pintD) & pres0<=(pres(l)+0.5*pintD) );
     
            W(l) = (pintD)/(length(j))*fs;
            
            if j(end)>length(pres0)
                break;
            end
            if pres0(j(end))>pres(end)+0.5*pintD
                break;
            end
            
            
            visco=viscosity(nanmean(T0(j)));
            
            dof = round(2*((length(j)-sL)/sOV +1));
            MADop(l) = (2/dof)^0.5;
            if plt ==1
            fprintf('\n Shear Sensor');
            end
            [epsilon(l),~,epsN(l),~,MADs(l),MADfs(l),MADcs(l), SPECsh]=dis_spec(pres0(j), shear0(j),pshear(j),K1,K2,Kmax,visco, sL, sOV, plt);
            %fprintf('\n Pseudoshear');
            if plt ==1
            fprintf('\n Pseudoshear');
            end
            [peps(l),Kcp]=dis_spec(pres0(j),pshear(j),zeros(length(j)), K1,K2,Kmax,visco, sL, sOV, plt);
            varpsh(l) = var(pshear(j)); %pseudoshear variance
  
            
                
            if Tdis~= 'n'
                 if plt ==1
                     fprintf('\n FP07 Temperature Sensor');
                 end
                 if sum(isfinite(ntcf(j)))==length(j) 
                     KB=1/(2*pi())*(epsilon(l)/(visco*k^2))^(1/4);  
                         [Xiv(l),Xi(l),Xif(l),KBT(l),X_fit_flag(l),~,~, MADx(l), MADxf(l), MADop(l),~,LKHR(l), maxKT(l), SPECT]=Xi_spec(pres0(j),ntcf(j),K1,fmaxT,KB,W(l),Tdis,sL, sOV, fit_noise, plt);
                         epsT(l) = visco*k^2*(2*pi()*KBT(l))^4;
                 end
            end
            l = l+1;
            
            if save_spec
                Ssh1(l,:) = interp1(SPECsh.K, SPECsh.S,K);
                Ssh10(l,:) = interp1(SPECsh.K, SPECsh.Snc,K);
                SgrT(l,:) = interp1(SPECT.K,SPECT.S,K);
                SgrT0(l,:) = interp1(SPECT.K,SPECT.Snc,K);
            end

        end
            
        
    end    
    
    N2=-g./pden.*grpden;
    Rden=alpha.*grT./ (beta.*grS); 
    Gamma_t=0.2;
    KOsb=Gamma_t*epsilon./N2;
    
    LTT=sqrt(pres_av(pres0,disT.^2,pres,pintD,2.7));
    LTden=sqrt(pres_av(pres0,disden.^2,pres,pintD,2.7));
    
    Gamma=N2.*Xi./ (2*epsilon.*grT.^2); 
    KT=Xi./(2*grT.^2);
    LTntc=sqrt(pres_av(pres0,disNTC.^2,pres,pintD,2.7));


    
    %OUTPUT variable
    OUTPUT.date = date;
    OUTPUT.pres = pres;
    OUTPUT.T = T;
    OUTPUT.S = S;
    OUTPUT.pden = pden;
    OUTPUT.N2 = N2;
    OUTPUT.grT = grT;
    OUTPUT.grS = grS;
    OUTPUT.Rden = Rden;
    OUTPUT.W = W;
    
    OUTPUT.epsilon = epsilon;
    OUTPUT.epsN = epsN;
    OUTPUT.peps = peps;
    OUTPUT.varpsh = varpsh;
    if shsens == 12
        OUTPUT.eps1 = eps1;
        OUTPUT.epsN1 = epsN1;
        OUTPUT.MADs1 = MADs1;
        OUTPUT.MADfs1 = MADfs1;
        OUTPUT.MADfops1 = MADcs1;

        OUTPUT.eps2 = eps2;
        OUTPUT.epsN2 = epsN2;
        OUTPUT.MADs2 = MADs2;
        OUTPUT.MADfs2 = MADfs2;
        OUTPUT.MADfops2 = MADcs2;

    else
        OUTPUT.MADs = MADs;
        OUTPUT.MADfs = MADfs;
        OUTPUT.MADfops = MADcs;

    end
 
    OUTPUT.Xi = Xi;
    OUTPUT.Xif = Xif;
    OUTPUT.Xiv = Xiv;
    OUTPUT.epsT = epsT;
    OUTPUT.MADx = MADx;
    OUTPUT.MADxf = MADxf;
    OUTPUT.MADxop = MADop;
    OUTPUT.LKHR = LKHR;
    OUTPUT.X_fit_flag = X_fit_flag;
    OUTPUT.KT = KT;
    OUTPUT.LTntc = LTntc;
    OUTPUT.KOsb = KOsb;
    OUTPUT.Gamma = Gamma;
    OUTPUT.LTden = LTden;
    OUTPUT.LTT = LTT;
    OUTPUT.LTuT = LTntc;
    OUTPUT.fluo = fluo;
    OUTPUT.turb = turb;
    
    if save_spec
       SPEC.K = K;
       SPEC.SgrT = SgrT;
       SPEC.SgrT0 = SgrT0;
       if shsens == 12
           SPEC.Ssh1 = Ssh1;
           SPEC.Ssh10 = Ssh10;
           SPEC.Ssh2 = Ssh2;
           SPEC.Ssh20 = Ssh20;
       else
           SPEC.Ssh = Ssh1;
           SPEC.Ssh0 = Ssh10;
       end
    end
    
end
