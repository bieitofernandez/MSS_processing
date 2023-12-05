clc
clear
warning('off','all');


%addpath(genpath('./scripts_to_share/'))
addpath(genpath('/home/bieito/Documents/SCIENCE/FUNCTIONS_DATABASES/FUNCIONSMATLAB/'))

p1=2;
p2=500;

pint=1;
pintD=4;
pintGR=4;
%sL = 0;
sL=1024;
sOV=sL/2;
fit_noise = 1;
tau = 7/1000.;
temperature_spectrum = 'K';

nL=52;%header size

ND = round((p2-p1)/pint)+1;



%directory='/home/bieito/Documents/datos_remedios/';
directory='/media/bieito/Elements1/SCIENCE_DATA/BOCATS/3shear/'; 
%warning('off','MATLAB:dispatcher:InexactCaseMatch');

files=dir([directory,'she*.tob']);
% for i=1:length(files);
%     filename(i)=files(i).name;
% end
filename=files(1).name;

file = nan(length(files),12);
profN = nan(length(files),1);
T=nan(ND,length(files));
theta=nan(ND,length(files));
S=nan(ND,length(files));
sigma0=nan(ND,length(files));
N2=nan(ND,length(files));
N2_s=nan(ND,length(files));
grT=nan(ND,length(files));
grS=nan(ND,length(files));
eps=nan(ND,length(files));
eps1=nan(ND,length(files));
eps2=nan(ND,length(files));
Xic=nan(ND,length(files));
Rden=nan(ND,length(files));
Xif=nan(ND,length(files));
Xiv=nan(ND,length(files));
epsT=nan(ND,length(files));
LTntc=nan(ND,length(files));
Gamma=nan(ND,length(files));
Rif = nan(ND,length(files));
KT=nan(ND,length(files));
KOsb=nan(ND,length(files));
LTT=nan(ND,length(files));
LTd=nan(ND,length(files));
flu=nan(ND,length(files));
turb=nan(ND,length(files));
peps=nan(ND,length(files));
varpsh = nan(ND,length(files));
MAD=nan(ND,length(files));
MADf=nan(ND,length(files));
MADop=nan(ND,length(files));
LKHR=nan(ND,length(files));
eps_fit_flag=nan(ND,length(files));
Xi_fit_flag=nan(ND,length(files));

K = nan(1024,length(files));
SgrT = nan(ND,1024,length(files));
Ssh1 = nan(ND,1024,length(files));
Ssh2 = nan(ND,1024,length(files));

DATES = nan( length(files),6);
DATET = nan(length(files),1);

Corrects_Next = 0;


iFIRST =1;
for i=iFIRST:length(files)
    filename=files(i).name;
    filename
    file(i,:) = filename;
    profN(i,:) = str2num(filename(5:8));
    filename=[directory,filename];
    data=importdata(filename,' ',nL);
   
    a=data.textdata{3};
    i0 = strfind(a,',')+2;
    a =a(i0:end);
    if contains(a,'de julio de')
        i0 = strfind(a,'de julio de');
        a = strcat(a(1:i0-1),' July',a(i0+length('de julio de'):end));
    elseif contains(a,'de junio de')
        i0 = strfind(a,'de junio de');
        a = strcat(a(1:i0-1),' June',a(i0+length('de junio de'):end))
    end

    DATES(i,:) = datevec(a,'dd mmmm yyyy HH:MM:SS');
    
    DATET(i,1) = datenum(DATES(i,:));
    
  
    fprintf('\nProcessing File %i/%i, profile: %i, time: %s',i,length(files),profN(i,1),datestr(DATES(i,:)))
   
    [OUT,SPEC]=COMPUTES_EVERYTHING(data,nL,p1,pint,p2,pintD,pintGR,sL,sOV,'B',12,tau,fit_noise,0);
 
    T(:,i) = OUT.T;
    theta(:,i) = OUT.theta;
    S(:,i) = OUT.S;
    sigma0(:,i) = OUT.pden;
    grT(:,i) = OUT.grT;
    grS(:,i) = OUT.grS;
    N2(:,i) = OUT.N2;
    N2_s(:,i) = OUT.N2_s;
    Rden(:,i) = OUT.Rden;
    eps(:,i) = OUT.epsilon;
    eps1(:,i) = OUT.eps1;
    eps2(:,i) = OUT.eps2;
    Xic(:,i) = OUT.Xi;
    epsT(:,i) = OUT.epsT;
    Xif(:,i) = OUT.Xif;
    Xiv(:,i) = OUT.Xiv;
    LTntc(:,i) = OUT.LTntc;
    Gamma(:,i) = OUT.Gamma;
    Rif(:,i) = OUT.Rif;
    KT(:,i) = OUT.KT;
    KOsb(:,i) = OUT.KOsb;
    LTT(:,i) = OUT.LTT;
    LTd(:,i) = OUT.LTden;
    flu(:,i) = OUT.fluo;
    turb(:,i) = OUT.turb;
    peps(:,i) = OUT.peps;
    varpsh(:,i) = OUT.varpsh;
    MAD(:,i) = OUT.MADx;
    MADf(:,i) = OUT.MADxf;
    MADop(:,i) = OUT.MADxop;
    LKHR(:,i) = OUT.LKHR;
    Xi_fit_flag(:,i) = OUT.X_fit_flag;
    
    K(:,i) = SPEC.K;
    SgrT(:,:,i) = SPEC.SgrT;
    Ssh1(:,:,i) = SPEC.Ssh1;
    Ssh2(:,:,i) = SPEC.Ssh2;
         
     %save('procesando.mat');
end

pres=OUT.pres;
      

%grT=sqrt(2*KT./Xi);
[Reb,GammaBt,KBt,regBt]=Bouffard_model(eps,N2,T,7);    

 
%KT = Xif./(2*grT.^2);
%eps(:,6:7) = NaN;
%KOsb(:,6:7)= NaN;

save profiler_all_profiles_MATLAB_BOCATS_all.mat DATES DATET file profN pres T theta ...
    S sigma0 N2 N2_s grT grS Rden eps eps1 eps2 Xic Xiv Xif epsT peps LTntc Rif Gamma KT ...
    KOsb LTT LTd varpsh MAD MADf MADop LKHR Xi_fit_flag Reb GammaBt KBt regBt flu turb 

save spectra_SINES.mat DATES DATET file profN pres K SgrT Ssh1 Ssh2




