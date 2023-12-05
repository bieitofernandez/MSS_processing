function SK=Tspec(spc,a,b,c)
%     q=3.4;
%     k=1.4e-7;
%     x=K / KB * sqrt(2*q);
%     f=x.*(exp(-x.^2/2)-x.*sqrt(pi()/2).*erfc(x/sqrt(2)));
%     SB=(q/2)^(1/2).*Xi/KB /k .*f;
%     SB=2*pi()*SB;
    %     kw=K/KB;
%     x=sqrt(2*q)*kw;
%     phik2=sqrt(q/2).*x.*exp(-x.^2/2)-x.^2.*sqrt(pi()/2).*erfc(x/sqrt(2));
%     SB=phik2./kw.^2;
%     SB=SB*Xi/(k*KB^3);

%Roget 1
% q=7.5;
% k=1.4e-7;
% y=K/KB*sqrt(q);
% 
% SK=exp(-sqrt(6)*y)./y*Xi*q^(3/2)/(k*KB^3);

if nargin==4
    Xi=a;
    KB=b;%*2*pi();
    K=c;%*2*pi();
elseif nargin==3
   Xi=a(1);
    KB=a(2);%*2*pi();
    K=b;%*2*pi();   
end
 
if spc=='K'    
    %Roget 2
    q=5.26;
    k=1.44e-7;
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*exp(-sqrt(3)*y); %Kraichnan
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
elseif spc=='B'
    %Roget 2
    q=3.9;
    k=1.44e-7;
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*(exp(-y.^2/2)-y.*sqrt(pi()/2).*(1-erf(y/sqrt(2)))); %Batchelor
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
end
SK = SK;%/2/pi(); %hai algo mal en todo o tema dos 2pi, creo que así
                %é o correcto
end
