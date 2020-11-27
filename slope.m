function [b,sb,a,sa]=slope(x,y)

    i_fin=find(isfinite(y));
         
    x=x(i_fin);
    y=y(i_fin);
    
    %n=length(x);
    %b=(n*sum(x.*y)-sum(x)*sum(y)) ./ (n*sum(x.^2)-sum(x)^2);
    
   
     whichstats={'beta','covb'};
     stats=regstats(y,x,'linear',whichstats);
     b0=stats.beta;
     covb=stats.covb;
     sb0=sqrt(diag(covb));
     
     b=b0(2);
     sb=sb0(2);
     
     a=b0(1);
     sa=sb0(1);
end
