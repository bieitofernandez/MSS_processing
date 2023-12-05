function [X,sX]=pres_av(pres0,x,pres,pint,fact,show)
    if nargin<6
       show=false;
    elseif nargin<5
        fact=0;
    end
    X=nan(size(pres));
    sX=nan(size(pres));
    
    iend = find(isfinite(pres0),1,'last');
    nn=find(pres>=pres0(iend),1,'first')-1;
    if isempty(nn)
        nn = length(pres);
    end
    for i=1:nn-1
         x0=x(pres0>=pres(i)-0.5*pint & pres0<=pres(i)+0.5*pint);
         X(i)=nanmean(x0);
         sX(i)=nanstd(x0);
         
          if fact>0
            iiout=find(x0>X(i)+fact*sX(i) | x0<X(i)-fact*sX(i));
            x0(iiout)=NaN;
            X(i)=nanmean(x0);
            sX(i)=nanstd(x0);
            if show
                fprintf('\n Depth: %1.2f m, ndel/ntot: %d/%d',pres(i),numel(iiout),numel(x0));
            end
         end
    end

end