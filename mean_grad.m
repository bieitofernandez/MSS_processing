function grX=mean_grad(pres0,x,pres,pint)
    
    grX=nan(size(pres));
    sgrX=nan(size(pres));

    nn=find(pres>=pres0(end),1,'first')-1;
    if isempty(nn)
        nn = length(pres);
    end

    for i=1:nn
         jj=find(pres0>=pres(i)-0.5*pint & pres0<=pres(i)+0.5*pint);
         if length(jj)>512
             x0=pres0(jj);
             y0=x(jj);

             [grX(i),sgrX(i)] = slope(-x0,y0);
         end
         
    end
end