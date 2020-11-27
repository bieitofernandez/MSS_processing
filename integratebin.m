function intX=integratebin(pres0,x,pres)

    intX=nan(size(pres));

    nn=find(pres>=pres0(end),1,'first')-1;

    for i=1:nn
         jj=find(pres0>=pres(i) & pres0<=pres(i+1));
         x0=pres0(jj);
         y0=x(jj);
         
         intX(i)= 0.5*sum( (y0(2:end)+(y0(1:(end-1)))) .* (x0(2:end)-(x0(1:(end-1)))) );
         intX(i)=intX(i)/(max(x0)-min(x0));
         
    end
    
