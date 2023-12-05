function yd=detrendp(x0,y0,p)
%habia un erro nesta
yd=nan(size(y0));
pint=p(2)-p(1);
    for i=1:length(p)
        
        ii=find(x0>=p(i)-0.5*pint & x0<=p(i)+0.5*pint);
        
        x1=x0(ii);
        y1=y0(ii);

        y=y1(isfinite(y1));
        x=x1(isfinite(y1));
        
        n=numel(y);
        
        if n>50
        
            b=(n*sum(x.*y)-sum(x)*sum(y)) ./ (n*sum(x.^2)-sum(x)^2);
            a=mean(y)-b*mean(x);

            yaj=a+b*x1;
            yd(ii)=y1-yaj;
        end
    end


    
end