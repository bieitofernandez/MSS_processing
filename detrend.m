function yd=detrend(x,y,order)

  if nargin<3
    order = 2;
  end

  yd=nan(size(y));
  n=length(x);
  %b=(n*sum(x.*y)-sum(x)*sum(y)) ./ (n*sum(x.^2)-sum(x)^2);
  %a=mean(y)-b*mean(x);

  %yaj=a+b*x;
  pol = polyfit(x,y,order);
  yaj = polyval(pol,x);
  yd=y-yaj;

    
end
