function xr=response(x,N,tau,deltat)
    xr=nan(size(x));
    for n=1:length(x)-N
        ii=[n:n+N]';
        a=1/N+tau/deltat*(12*ii-6*(N+1)) / (N*(N^2-1));
        a=a/sum(a);
        xr(n+floor(N/2)+1)=sum(a.*x(ii));
    end
    xr(1:floor(N/2))=x(1:floor(N/2));
    xr(floor(N/2)+1:end)=x(floor(N/2)+1:end);
end


