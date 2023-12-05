function dydx=grad_II(x,y)
    dydx=nan(size(y));
    
    for i=2:length(y)
        dydx(i,1)=(y(i)-y(i-1))/(x(i)-x(i-1));
    end

end