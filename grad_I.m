
function dydx=grad_I(x,y)
    dydx=nan(size(y));
    
    for i=3:length(y)
        dydx(i,1)=(-y(i)+4*(i-1)-3*y(i-2))/(x(i)-x(i-2));
    end

end