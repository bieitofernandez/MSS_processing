function x0=merge_sensors(x1,x2,fact)

    if nargin<3
        fact=5;
    end
    
    x0=nan(size(x1));
    for i=1:length(x1)
       if (x2(i)>fact*x1(i) || x1(i)>fact*x2(i))
          x0(i)=min( [x1(i),x2(i)]); 
       else
           x0(i)=mean([x1(i),x2(i)]);
       end 
    end

end