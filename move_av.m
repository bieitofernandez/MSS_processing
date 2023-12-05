function mx=move_av(x,n)
    
    if mod(n,2)==0;
        n=n+1;
    end
        mx=-1*ones(size(x));
    
    for j=1:length(x)-n
       mx(j+floor(n/2))=nanmean(x(j:j+n-1)); 
    end

    mx(1:floor(n/2)) = mx(1+floor(n/2));
    mx(length(x)-floor(n/2):end) = mx(length(x)-floor(n/2)-1);
end
