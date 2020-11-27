function mx=move_av(x,n)
    
    if mod(n,2)==0;
        n=n+1;
    end
        mx=nan(size(x));
    
    for j=1:length(x)-n
       mx(j+floor(n/2))=nanmean(x(j:j+n-1)); 
    end
   
end