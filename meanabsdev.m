function MAD = meanabsdev(Sobs, Steo, Sn)
    if nargin<3
        Sn = zeros(size(Steo));
    end
    
    MAD = mean( abs( (Sobs./(Steo+Sn)) - mean(Sobs./(Steo+Sn)) )  );