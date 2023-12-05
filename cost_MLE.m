function C11 = cost_MLE(Sobs, Steo, dof, Sn)
    if nargin<4
        Sn = 0;
    end
    if nargin<3
        dof = 6;
    end
    
    
    C11 = log( dof./(Steo + Sn+eps).*chi2pdf( dof.*Sobs./(Steo+Sn+eps) , dof ));
    C11 = - sum( C11(isfinite(C11)));
