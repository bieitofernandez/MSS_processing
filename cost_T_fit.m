function [C11, XiC] = cost_T_fit(K,Sobs, KB, iKfit, Tdis, dof, Sn)
    if nargin<7
        Sn = zeros(size(Sobs));
    end
    if nargin<6
        dof = 6;
    end
    
    iKB = find(K<=KB,1,'last');
    iKnoi= find(Sobs(iKfit(1):end)<2*Sn(iKfit(1):end),1,'first')+iKfit(1)-1;
    iKcorr = iKfit(1):min([iKB,iKnoi,iKfit(end)]);
	
    %initial Xi guess
    Xi0 = 6*1.4e-7*sum( (Sobs(iKcorr(2:end))+Sobs(iKcorr(1:end-1))).* (K((iKcorr(2:end)))-K((iKcorr(1:end-1)))) )/2;
    XiN0 = 6*1.4e-7*sum( (Sn(iKcorr(2:end))+Sn(iKcorr(1:end-1))).* (K(iKcorr(2:end))-K((iKcorr(1:end-1)))) )/2;
    Xi0 = Xi0- XiN0;
    
    %corrects for loss variance
    Steo = Tspec(Tdis, Xi0, KB, K );
    XiT = 6*1.4e-7*sum( (Steo(iKcorr(2:end))+Steo(iKcorr(1:end-1))).* (K((iKcorr(2:end)))-K((iKcorr(1:end-1)))) )/2;
    XiC = Xi0*Xi0/XiT;
    
    Steo = Tspec(Tdis, XiC, KB, K );
    C11 = cost_MLE(Sobs(iKfit), Steo(iKfit), dof, Sn(iKfit));
    %C11 = cost_MLE(Sobs, Steo, dof, Sn);
    
    
    
    
