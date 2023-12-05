function [Q,Gamma,K,reg] = Bouffard_model(eps,N2,T,Pr)
         
 	 visco = viscosity(T);
	 [a,b] = size(eps);
         Q=eps./ (visco.*N2);
         K=nan(a,b);
         reg=nan(a,b);
	 for i=1:b
	     for j=1:a
                 if Q(j,i)<10^(1/3)*Pr^(-1/2)
                     K(j,i)=1e-7;
                     reg(j,i)=1;
                 elseif Q(j,i)<(3*log(sqrt(Pr)))^2
                     K(j,i)=0.1*Pr^(-1/4)*visco(j,i)*Q(j,i)^(3/2);
                     reg(j,i)=2;
                 elseif Q(j,i)<100
                     K(j,i)=0.2*visco(j,i)*Q(j,i);
                     reg(j,i)=3;
                 elseif ~isnan(Q(j,i))
                     K(j,i)=2*visco(j,i)*Q(j,i)^0.5;
                     reg(j,i)=4;
                 end
                 
      end
      Gamma = K .*N2 ./eps;
  end 
  
   