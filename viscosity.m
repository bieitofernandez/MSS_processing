function v=viscosity(T)
%seawater viscosity
v=1.792747-0.052126103*T+0.0005918645.*T.^2;
v=v*10^-6;
end