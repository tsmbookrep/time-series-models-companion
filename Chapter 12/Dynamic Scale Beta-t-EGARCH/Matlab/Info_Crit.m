function [BIC,AIC] = Info_Crit (Logl,par,n)

% Logl=-Logl;

BIC=-2*Logl+size(par,1)*log(n);

AIC=2*size(par,1)-2*Logl;

end