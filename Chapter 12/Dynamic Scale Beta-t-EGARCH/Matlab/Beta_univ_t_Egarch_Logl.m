function Logl= Beta_univ_t_Egarch_Logl (par,y,inf,I)

if I==1
    
    mu = par(1); %mean
    omega = par(2); %unconstrained mean lambda scale
    kapa= par(3); %dinamic cond score par
    % vega_b = par(4); %inverse shape parameter, between  0 and 1
    vega_l = par(4); %log shape parameter
    
    % vega=1/vega_b;
    vega=exp(vega_l);

    if size(par,1)==5
        kapa_2=par(5);
    end
    
else
    
    mu = par(1); %mean
    omega = par(2); %unconstrained mean lambda scale
    phi = par(3); %dinamic AR parameter scale
    kapa= par(4); %dinamic cond score par
    % vega_b = par(5); %inverse shape parameter, between  0 and 1
    vega_l = par(5); %log shape parameter

    % vega=1/vega_b;
    vega=exp(vega_l);

    if size(par,1)==6
        kapa_2=par(6);
    end
    
end

% nab_b = 0; %inverse shape parameter, between  0 and 1
% vega = 2; %other shape parameter, greater than 1

T=size(y,1);

lam=zeros(T,1);
res=zeros(T,1);
u=zeros(T,1);
Beta=zeros(T,1);

%set up info matrix for dynamic component (scale). 
if inf==1
    %if normalised for it
    Inf=2*vega/(vega+3);
else
    %if don't normalise
    Inf=1;
end

%start AR iterations for dynamic scale
lam_=omega;

for i=1:T
   
	lam(i)=lam_;
	res(i)=(y(i)-mu)*exp(-lam(i));

    Beta(i)=((res(i)^2)/vega)/(1+((res(i)^2)/vega));
    u(i)=((vega+1)*Beta(i)-1)/Inf;
    
    if I==1
        
        if size(par,1)==5
            lam_=lam(i)+kapa*u(i)+kapa_2*sign(mu-y(i))*(u(i)+1);
        else
            lam_=lam(i)+kapa*u(i);
        end
   	
    else
        
        if size(par,1)==6
            lam_=omega*(1-phi)+phi*lam(i)+kapa*u(i)+kapa_2*sign(mu-y(i))*(u(i)+1);
        else
            lam_=omega*(1-phi)+phi*lam(i)+kapa*u(i);
        end

    end

end

LoglresSum=log(ones(T,1)+1/vega*res.^2);

Logl=T*log(gamma((vega+1)/2)/(gamma(vega/2)*sqrt(pi*vega)))-sum(lam)-(vega+1)/2*sum(LoglresSum);

Logl=-Logl;




   



