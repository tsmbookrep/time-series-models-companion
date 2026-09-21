function [Logl,mu,res]= Bernoulli_Logl(par,y,tsmc,trend)

if tsmc==1
    
    omega = par(1); %unconstrained mean 
    kapa= par(3); %dinamic cond score par
  
else
    
    omega = par(1); %unconstrained mean 
    phi = par(2); %dinamic AR parameter scale
    kapa= par(3); %dinamic cond score par

end


T=size(y,1);

mu=zeros(T,1);
res=zeros(T,1);

%start AR iterations for dynamic location
mu_=mean(y);

for i=1:T
   
	mu(i)=mu_;
	res(i)=(y(i)-mu_); %the score is simply the residual

    if tsmc==1
            
        mu_=omega+kapa*y(i);

    else
        
        mu_= omega + phi*mu_ + kapa*res(i);
        

    end

end

Logl = -sum(y.*log(mu) + (1-y).*log(1-mu));





   



