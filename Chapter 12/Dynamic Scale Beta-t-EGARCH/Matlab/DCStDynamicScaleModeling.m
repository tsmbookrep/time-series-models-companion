clear all
close all

% Resolve all folders relative to this script, so the repository can be run from any working directory.
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
dataDir = fullfile(scriptDir,'..','Data');
resultsDir = fullfile(scriptDir,'..','Results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This code sets up three models for modelling scale and degrees of freedom
%(df) in time series data using the Beta-t-EGARCH framework of Harvey (2013)

%The Beta-t-EGARCH model is a model for heavy tailed time series data modelled with a 
%student's t distributions which assumes a score-driven dynamics for either 
%the scale parameter \sigma, the df parameter \vega, or both.

%The estimation of the parameters of the models are through Maximum
%Likelihood

%A score-driven parameter for the scale parameter \sigma assumes a
%location scale model structure such as

% y_t=\mu+\sigma_t*\epsilon_t , where \epsilon_t is conditionally distributed with a student t distribution with parameters (0,1,\vega) scale is 1.

%then \sigma_t=exp(\lambda_t) where \lambda_t has a score-diven dynamics
%such as

% \lambda_t+1=\omega_s*(1-\phi_s)+\phi_s*\lambda_t+\kappa_s*u_st (This is an AR(1) specification)

%where u_st is the score with respect to \lambda_t of the conditional
%density of y_t, u_st=d f(y_t|\lambda_t,\vega) / d \lambda_t *S_t, where
%S_t can be wither 1 or the information matrix with respect to \lambda_t

%The recursion is initialised at \lambda_1=\omega_s

%THESE FOLLOWING PARTS ARE NOT IMPLEMENTED!!!!!!!!!!!!!!!

%A score-driven parameter for the scale parameter \vega assumes a
%location scale model structure such as

% y_t=\mu+\sigma*\epsilon_t , where \epsilon_t is conditionally distributed with a student t distribution with parameters (0,1,\vega_t) scale is 1.

%then \vega_t=exp(\upsilon_t) where \upsilon_t has a score-diven dynamics
%such as

% \upsilon_t+1=\omega_v*(1-\phi_v)+\phi_v*\upsilon_t+\kappa_v*u_vt (This is an AR(1) specification)

%where u_vt is the score with respect to \upsilon_t of the conditional
%density of y_t, u_vt=d f(y_t|\sigma,\vega_t) / d \upsilon_t *V_t, where
%V_t can be wither 1 or the information matrix with respect to \upsilon_t

%The recursion is initialised at \upsilon_1=\omega_v

%A model with both the parameters time varying is a combination of both the
%specifications simultaneously.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This parts loads the dtime series data

% X_=xlsread(fullfile(dataDir,'FTSE Close---.xlsx')); %FTSE
% X_=xlsread(fullfile(dataDir,'DJ Close---.xlsx')); %Dow Jones
X_=xlsread(fullfile(dataDir,'NSX Close---.xlsx')); %NSX
% X_=xlsread(fullfile(dataDir,'Hang Seng Close.xlsx'));
% X_=xlsread(fullfile(dataDir,'Nikkei 225 Close.xlsx'));
% X_=xlsread(fullfile(dataDir,'DAX Close.xlsx'));
% X_=xlsread(fullfile(dataDir,'CAC 40 Close.xlsx'));
% X_=xlsread(fullfile(dataDir,'S&P 100 Close.xlsx'));

%THis following part Set up price data in terms of returns
ss=size(X_,1);
X_1=X_([2:ss],:);
X_2=X_([1:ss-1],:);
%Setup as standard returns
% X=(X_1./X_2)-1;
%Setup as log returns
X=log(X_1./X_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This part defines the optimization options for the optimization
%algorithms, either fmincon for restricted optimization, or fminsearch for
%unrestricted

conoptions = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','MaxIter',1000*67,'MaxFunEvals',3000*67,'TolFun',0.00000001,'TolX',0.0000000000000001,'TolCon',0.00000000000001)
% 'TolX',0.0000000000000001    'TolFun',0.00000001
optionsunc = optimoptions('fminunc','Algorithm','quasi-newton','TolX',0.00001 ,'TolFun',0.00001, 'Display', 'iter-detailed', 'MaxIter', 1000*67, 'MaxFunEvals', 3000*67, 'HessUpdate', 'bfgs')

options = optimset('TolX',0.0000000000000001 ,'TolFun',0.00000001, 'Display', 'iter', 'MaxIter', 1000*67, 'MaxFunEvals', 3000*67, 'HessUpdate', 'bfgs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% options = optimset('TolX', 0.001, 'Display', 'off', 'MaxIter', 100000*67, 'MaxFunEvals', 30000*67, 'HessUpdate', 'bfgs') ;

%Set of controls set up the modelling specification for the scale of the series in the DCS-t-EGARCH
%If set to 1 the dynamics is of a randomw walk with phi=1, I(1), otherwise
%is an AR(1)
I_s=0;
%If this is 1 we standardise the score by its information matrix
Inf_s=1;
%if this is set to 1 we add a leverage effect into the scale dynamics as,
% +kappa_l*sign(-y_t)(u_t+1)
lev_s=1;

%Choose the optimization algorithm for dynamic scale, Constr=1 fmincon, Constr=0
%fminsearch
Constr_s=1;

%Choose the starting values for the kappa parameter for the optimnization
kk_s=0.01;

%THESE FOLLOWING PARTS ARE NOT IMPLEMENTED!!!!!!!!!!!!!!!

%Set of controls set up the modelling specification for the df of the series in the DCS-t-EGARCH
%If set to 1 the dynamics is of a randomw walk with phi=1, I(1), otherwise
%is an AR(1)
I_v=0;
%If this is 1 we standardise the score by its information matrix
Inf_v=1;
%If the below is set to 1, before etimating the df model we standardised
%the data by the previously fitted dynamic scale
Stan=0;

%Set of controls set up the modelling specification for both the scale and the df of the series in the DCS-t-EGARCH
%If set to 1 the dynamics of the scale is of a randomw walk with phi=1, I(1), otherwise
%is an AR(1)
I_sv_s=0;
%If this is 1 we standardise the score with respect to scale by its information matrix
Inf_sv_s=1;
%if this is set to 1 we add a leverage effect into the scale dynamics as,
% +kappa_l*sign(-y_t)(u_t+1)
lev_sv_s=0;
%If set to 1 the dynamics of the df is of a randomw walk with phi=1, I(1), otherwise
%is an AR(1)
I_sv_v=0;
%If this is 1 we standardise the score with respect to df by its information matrix
Inf_sv_v=1;

%Optimization algorithm for dynamic scale, Constr=1 fmincon, Constr=0
%fminsearch
Constr_v=0;
Constr_sv=0;

%Perform test if Tst=1
Tst=1;

%Choose the starting values for the kappa parameter for the optimnization
kk_v=0.9; %0.9 for standardized data after fitting Egarch scale
kk_sv=0.8; %0.7 better so far

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t Egarch WITH DYNAMIC SCALE

% mu = par(1); %mean
% omega = par(2); %unconstrained mean lambda scale
% phi = par(3); %dinamic AR parameter scale
% kapa= par(4); %dinamic cond score par
% vega_b = par(5); %log shape parameter (Could be the inverse of shape parameterbetween 0 and 1, maybe easier)

Y=X;

%This gives you a plot of the unconditional distribution of the data with
%respect to a Gaussian standard normal and few ACF figures
figure
%For a t distribution
[ff_x,r_x]=ecdf(Y);
ecdfhist(ff_x,r_x,100);
hold on
xx_x=min(Y):.001:max(Y);
yyn_x=1/sqrt(2*pi*std(Y)^2).*exp(-1/2*((xx_x-mean(Y))/std(Y)).^2);
plot(xx_x,yyn_x,'g-')
hold off
title('Unconditional Distr of the Data')

figure
autocorr(Y,100)
title('ACF Data')

figure
autocorr(Y.^2,100)
title('ACF Squared Data')

figure
qqplot(Y)
title('Empirical Density Data vs Gaussian Distr QQPlot')

% Y=Y(:,1);
% Y=cat(1,Ret1,Ret2);

%Size of the dataset
t=size(Y,1);

%Some initial values for the optimization
if I_s==1
    
    if lev_s==1
        %This is we include also leverage
        rS=ones(5,1)*kk_s;
    %     rS(1,1)=0.02;
        rS(1,1)=mean(Y);
        rS(2,1)=log(std(Y));
        rS(4,1)=log(8); %initial values of the df around 8
        
        if Constr_s==1
            %Setting up constrains in the parameters
            lbS=-inf(5,1);
            ubS=inf(5,1);
            
            lbS(3,1)=0; %kapa dynamic score
%             lbS(5,1)=0; %kapa leverage, positive
        end
    
    else
        %if we don't include leverage
        rS=ones(4,1)*kk_s;
    %     rS(1,1)=0.02;
        rS(1,1)=mean(Y);
        rS(2,1)=log(std(Y));
        rS(4,1)=log(8);

        if Constr_s==1
            lbS=-inf(4,1);
            ubS=inf(4,1);

            lbS(3,1)=0; %kapa dynamic score
        end
    
    end
    
else

    if lev_s==1
        %This is we include also leverage
        rS=ones(6,1)*kk_s;
    %     rS(1,1)=0.02;
        rS(1,1)=mean(Y);
        rS(2,1)=log(std(Y));
        rS(3,1)=0.9;
        rS(5,1)=log(8);

        if Constr_s==1
            lbS=-inf(6,1);
            ubS=inf(6,1);

            ubS(3,1)=1; %phi stationary between -1 and 1
            lbS(3,1)=0; %set phi to be positive
            lbS(4,1)=0; %kapa dynamic score
%             lbS(6,1)=0; %kapa leverage, positive
        end
    
    else
        %if we don't include leverage
        rS=ones(5,1)*kk_s;
    %     rS(1,1)=0.02;
        rS(1,1)=mean(Y);
        rS(2,1)=log(std(Y));
        rS(3,1)=0.9;
        rS(5,1)=log(8);

        if Constr_s==1
            lbS=-inf(5,1);
            ubS=inf(5,1);

            ubS(3,1)=1; %phi stationary between -1 and 1
            lbS(3,1)=0; %set phi to be positive 
            lbS(4,1)=0; %kapa dynamic score
        end
    
    end
    
end

%Check the value of the logl at the initial parameter vector rS
Beta_univ_t_Egarch_Logl (rS,Y,Inf_s,I_s)

%Construct the Logl
f_s=@(x)Beta_univ_t_Egarch_Logl (x,Y,Inf_s,I_s)

%Optimise
if Constr_s==1
    [x_s,fval,exitflag,output,v_,w_,hessian] = fmincon(f_s,rS,[],[],[],[],lbS,ubS,[],conoptions)    
else
    [x_s,fval,exitflag,output] = fminsearch(f_s,rS,options)    
end

%Retrive and save the parameters
estPar=x_s;

if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
xlswrite(fullfile(resultsDir,'estPar.xlsx'),x_s);

% x_s=xlsread('estPar_x_s.xlsx'); %This is from standard returnd

%Obtain the fitted values from the estimation
[Logl_s,lam_s,res_s,u_s,Beta_s,fit_s]= Beta_univ_t_Egarch_Logl_fact (x_s,Y,Inf_s,I_s);
%Get the info criteria
[BIC_s,AIC_s]=Info_Crit (Logl_s,x_s,t);

% %This is wrong
% fitY=zeros(t,1);
% 
% fitY(1,1)=undif(1,1);
% 
% for i=1:t-1
%     
%     fitY(i+1,1)=fit(i,1)+fitY(i,1);
% end

Ysq=Y.^2;

par_s=x_s;

%These are the estimated parameters individually
if I_s==1
    
    mu_s = par_s(1); %mean
    omega_s = par_s(2); %unconstrained mean lambda scale
    kapa_s= par_s(3); %dinamic cond score par
    % vega_b_s = par_s(4); %inverse shape parameter, between  0 and 1
    vega_l_s = par_s(4); %log shape parameter
    
    if lev_s==1
        kapa_s_l=par_s(5);
    end
    
else

    mu_s = par_s(1); %mean
    omega_s = par_s(2); %unconstrained mean lambda scale
    phi_s = par_s(3); %dinamic AR parameter scale
    kapa_s= par_s(4); %dinamic cond score par
    % vega_b_s = par_s(5); %inverse shape parameter, between  0 and 1
    vega_l_s = par_s(5); %log shape parameter
    
    if lev_s==1
        kapa_s_l=par_s(6);
    end
    
end

% vega_s=1/vega_b_s;
vega_s=exp(vega_l_s);

%Fitted values for sigma over the data
figure
subplot(2,3,[1:3])
hold on
plot(Y-mu_s)
plot(exp(lam_s),'r')
box on
legend('$Y-\mu$ Returns','$\hat{\sigma}_{t|t-1}=\exp\left(\hat{\lambda}_{t|t-1}\right)$','Interpreter','latex')
title('Scale Fit - Returns')
subplot(2,3,[4:6])
hold on
plot(abs(Y-mu_s))
plot(exp(lam_s),'r')
box on
legend('$\left|Y-\mu\right|$ Returns','$\hat{\sigma}_{t|t-1}=\exp\left(\hat{\lambda}_{t|t-1}\right)$','Interpreter','latex')
title('Scale Fit - Abs Returns')

%Fitted distribution of the residuals
figure
[ff_s,r_s]=ecdf(res_s);
ecdfhist(ff_s,r_s,100);
hold on
xx_s=min(res_s):.01:max(res_s);
yy_s=gamma((vega_s+1)/2)/(gamma(vega_s/2)*sqrt(pi*vega_s))*(ones(size(xx_s,1),1)+(xx_s.^2)/vega_s).^(-(vega_s+1)/2);
yyn_s=1/sqrt(2*pi).*exp(-1/2*xx_s.^2);
plot(xx_s,yy_s,'r-')
plot(xx_s,yyn_s,'g-')
legend('Empirical Density Residuals','Fitted t','Standard Gaussian')
hold off
box on
title('Empirical Density Fitted Residuals vs Fitted t Distr')

%ACF of fitted quantities
figure
autocorr(res_s,100)
title('ACF Residuals')

figure
autocorr(res_s.^2,100)
title('ACF Squared Residuals')

figure
autocorr(u_s,100)
title('ACF Fitted Scores')

%Compute the probability integral transform (PIT) of the residuals, if they
%are uniform the distributional assumption of a t distribution is correct.
fs=@(e)gamma((vega_s+1)/2)/(gamma(vega_s/2)*sqrt(pi*vega_s))*(1+(e.^2)/vega_s).^(-(vega_s+1)/2);
pit_s=@(g)integral(fs,-inf,g);
PIT_s=zeros(size(res_s,1),1);
ord_res_s=sort(res_s);
for i=1:size(res_s,1)
    PIT_s(i)=pit_s(ord_res_s(i));
end
PIT_uni_comp=linspace(0,1,size(PIT_s,1));

figure 
hold on
plot(PIT_s)
plot(PIT_uni_comp,'k--')
box on
title('Ordered PIT DCS Model vs Uniform')
xlabel('Residuals')
ylabel('PIT')
hold off


