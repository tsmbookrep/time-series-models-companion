clear all
close all

% Resolve all folders relative to this script, so the chapter runs from any working directory.
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
dataDir = fullfile(scriptDir,'..','Data');
resultsDir = fullfile(scriptDir,'..','Results');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This code sets up three models for modelling the Bernoulli score driven
%model as in Time Series Models: Theory and Applications (Harvey)


%The estimation of the parameters of the models are through Maximum
%Likelihood

%A score-driven parameter for the location/probability parameter \mu (\pi) assumes a
%location  model structure such as


% \mu_t+1=\omega_s*(1-\phi_s)+\phi_s*\mu_t+\kappa_s*u_mt (This is an AR(1) specification)

%where u_mt is the score with respect to \mu_t (\pi) of the conditional
%density of y_t

%The recursion is initialised at \lambda_1=\omega_s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This parts loads the boat race outcomes time series

X=xlsread(fullfile(dataDir,'BoatRace46.xlsx'));

for i=1:size(X,1)
    if isnan(X(i,1))
        X(i,1)=X(i-1,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This part defines the optimization options for the optimization
%algorithms, either fmincon for restricted optimization, or fminsearch for
%unrestricted

conoptions = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','MaxIter',1000*67,'MaxFunEvals',3000*67,'TolFun',0.00000001,'TolX',0.0000000000000001,'TolCon',0.00000000000001)
% 'TolX',0.0000000000000001    'TolFun',0.00000001
optionsunc = optimoptions('fminunc','Algorithm','quasi-newton','TolX',0.00001 ,'TolFun',0.00001, 'Display', 'iter-detailed', 'MaxIter', 1000*67, 'MaxFunEvals', 3000*67, 'HessUpdate', 'bfgs')

options = optimset('TolX',0.0000000000000001 ,'TolFun',0.00000001, 'Display', 'iter', 'MaxIter', 1000*67, 'MaxFunEvals', 3000*67, 'HessUpdate', 'bfgs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% options = optimset('TolX', 0.001, 'Display', 'off', 'MaxIter', 100000*67, 'MaxFunEvals', 30000*67, 'HessUpdate', 'bfgs') ;

%Set of controls set up the modelling specification for the location of in the DCS Bernoulli model
%If set to 1 simple Markov Chain model (no score)
tsmc=0;


Y=X;


%Size of the dataset
t=size(Y,1);

%Some initial values for the optimization
if tsmc==1

    omega_init = 0.1;

    kappa_init = 0.1;

    par = [omega_init; kappa_init];

    C = [1  0  0;    % delta
         1  0  1]    % delta + kappa
        
    A = [C;-C];
    
    b = [ones(2,1);
         zeros(2,1)];

    lb = [   -Inf; -1; -Inf];
    ub = [    Inf;  1;  Inf];

   
else

    omega_init = 0;
    phi_init = 0.99;
    kappa_init = 0.1;

    
    par = [omega_init; phi_init; kappa_init];

    C = [1  0  0;     % delta
         1  0  1;     % delta + kappa
         1  1  0;     % delta + phi
         1  1 -1];    % delta + phi - kappa
    
    
    % Previous constraints: tol <= C*par <= 1-tol
    A = [ C;
         -C];
    
    b = [ones(4,1);
         zeros(4,1)];

    lb = [   -Inf; -1; -Inf];
    ub = [    Inf;  1;  Inf];
    

end
    

%Check the value of the logl at the initial parameter vector rS
[Logl_start,mu_start,res_start]=Bernoulli_Logl(par,Y,tsmc)

%Construct the Logl
f_s=@(x) Bernoulli_Logl(x,Y,tsmc);

%Optimise

[x_s,fval,exitflag,output,v_,w_,hessian] =  fmincon(f_s, par, A, b, [], [], lb, ub, [], conoptions); 

%Retrive and save the parameters
estPar=x_s;

[Logl,mu,res]= Bernoulli_Logl(x_s,Y,tsmc)

xlswrite(fullfile(resultsDir,'estPar.xlsx'),estPar);


%Plot (for boat race example)

figure
years = 1946:2026; %change years if needed

plot(years, Y, '.', 'MarkerSize', 20);

hold on

plot(years, mu, '-', ...
    'LineWidth', 1.2, ...
    'Color', [0.8500 0.3250 0.0980]);


% Show a labelled tick every 10 years
xticks(1950:10:2020);
xlim([1946 2026]);
legend({'Winner (0 Oxf, 1 Cam)', 'Dyn. Prob. (AR1)'}, ...
       'Location', 'northeast', 'FontSize', 9);



xlabel('Year');
ylabel('Result and prob');




