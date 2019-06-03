function Transient_Params

% Set up parameters for transient protogynous population models

% General life history
Amax = 20; % Number of age classes
isfished = nan(Amax,1) ; % vector indicatign which age classes are in the fishery

% Natural mortality
UA = 0.35; % this is named UA for reasons that used to make but don't anymore.

% Reproduction (see Easter et al. 2016 Mar Ecol Prog Ser for sources)
a = 1 ;
c = 7.04 ;  % Constant in fecundity relationship
e = 2.95 ;  % Exponent in allometric relationship
K = 0.000003 ;  % Slope of fertilization fxn parameter
X = 0.09 ;  % Intercept of fertilization fxn parameter
EggProd = c.*(D(2,:)').^e ; % Potential egg production at each age
Lm = 20 ; % Length at which 50% of fish mature 
q = 1 ;  % Shape parameter in maturity function 
M_tmp = 1./(1+exp(-q.*((D(2,:)')-Lm)));   % Probability of Maturity 
Lambda = 1.02; % This is the fixed population growth rate in the closed population case

% von Bertalanffy Growth
D = zeros(2,Amax) ; 
D(1,:) = (0:(Amax-1)) ;
D(2,1) = 8 ;
k = 0.05 ;
Linf = 90 ;
T0 = -1.875; % Age at size 0 (gives it 8cm at size 0)
D(2,:) = Linf.*(1-exp(-k.*((0:(Amax-1))-T0))) ;



% Adult survival
r = 1 ;  % Steepness of selectivity curve (1,0.1)   
Lf = 30 ;  % Length at which 50% chance a fish will be removed

save transient_params.mat

clear all

% Scenario-specific values:
%--------------------------------------------------------------------------
% Sex-Change cue type 1 ('fixed') - Based on absolute length 


Lm = 20 ; % STANDARD Length at which 50% of fish mature 
Lc = 30 ;

p = 1 ;
q = 1 ;  % Shape parameter in maturity function 

save transient_SC1_params.mat

clear all

%--------------------------------------------------------------------------
% Sex-Change cue type 2 ('flexible') - Based on population length
% distribution

% Difference from the mean size at which prob. of maturity is 0.5
Lm = 4 ;
% Difference from the mean size at which prob. sex change is 0.5
Lc = 14 ;

p = 1.0 ;  % Shape parameter in sex change fxn
q = 1.0 ;  % Shape parameter in maturity function

save transient_SC2_params.mat

clear all

%--------------------------------------------------------------------------
