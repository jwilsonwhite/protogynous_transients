function [Ntotal, FishedNtotal, AgeDist, SizeDist, BiomSRatio, NumSRatio, LambInit, Growth, Theta2, TimeConv, FertEggs] = Transient_Model(S,F,PHI,Lff)

load transient_params.mat
Lf = Lff;
% ------------------------------------------------------------------------
% Find Stable Age Distribution (SAD)

Z = UA ;
surv = exp(-Z) ;  % Instantaneous mortality rate... e^-M
Surv2 = repmat(surv/Lambda,[20,1]) ;

Ninit = [1; cumprod(Surv2(1:end-1))] ;
N(:,1) = Ninit ;


switch S 
% Different reproductive scenarios
%-------------------------------------------------------------------------
case 'GON'  %Gonochores
 
%-------------------------------------------------------------------------    
case 'SC1'   % Absolute Length
load transient_SC1_params.mat
% Probability of Maturity
M(:,1) = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Probability of Sex Change
SC(:,1) = 1./(1+exp(-p.*((D(2,:)')-Lc))) ; %prob. of being male

%------------------------------------------------------------------------
case 'SC2'   % Mean Length
load transient_SC2_params.mat
% Probability of Maturity
M(:,1) = 1./(1+exp(-q.*(D(2,:)'-(sum(D(2,:)'.*N(:,1)/sum(N(:,1)))+Lm)))) ;
% Probability of Sex Change
SC(:,1) = 1./(1+exp(-p.*(D(2,:)'-(sum(D(2,:)'.*N(:,1)/sum(N(:,1)))+Lc)))) ; %prob of being male

end % end switch

%-----------------------------------------------------------------------
switch S
    case 'GON' 
% Reproduction 
E(1)=sum(EggProd.*M_tmp.*N(:,1)./2) ; 
FertEggs(1) = E(1) ;
B2(1) = 0.5 ; 
% Calc. for instantaneous growth rate
X = FertEggs(1) ;
G = Lambda/X ;

    case {'SC1','SC2'} 
    
% Reproduction        
E(1) = sum((c.*D(2,:)'.^e).*N(:,1).*(1-SC(:,1)).*M(:,1)) ; % Egg production
TB(1) = sum((N(:,1).*M(:,1).*D(2,1:end)'.^3)) ; % Total mature biomass
MB(1) = sum((N(:,1).*SC(:,1).*M(:,1).*D(2,1:end)'.^3)) ; % Mature male biomass
B(1) = MB(1)/TB(1) ; % Biomass sex ratio
B2(1)= sum(N(:,1).*M(:,1).*SC(:,1))./sum(N(:,1).*M(:,1)); % Numerical sex ratio

Pf(1) = (betacdf(B(1),a,PHI)) ; % Proportion of eggs fertilized
FertEggs(1) = E(1)*Pf(1) ;

% Calc. for instantaneous growth rate
X = FertEggs(1) ;
G = Lambda/X ;

end % end switch



% ------------------------------------------------------------------------
% Start fishing population
TF = 50 ;
isfished(:,1) = 1./(1+exp(-r.*((D(2,:)')-Lf))) ;
Z = (UA+(F.*isfished)) ; % Annual adult survival per age class
surv = exp(-Z) ;  % Instantaneous mortality rate... e^-M

A = diag(surv(1:end-1))  ;
A = [A,zeros(size(A,1),1)]  ;
A = [zeros(1,size(A,2));A]  ; % Creates survival matrix

for tf = 2:TF

    N(:,tf) = A*N(:,tf-1) ;
    
%-----------------------------------------------------------------------
switch S
    case 'GON' ;
% Reproduction 
E(tf)=sum(EggProd.*M_tmp.*N(:,tf-1)./2) ; 
FertEggs(tf) = E(tf)*G ;
B(tf) = 1 ; % filler for sex ratio calcs
B2(tf) = 0.5 ;     

    case {'SC1','SC2'} ;
    
% Reproduction 
E(tf) = sum((c.*D(2,:)'.^e).*N(:,tf-1).*(1-SC(:,tf-1)).*M(:,tf-1)) ; % Egg production
TB(tf) = sum((N(:,tf-1).*M(:,tf-1).*D(2,1:end)'.^3)) ; % Total mature biomass
MB(tf) = sum((N(:,tf-1).*SC(:,tf-1).*M(:,tf-1).*D(2,1:end)'.^3)) ; % Mature male biomass
B(tf) = MB(tf-1)/TB(tf-1) ; % Biomass sex ratio
B2(tf)= sum(N(:,tf-1).*M(:,tf-1).*SC(:,tf-1))./sum(N(:,tf-1).*M(:,tf-1)); % Numerical sex ratio
Pf(tf) = (betacdf(B(tf),a,PHI)) ; % Proportion of eggs fertilized
FertEggs(tf) = E(tf)*Pf(tf)*G ;

end % end switch

% Recruits entering A0 
N(1,tf) = FertEggs(tf) ;

switch S 

%-------------------------------------------------------------------------
case 'GON'  %Gonochores
 
%-------------------------------------------------------------------------    
case 'SC1'   % Absolute Length
load transient_SC1_params.mat
% Probability of Maturity
M(:,tf) = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Probability of Sex Change
SC(:,tf) = 1./(1+exp(-p.*((D(2,:)')-Lc))) ; %prob. of being male

%------------------------------------------------------------------------
case 'SC2'   % Mean Length
load transient_SC2_params.mat
% Probability of Maturity
M(:,tf) = 1./(1+exp(-q.*(D(2,:)'-(sum(D(2,:)'.*N(:,tf)/sum(N(:,tf)))+Lm)))) ;
% Probability of Sex Change
SC(:,tf) = 1./(1+exp(-p.*(D(2,:)'-(sum(D(2,:)'.*N(:,tf)/sum(N(:,tf)))+Lc)))) ; %prob of being male

end % end switch

end

Npost= N(:,TF) ;
PostFishBiomSRatio = B(TF) ;
PostFishNumSRatio = B2(TF) ;   


%-----------------------------------------------------------------------
% Stop fishing
T = 50 ;

Z = (UA) ; % Annual adult survival per age class
surv = exp(-Z) ;  % Instantaneous mortality rate... e^-M
Surv = repmat(surv,[20,1]) ;

A = diag(Surv(1:end-1))  ;
A = [A,zeros(size(A,1),1)]  ;
A = [zeros(1,size(A,2));A]  ; % Creates survival matrix

for t = TF+1:(TF+T)
    
    N(:,t) = A*N(:,t-1) ;
    
%-----------------------------------------------------------------------
switch S
    case 'GON' ;
% Reproduction 
E(t)=sum(EggProd.*M_tmp.*N(:,t-1)./2) ; 
FertEggs(t) = E(t)*G ;
B(t) = 1.0 ; % filler for sex ratio calcs
B2(t) = 0.5 ;     

    
case {'SC1','SC2'} ;
    
% Reproduction        
E(t) = sum((c.*D(2,:)'.^e).*N(:,t-1).*(1-SC(:,t-1)).*M(:,t-1)) ; % Egg production
TB(t) = sum((N(:,t-1).*M(:,t-1).*D(2,1:end)'.^3)) ; % Total mature biomass
MB(t) = sum((N(:,t-1).*SC(:,t-1).*M(:,t-1).*D(2,1:end)'.^3)) ; % Mature male biomass
B(t) = MB(t-1)/TB(t-1) ; % Biomass sex ratio
B2(t)= sum(N(:,t-1).*M(:,t-1).*SC(:,t-1))./sum(N(:,t-1).*M(:,t-1)); % Numerical sex ratio

Pf(t) = (betacdf(B(t),a,PHI)) ; % Proportion of eggs fertilized
FertEggs(t) = E(t)*Pf(t)*G ;

end % end switch

% Recruits entering A0 
N(1,t) = FertEggs(t) ;

switch S 

%-------------------------------------------------------------------------
case 'GON'  %Gonochores

%-------------------------------------------------------------------------    
case 'SC1'   % Absolute Length
load transient_SC1_params.mat
% Probability of Maturity
M(:,t) = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Probability of Sex Change
SC(:,t) = 1./(1+exp(-p.*((D(2,:)')-Lc))) ; %prob. of being male

%------------------------------------------------------------------------
case 'SC2'   % Mean Length
load transient_SC2_params.mat
% Probability of Maturity
M(:,t) = 1./(1+exp(-q.*(D(2,:)'-(sum(D(2,:)'.*N(:,t)/sum(N(:,t)))+Lm)))) ;
% Probability of Sex Change
SC(:,t) = 1./(1+exp(-p.*(D(2,:)'-(sum(D(2,:)'.*N(:,t)/sum(N(:,t)))+Lc)))) ; %prob of being male


end % end switch


Growth(t) = (sum(N(:,t)))/(sum(N(:,t-1))) ;


SizeDist(:,t) = N(:,t).*D(2,:)' ;
BiomSRatio(t) = B(t) ; % Biomass sex ratio
NumSRatio(t) = B2(t) ; % Numerical sex ratio

end

% Calculate several transient statistics

Ntotal = sum(N) ;
T = 100 ;
for t = 1:T
    FishedNtotal(t) = sum(N(:,t).*isfished) ;
    
end

AgeDist = N ;
LambInit = (sum(N(:,51)))/(sum(N(:,50))) ; % Initial growth rate

A = N(:,1) ;
B = N(:,50) ;
CosTheta = dot(A,B)./norm(A)./norm(B) ;
Theta = acos(real(CosTheta)) ;
Theta2 = rad2deg(Theta) ; % Deviation from SAD


GR = Growth(51:100)' ;
GR2 = GR >= Lambda-0.001 & GR <= Lambda+0.001 ;

% Quick algorithm to find the time at which the population growth rate is at 99% of
% lambda
if isempty(find(flipud(GR2)==0,1,'first'))
    Conv = 1     ;
else Conv = (find(flipud(GR2)==0,1,'first')) ;
end

if Conv == 1  ;
    TimeConv = 1  ; 
else 
TimeConv = 50-Conv+1 ;
end



