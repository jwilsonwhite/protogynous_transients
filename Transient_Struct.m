function Transient_Struct
Transient = struct([]); % Create empty structure
Transient_Params ; % Run parameter creation file
load('transient_params.mat')

savename = strcat('transient_feb2019.mat') ;
 

S = {'GON','SC1','SC2'}; % Scenarios
F = [0:0.1:2] ; % Fishing rate
PHI = [1:1:20] ; % values of the male importance parameter
Lfs = [10 20 30]; %Length at which fishing starts

for s = 1:length(S) 
    
for f = 1:length(F)
    Transient(1).(S{s}).F(f).F=F(f) ;

for l = 1:length(Lfs)
    Transient(1).(S{s}).F(f).Lf(l).Lf=Lfs(l) ;
  
for phi = 1:length(PHI)
    
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).PHI=PHI(phi) ;

    [Ntotal, FishedNtotal, AgeDist, SizeDist, BiomSRatio, NumSRatio, LambInit, Growth, Theta2, TimeConv, FertEggs] = Transient_Model(S{s},F(f),UA,PHI(phi),Lfs(l)) ;
    
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).Ntotal = Ntotal ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).FishedNtotal = FishedNtotal ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).AgeDist = AgeDist ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).SizeDist = SizeDist ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).BiomSRatio = BiomSRatio ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).NumSRatio = NumSRatio ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).LambInit = LambInit ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).Growth = Growth ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).Theta2 = Theta2 ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).TimeConv = TimeConv ;
    Transient(1).(S{s}).F(f).Lf(l).PHI(phi).FertEggs = FertEggs ;


    
    
end
end
end
end

save(savename)



