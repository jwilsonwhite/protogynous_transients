function beta_fert_curve

% Figure to illustrate some different possible shapes of the fertilization
% curve

% Results from Pang's data (beta parameter):
% BEG = 40.88
% BBG = 200 (effectively Inf)

figure(11)
clf
set(gcf,'position',[200 500 400 300])

hold on

R = linspace(0,0.5,100);
a = [1 1 1];
b = [1.2823, 40.88, 20];
LS = {'-','--','-.'};

for bb = 1:length(b)
    
    plot(R,betacdf(R*2,a(bb),b(bb)),'linestyle',LS{bb},'color','k')
    
end

set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:0.1:1,'ytick',0:0.2:1);
set(gca,'xgrid','on')
set(gca,'xcolor','k','ycolor','k')
set(gca,'fontsize',12)
ylabel('Fertilization rate','fontsize',14)
xlabel('Sex ratio (proportion male)','fontsize',14)
axis square




