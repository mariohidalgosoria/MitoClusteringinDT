%%fig2a cluster size distribution

clear
CS =0;
%LU=1.0;
LU=0.5;
Nruns=100;
for i=1:Nruns
    %gamma 2.1
    %name=sprintf('../results/NO_FUS_B/SubTreeCV11/STC.%dclsnap.txt',i);
    %gamma=2
    %name=sprintf('../results/NO_FUS_B/SubTreeCV12/STC.%dclsnap.txt',i);
    %gamma=1.5
    name=sprintf('../results/NO_FUS_B/SubTreeCV13/STC.%dclsnap.txt',i);
    %gamma=1.9
    %name=sprintf('../results/NO_FUS_B/SubTreeCV14/STC.%dclsnap.txt',i);
    CSS = load(name);
    CSS = CSS';
    CS = [CS; CSS];
end
CS(1,:)=[];
CS = CS/LU; 
hist(CS)
mean(CS)

MM=max(CS);
MM= int32(MM);

hM=histogram(CS,MM);
pM=histcounts(CS,MM,'Normalization','pdf');
binCentersM = hM.BinEdges + (hM.BinWidth/2);
binCentersM = binCentersM';
pM = pM';
binM = binCentersM(1:end-1);
MP = [binM,pM];
%hold on
loglog(binM,pM,'.r')
xlabel('i')
ylabel('Cluster size dist')
xlim([1 10])
hold off
%gamma=2.1
%save("../results/NO_FUS_B/SubTreeCV11/STCSIMCLSDG21.dat","MP","-ascii")
%gamma=2
%save("../results/NO_FUS_B/SubTreeCV12/STCSIMCLSDG20.dat","MP","-ascii")
%gamma=1.5
save("../results/NO_FUS_B/SubTreeCV13/STCSIMCLSDG15.dat","MP","-ascii")
%gamma=1.9
%save("../results/NO_FUS_B/SubTreeCV14/STCSIMCLSDG19.dat","MP","-ascii")
%% mean field 


%%%%%%%%%%%%functions taken from lenafabr%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../networktools/')
addpath('../../dendriticmito/dendriticTrees_public/')%% load in realistic tree structure
% TODO: check input data
load('../examples/MCFO_HSE1_ST67.mat')
% switch to using um
NT.scaleCoords(umperpx);
NT.setCumEdgeLen();
figure
NT.plotNetwork(struct('datatipindex',true));
% rescale radii so that the trunk matches the experimental tree trunk
% radius
rtrunkexpt = (mean(NT.edgewidth{trunkedge}(:,1))/2);
radii = radii*rtrunkexpt;
%save net file for simulations
NT.outputNetwork(radii,'../results/NO_FUS_B/MCFO_HSE1_ST67_MU.net',struct('WRITEEVS',true))
NT0 = copy(NT);
%check ma radii/trunk in um 1.965
max(radii)
Vtree=pi*sum(NT.edgelens.*radii'.^2)%687.5



%% get mean-field expected total mass
v = 0.6;
gamma = 2.1;
mitovol=0.5;
%kb = 0.010032;
%kp = 0.3205;
%Au = 0.067335;
%kb = 0.0089;
%kp = 0.225;
%Au = 0.0953;
%unitlen = 1;
kb = 0.0115;
kp = 0.3205;
Au = 0.0772;
unitlen = 0.5;

% compute motile mito density, first without 2*kp/v prefactor
rhoWtrunk1 = 1;
[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
rhoWvals = rhoWvals*2*kp/v; % density in real units (per um)


% alpha parameter on each edge
alpj = v*Au./(2*kb*radii.^gamma);

% get densities of each different size stationary cluster
[M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvals,alpj,100);

% get volume density, unitless
voldens = mitovol*(M1vals.*NT.edgelens')./(pi*NT.edgelens'.*radii.^2);

% metrics (for comparison to expt)
avgvoldens = sum(voldens.*NT.edgelens'.*radii.^2)/sum(NT.edgelens'.*radii.^2)
totmass = sum(M1vals.*(NT.edgelens'));
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize;

expectedgemass = M1vals.*(NT.edgelens');
totmass = sum(expectedgemass)

% mean field cluster size distribution
[ivals, dcs] = getClusterSizeDist2(NT,rhoWvals,alpj,MM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2b vol density per branch
%% load in many runs and get average density on each edge
figure
clear
hold off
addpath('../../networktools/')
%% load in network structure (include curved edge paths)
load('../examples/MCFO_HSE1_ST67.mat','img','NT','radii','trunkedge','umperpx')
NT.plotNetwork(struct('labels',true))

% switch to using um
NT.scaleCoords(umperpx);
NT.setCumEdgeLen();
NT.plotNetwork(struct('datatipindex',true));
% rescale radii so that the trunk matches the experimental tree trunk
% radius
rtrunkexpt = mean(NT.edgewidth{trunkedge}(:,1))/2;
radii = radii*rtrunkexpt;
%gamma=2.1
%dirname = '../results/NO_FUS_B/SubTreeCV11/';
%gamma=2.0
%dirname = '../results/NO_FUS_B/SubTreeCV12/';
%gamma=1.5
dirname = '../results/NO_FUS_B/SubTreeCV13/';
%gamma=1.9
%dirname = '../results/NO_FUS_B/SubTreeCV14/';

edgedensity = [];

%unitlen = 1;
unitlen = 0.5;
Nruns=100;
for fc = 1:Nruns
    filename = [dirname sprintf('STC.%dsnap.txt',fc)];
    snapshots = parseParticleSegSnapshots(filename);
    nsnaps(fc) = length(snapshots);
    for sc = 1:1:length(snapshots)
        [clustermass,edgemass,parentlist] = getMassDistrib(snapshots(sc),NT);    
        edgedensity(:,end+1) = (edgemass*0.5)/unitlen./(pi*NT.edgelens.*radii'.^2);
    end    
end
avgedgedensity = mean(edgedensity,2);
STEedgedensity = std(edgedensity')'/sqrt(size(edgedensity,2));

%% predicted mean field edge density
v = 0.6;
%gamma = 2.1;
%kb =0.0115;
%kp = 0.3205;
%Au = 0.0772;
%gamma = 2.0;
%kb =0.01;
%kp = 0.3191;
%Au = 0.0679;
gamma = 1.5;
kb =0.01;
kp = 0.3447;
Au = 0.0539;
%gamma = 1.9;
%kb =0.01;
%kp = 0.3205;
%Au = 0.0672;



unitlen = 0.5;
mitovol=0.5;

% compute motile mito density, first without 2*kp/v prefactor
rhoWtrunk1 = 1;
[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
rhoWvals = rhoWvals*2*kp/v; % density in real units (per um)


% alpha parameter on each edge
alpj = v*Au./(2*kb*radii.^gamma);

% get densities of each different size stationary cluster
[M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvals,alpj,100);

% get volume density, unitless
voldens = mitovol*(M1vals.*NT.edgelens')./(pi*NT.edgelens'.*radii.^2);

% metrics (for comparison to expt)
avgvoldens = sum(voldens.*NT.edgelens'.*radii.^2)/sum(NT.edgelens'.*radii.^2)
totmass = sum(M1vals.*(NT.edgelens'));
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize;

expectedgemass = M1vals.*(NT.edgelens');
totmass = sum(expectedgemass)

%% adjusted motile density
%% find corrections in motile density to avoid finite lenght effects
% coefficients for 4th order polynomial, strting with highest order term
a0 = v*Au./(2*kb*radii(trunkedge).^gamma);

%polycoeff = [unitlen*a0^3/2, 0, a0*unitlen/2 - a0*v/2/kb, a0*kp/kb + v/2/kb,-kp/kb];
polycoeff = [a0*kb*unitlen - v*a0/2, a0*kp + v/2, - kp]
rho1mtrunk = roots(polycoeff);
[~,ind] = min(abs(rho1mtrunk-2*kp/v));
rho1mtrunk = rho1mtrunk(ind);
trunkmotmass = rho1mtrunk*(NT.edgelens(trunkedge)-0.5*unitlen)
%x = rho1mtrunk
%test = (kp/kb - v*x/2/kb)*(1-a0*x) - a0*unitlen*x^2/2 - unitlen*a0^3/2*x^4

% alpha parameter on each edge
alpj = v*Au./(2*kb*radii.^gamma);

% get densities of each different size stationary cluster
rhoWvalscl = rhoWvals*v/2/kp.*rho1mtrunk;

[M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvalscl,alpj,100);

% get volume density, unitless
voldens = mitovol*(M1vals.*NT.edgelens')./(pi*NT.edgelens'.*radii.^2);

% metrics (for comparison to expt)
avgvoldens = sum(voldens.*NT.edgelens'.*radii.^2)/sum(NT.edgelens'.*radii.^2)
totmass = sum(M1vals.*(NT.edgelens'));
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize;

expectedgemass = M1vals.*(NT.edgelens');
totmass = sum(expectedgemass)



modeldensity = ((M1vals*0.5)./(pi*radii.^2))';

scatter(modeldensity,avgedgedensity)
hold all
xlist = linspace(0.09,0.28);
plot(xlist,xlist,'k')
hold off

xlabel('predicted c_j')
ylabel('simulated c_j')
VVDMMF = [modeldensity,avgedgedensity];
%gamma=2.1
%save("../results/NO_FUS_B/SubTreeCV11/LLVDMFM.dat","VVDMMF","-ascii")
%gamma=2.0
%save("../results/NO_FUS_B/SubTreeCV12/LLVDMFM.dat","VVDMMF","-ascii")
%gamma=1.5
save("../results/NO_FUS_B/SubTreeCV13/LLVDMFM.dat","VVDMMF","-ascii")
%gamma=1.9
%save("../results/NO_FUS_B/SubTreeCV14/LLVDMFM.dat","VVDMMF","-ascii")

%% fig 2c mean field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create mesh with sub-cells
%of the edge segments
figure
clear
hold off
addpath('../../networktools/')
%% load in network structure (include curved edge paths)
load('../examples/MCFO_HSE1_ST67.mat','img','NT','radii','trunkedge','umperpx')
% switch to using um
NT.scaleCoords(umperpx);
NT.setCumEdgeLen();
NT.plotNetwork(struct('datatipindex',true));
% rescale radii so that the trunk matches the experimental tree trunk
% radius
rtrunkexpt = mean(NT.edgewidth{trunkedge}(:,1))/2;
radii = radii*rtrunkexpt;
NT.plotNetwork(struct('labels',true))



M=zeros(NT.nedge,16);

count = 0;
m = min(NT.edgelens)/3;
for ec = 1:NT.nedge
	   sec = NT.edgelens(ec)/m;
           sec = floor(sec);
           scc = NT.edgelens(ec) / sec;
      for ecc=1:sec
	  M(count+ecc,1)=ec;
          M(count+ecc,2)=ecc;
          M(count+ecc,3)=(ecc-1)*scc;        
          M(count+ecc,4)= ecc*scc;
          
	
	 newpos = interp1(NT.cumedgelen{ec},NT.edgepath{ec},(ecc-1)*scc);
         M(count+ecc,5)=newpos(1,1);
         M(count+ecc,6)=newpos(1,2);
         %store connectivity
         M(count+ecc,7)=NT.edgenodes(ec,1);
         M(count+ecc,8)=NT.edgenodes(ec,2);
         M(count+ecc,9)=NT.nodenodes(M(count+ecc,7),1);
         M(count+ecc,10)=NT.nodenodes(M(count+ecc,7),2);
         M(count+ecc,11)=NT.nodenodes(M(count+ecc,7),3);
         M(count+ecc,12)=NT.nodeedges(M(count+ecc,7),1);
         M(count+ecc,13)=NT.nodeedges(M(count+ecc,7),2);
         M(count+ecc,14)=NT.nodeedges(M(count+ecc,7),3);
         M(count+ecc,15)=NT.degrees(M(count+ecc,7),1);
      end
      count = count + sec;
end

hold off
plotopt = struct('plotnodes',[]);
NT.plotNetwork(plotopt)
hold on
plot(M(:,5),M(:,6),'.k')

%second find an array for the predicted mean field mass for each site
v=0.6;
%gamma=2.1;
%kb=0.0115;
%Au=0.0772;
%kp=0.3205;
%gamma=2.0;
%kb=0.01;
%Au=0.0679;
%kp=0.3191;
gamma = 1.5;
kb =0.01;
kp = 0.3447;
Au = 0.0539;
%gamma = 1.9;
%kb =0.01;
%kp = 0.3205;
%Au = 0.0672;


mdp = zeros(length(M),1);


%%%%%%%%%%%%%%%%%OLD
rtrunk=radii(trunkedge);
rhoWtrunk1 = 1;
[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
alpj = v*Au./(2*kb*radii.^gamma);
[M0vals, M1vals] = getMitoClusterDensity(rhoWvals*2*kp/v,alpj);
%%check mass and cluster size
totmass = sum(M1vals.*(NT.edgelens'))
getTreeAvgClusterSize(NT,rhoWvals*2*kp/v,alpj)
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitlen = 0.5;
mitovol=0.5;

% compute motile mito density, first without 2*kp/v prefactor
rhoWtrunk1 = 1;
[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
rhoWvals = rhoWvals*2*kp/v; % density in real units (per um)


% alpha parameter on each edge
alpj = v*Au./(2*kb*radii.^gamma);

% get densities of each different size stationary cluster
[M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvals,alpj,100);

% get volume density, unitless
voldens = mitovol*(M1vals.*NT.edgelens')./(pi*NT.edgelens'.*radii.^2);

% metrics (for comparison to expt)
avgvoldens = sum(voldens.*NT.edgelens'.*radii.^2)/sum(NT.edgelens'.*radii.^2)
totmass = sum(M1vals.*(NT.edgelens'));
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize;

expectedgemass = M1vals.*(NT.edgelens');
totmass = sum(expectedgemass)

motilemassU=sum((2*kp/v)*rhoWvals.*NT.edgelens')
%fraction of motile monomers mass wr to the total mass, should be 
motilemassU/totmass
%at the primary trunk the fraction of motile mito should be 0.02
((2*kp/v)*rhoWvals(trunkedge)*NT.edgelens(trunkedge))/(M1vals(trunkedge)*NT.edgelens(trunkedge))

%% adjusted motile density
% coefficients for 4th order polynomial, strting with highest order term
a0 = v*Au./(2*kb*radii(trunkedge).^gamma);
%polycoeff = [unitlen*a0^3/2, 0, a0*unitlen/2 - a0*v/2/kb, a0*kp/kb + v/2/kb,-kp/kb];
polycoeff = [a0*kb*unitlen - v*a0/2, a0*kp + v/2, - kp]
rho1mtrunk = roots(polycoeff);
[~,ind] = min(abs(rho1mtrunk-2*kp/v));
rho1mtrunk = rho1mtrunk(ind);
trunkmotmass = rho1mtrunk*(NT.edgelens(trunkedge)-0.5*unitlen)
%x = rho1mtrunk
%test = (kp/kb - v*x/2/kb)*(1-a0*x) - a0*unitlen*x^2/2 - unitlen*a0^3/2*x^4

% alpha parameter on each edge
alpj = v*Au./(2*kb*radii.^gamma);

% get densities of each different size stationary cluster
rhoWvalscl = rhoWvals*v/2/kp.*rho1mtrunk;

[M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvalscl,alpj,100);

% get volume density, unitless
voldens = mitovol*(M1vals.*NT.edgelens')./(pi*NT.edgelens'.*radii.^2);

% metrics (for comparison to expt)
avgvoldens = sum(voldens.*NT.edgelens'.*radii.^2)/sum(NT.edgelens'.*radii.^2)
totmass = sum(M1vals.*(NT.edgelens'));
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize;

expectedgemass = M1vals.*(NT.edgelens');
totmass = sum(expectedgemass)


totmass = sum(M1vals.*(NT.edgelens'))
%should be 0.2
VMT= sum(M1vals.*NT.edgelens')*mitovol/Vtree
%should be 2
avgclustsize = totmass/sum(M0vals.*(NT.edgelens'))
expectnclust = totmass/avgclustsize
%get mass of motile monomers 
motilemassU=sum((2*kp/v)*rhoWvals.*NT.edgelens')
%fraction of motile monomers mass wr to the total mass, should be 0.02
motilemassU/totmass
%fraction of motile monomeres in the primary trunk
 (2*kp/v)*rhoWvals(trunkedge)*NT.edgelens(trunkedge)/M1vals(trunkedge)*NT.edgelens(trunkedge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:length(M)
mdp(l)=  (  M1vals(M(l,1)) )*NT.edgelens(M(l,1));
end
for l=1:length(M)
        sec = NT.edgelens(M(l,1))/m;
        sec = floor(sec);
        scc = NT.edgelens(M(l,1)) / sec;
	mdp(l)= mdp(l)/ceil(NT.edgelens(M(l,1))/scc);
end
for l=1:length(M)	
	mdp(l)= (mdp(l)*0.5)/(pi*(M(l,4)-M(l,3))*radii(M(l,1))^2);
end
hold off
plotopt = struct('plotnodes',[]);
NT.plotNetwork(plotopt)
hold on			  
scatter(M(:,5),M(:,6),40,mdp,'square','filled')
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap("parula");
h=colorbar;
set(h,'fontsize',14);
%title('Mean Field; \gamma=2.1','Position',[60,110],'FontSize',14)
%title('Mean Field; \gamma=2.0','Position',[60,110],'FontSize',14)
title('Mean Field; \gamma=1.5','Position',[60,110],'FontSize',14)
%title('Mean Field; \gamma=1.9','Position',[60,110],'FontSize',14)
ylabel(h,'volume density per site  [A.U.]','FontSize',16)
%clim([0.1 0.28])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2d simulations
figure
clear
hold off
addpath('../../networktools/')
%% load in network structure (include curved edge paths)
load('../examples/MCFO_HSE1_ST67.mat','img','NT','radii','trunkedge','umperpx')
% switch to using um
NT.scaleCoords(umperpx);
NT.setCumEdgeLen();
NT.plotNetwork(struct('datatipindex',true));
% rescale radii so that the trunk matches the experimental tree trunk
% radius
rtrunkexpt = mean(NT.edgewidth{trunkedge}(:,1))/2;
radii = radii*rtrunkexpt;
NT.plotNetwork(struct('labels',true))

%create mesh with sub-cells
%of the edge segments

M=zeros(NT.nedge,16);

count = 0;
m = min(NT.edgelens)/3;
for ec = 1:NT.nedge
	   sec = NT.edgelens(ec)/m;
           sec = floor(sec);
           scc = NT.edgelens(ec) / sec;
      for ecc=1:sec
	  M(count+ecc,1)=ec;
          M(count+ecc,2)=ecc;
          M(count+ecc,3)=(ecc-1)*scc;        
          M(count+ecc,4)= ecc*scc;
          
	
	 newpos = interp1(NT.cumedgelen{ec},NT.edgepath{ec},(ecc-1)*scc);
         M(count+ecc,5)=newpos(1,1);
         M(count+ecc,6)=newpos(1,2);
         %store connectivity
         M(count+ecc,7)=NT.edgenodes(ec,1);
         M(count+ecc,8)=NT.edgenodes(ec,2);
         M(count+ecc,9)=NT.nodenodes(M(count+ecc,7),1);
         M(count+ecc,10)=NT.nodenodes(M(count+ecc,7),2);
         M(count+ecc,11)=NT.nodenodes(M(count+ecc,7),3);
         M(count+ecc,12)=NT.nodeedges(M(count+ecc,7),1);
         M(count+ecc,13)=NT.nodeedges(M(count+ecc,7),2);
         M(count+ecc,14)=NT.nodeedges(M(count+ecc,7),3);
         M(count+ecc,15)=NT.degrees(M(count+ecc,7),1);
      end
      count = count + sec;
end

hold off
plotopt = struct('plotnodes',[]);
NT.plotNetwork(plotopt)
hold on
plot(M(:,5),M(:,6),'.k')

%%simualtions 
%LU=1;
NINT=100;
LU=0.5;
%NINT=60;
for k=1:NINT
   %gamma=2.1 
   %name=sprintf('../results/NO_FUS_B/SubTreeCV11/STC.%dsnap.txt',k);
   %gamma=2.0 
   %name=sprintf('../results/NO_FUS_B/SubTreeCV12/STC.%dsnap.txt',k);
   %gamma=1.5 
   name=sprintf('../results/NO_FUS_B/SubTreeCV13/STC.%dsnap.txt',k);
   %gamma=1.9 
   %name=sprintf('../results/NO_FUS_B/SubTreeCV14/STC.%dsnap.txt',k);
%% parse particle snapshots
   snapshots = parseParticleSegSnapshots(name);
   for tc =1:length(snapshots)
   grp = snapshots(tc);
   [segments,nsegpt] = interpolateParticleSegPos(grp,NT);

    for pc = 1:grp.npart 
      
	    for i =1:length(M)   
	     if (grp.edgeU(pc) == M(i,1) && grp.edgeD(pc) == M(i,1)) 	       
                 if (grp.edgeposU(pc) > M(i,3)  && grp.edgeposD(pc) < M(i,4))
                     M(i,16) = M(i,16) + (grp.edgeposD(pc) - grp.edgeposU(pc));
                     %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
	         elseif (grp.edgeposD(pc)>M(i,3) && grp.edgeposD(pc)< M(i,4))
	             M(i,16) = M(i,16) + (grp.edgeposD(pc) - M(i,3));
                     %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
	         elseif (grp.edgeposU(pc)>M(i,3) && grp.edgeposU(pc)< M(i,4))
	             M(i,16) = M(i,16) + (M(i,4) -grp.edgeposU(pc));
                     %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
                elseif (grp.edgeposU(pc)<M(i,3) && grp.edgeposD(pc)> M(i,4))
	            M(i,16) = M(i,16) + M(i,4)-M(i,3);
                     %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
             
                 end


	     elseif(grp.edgeD(pc) == M(i,1) && grp.edgeU(pc) == M(i,12))
	           if (grp.edgeposD(pc) < M(i,4) && grp.edgeposD(pc) > M(i,3) )
                       M(i,16) = M(i,16) + (grp.edgeposD(pc)-M(i,3));
                       %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
                elseif (grp.edgeposD(pc) > M(i,4)  )
                       M(i,16) = M(i,16) + M(i,4)-M(i,3);
                    %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
	        end
	     elseif(grp.edgeU(pc) == M(i,1) && (grp.edgeD(pc) == NT.nodeedges(M(i,8),2) || grp.edgeD(pc) == NT.nodeedges(M(i,8),3)  ))
	        
		    if (grp.edgeposU(pc) < M(i,4) && grp.edgeposU(pc) > M(i,3) )
		    M(i,16) = M(i,16) + (M(i,4)-grp.edgeposU(pc));
                       %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
                    elseif (grp.edgeposU(pc)<M(i,3) )
	             M(i,16) = M(i,16) + M(i,4)-M(i,3);
                       %M(i,16) = M(i,16)/(M(i,4)-M(i,3));
		    end  
         end    
	   
        end    
     
     end
   end
       
end

for ii=1:length(M)
       M(ii,16)=M(ii,16)/length(snapshots);
   end
%volumetric density
for i=1:length(M)	
	M(i,16)= (M(i,16)*0.5)/(pi*(M(i,4)-M(i,3))*radii(M(i,1))^2);
end
		
%%
figure
hold off
plotopt = struct('plotnodes',[]);
NT.plotNetwork(plotopt)
hold on			  
scatter(M(:,5),M(:,6),30,M(:,16)/(NINT*LU),'square','filled')
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap("parula");
h=colorbar;
set(h,'fontsize',14);
%title('Simulations; \gamma=2.1','Position',[70,110],'FontSize',14)
%title('Simulations; \gamma=2.0','Position',[70,110],'FontSize',14)
title('Simulations; \gamma=1.5','Position',[70,110],'FontSize',14)
%title('Simulations; \gamma=1.9','Position',[70,110],'FontSize',14)
ylabel(h,'volume density per site  [A.U.]','FontSize',16)
%clim([0.1 0.28])
%clim([0.1 0.22])
clim([0.1 0.3])
%clim([0.1 0.22])
