%Matlab script to generate asymmetry verus gamma
clear
hold off
addpath('../../networktools/')
addpath('../../dendriticmito/dendriticTrees_public/')
%% load in all dendritic tree networks
dirname = '/storage-server/shareddata/dendriticmito/cleanedNetWorkspaces/';
% load pixel calibration data
calibdata = readtable([dirname 'pixel_calibrations.csv'])
% load tree data
files = dir([dirname '*.mat'])
filenames = {files.name};

clear allNetworks networkNames origFileName usenetwork
allradii = {};
for fc = 1:length(filenames)
    %disp(filenames{fc})
    
    load([dirname filenames{fc}],'NT','parentnode','origimgfile');
    origFileNames{fc} = origimgfile;

    for pc = 1:height(calibdata)        
        if contains(origimgfile,calibdata.Cell{pc})
            disp(sprintf('%s: %s, %f',origimgfile,calibdata.Cell{pc},calibdata.Scale_micron_pixel_(pc)))
            umperpx(fc) = calibdata.Scale_micron_pixel_(pc);

            networkNames{fc} = calibdata.Cell{pc};            
        end
    end

    % scale the network to be in um units
    
    NT.scaleCoords(umperpx(fc));
    
    % get radii as the average of the saved width measurements along
    % that edge. The first column is the actual width, the second is
    % where along the edge it is measured
    radii = zeros(NT.nedge,1);
    for ec = 1:NT.nedge        
         radii(ec) = mean(NT.edgewidth{ec}(:,1));
    end

    %allradii{fc} = radii*umperpx(fc);
    allradii{fc} = radii;

    %% reorder edges and nodes arrays to follow directed tree
    isset = false(1,NT.nedge);
    wasreversed = false(1,NT.nedge);    
    rootnodes(fc) = parentnode;
    directedTreeEdges(NT,parentnode,isset,wasreversed);
    NT.rootnode = parentnode;
    allNetworks(fc) = NT;  

    trunk = NT.nodeedges(NT.rootnode,1);    
       usenetwork(fc) = (size(NT.edgewidth{trunk},2)==6);    
end
usenetwork
exptradii = allradii;

%% Assign radii according to Liao and Howard paper
allrat = [];
p = 2.05;
d0 = 0.24; % minimum diameter in microns
beta = 1980*1e-6; % area per tip in square microns

ntips = []; allrad = [];
for fc = 1:length(allNetworks)
    %%
    NT = allNetworks(fc);
    trunkedge = NT.nodeedges(NT.rootnode,1);

    % trunk radius is 4um
    rtrunk=mean(NT.edgewidth{trunkedge}(:,1));
    trunkradii(fc) = rtrunk*umperpx(fc);
    
    % get number of tips supported by each subtree
    [stL,stEta,stD,muvals,stN] = setSubtreeInfo(NT,trunkedge,2,'L/D');

    radii = (beta*stN+d0^p).^(1/p);
    howardradii{fc} = radii;
    
end


%% work with one particular network 
fc =10;
NT = allNetworks(fc);
disp(networkNames{fc})

plotopt = struct('plotnodes',[]);
NT.plotNetwork(plotopt)

NT.plotNetwork()
hold all
plot(NT.nodepos(NT.rootnode,1),NT.nodepos(NT.rootnode,2),'ro')
hold off

alph=2.0;
%alph=3/2;
%alph=3;
rm = 0.24;
doexptradii =false;
dohoward = false;
doscaledradwexptradii=true;

trunkedge = NT.nodeedges(NT.rootnode,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (doexptradii)
        % set experimental radii
        radii = exptradii{fc}';
        rtrunk = radii(trunkedge);
        stL = []; stV = []; stD = [];
        [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii);

        DR2all{fc} = stD.*radii.^2;
        stVall{fc} = stV;
    elseif (dohoward)
        %set LH scaling
        radii = howardradii{fc};
        rtrunk = radii(trunkedge);
        stL = []; stV = []; stD = [];
        [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii);
     elseif (doscaledradwexptradii)
        %%rm = 0; % no minimal radius               
        rm = 0.12; %min rad
        % rescale radii so that the trunk matches the experimental tree trunk
        % radius
        radii = exptradii{fc}';
        rtrunk = radii(trunkedge);
        [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'L/D');
        %[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'equal');
         
        % set radii with power law scaling
        %%radii = setRadiiFromSubtreeInfo(NT,trunkedge,alph,rtrunk,muvals);
        % set radii with rm parameter
        [radii,stV] = setRadiiWithRm(NT,trunkedge,alph,rm,rtrunk,stL./stD);
    
    else 
        %set radial scaling
        %trunkedge = NT.nodeedges(NT.rootnode,1);
        %rtrunk = mean(NT.edgewidth{trunkedge}(:,1))*umperpx(fc);
        rm = 0; % no minimal radius
        %%rm = 0.12;
        %%rtrunk = 1; % trunk radius set to 1
        rtrunk = 1;
        [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'L/D');
        %[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'equal');
         
        % set radii with power law scaling
        radii = setRadiiFromSubtreeInfo(NT,trunkedge,alph,rtrunk,muvals);
        % set radii with rm parameter
        [radii,stV] = setRadiiWithRm(NT,trunkedge,alph,rm,rtrunk,stL./stD);

    end

% set edgevals to radii^2
        % set edgewidths to radii
   NT.edgewidth = [];
   NT.edgevals = [];
   for ec = 1:NT.nedge
        NT.edgewidth(ec) = radii(ec);
        NT.edgevals(ec) = radii(ec).^2;
   end
Vtree=pi*sum(NT.edgelens.*radii'.^2)        
MitoUnits=Vtree*0.2/0.5
max(radii)        

sr=0;
countsr=0;
for ec=1:NT.nedge
    n2=NT.edgenodes(ec,2);
    deg=NT.degrees(n2,1);
    if(deg == 1 && ec~=trunkedge)
       sr= sr+ radii(ec);
       countsr=countsr+1;
    end
end
sr/countsr
   

voldensity = 0.2; % volume density (unitless)
clustsize = 2.0;
proxstoprate = 0.025; % stopping rate in primary trunk in s^-1
v = 0.6;% velocity in um/s
kb=0.01;
mitovol = 0.5;% volume of 1 mito per cubic micron
ngamma=20;
gammalist=linspace(0,3,ngamma);
asyg=zeros(1,length(gammalist));
for i=1:length(gammalist)
    %% old
    %Aufit= fitAuForClusterSize(NT,radii,kp,kb,gammalist(i),v,AC);
    %Au=Aufit;
    %alpj = v*Au./(2*kb*radii.^gammalist(i));
    %rhoWtrunk1 = 1;
    %[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
    %[M0vals, M1vals] = getMitoClusterDensity(rhoWvals*2*kp/v,alpj);
    %%

    %%mito vol and av cluster size  constant
    %[Au, kp, kb] = fitAukpkb(NT,radii,gammalist(i), v, mitovol, voldensity,clustsize,proxstoprate);
    [Au, kp] = fitAukp(NT,radii,gammalist(i), v, mitovol, voldensity,clustsize,kb);
    %% get densities of mito clusters
    % compute motile mito density, first without 2*kp/v prefactor
    rhoWtrunk1 = 1;
    [rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
    rhoWvals = rhoWvals*2*kp/v; % density in real units (per um)
    % alpha parameter on each edge
    alpj = v*Au./(2*kb*radii.^gammalist(i));
    % get densities of each different size stationary cluster
    [M0vals, M1vals,rhoSvals] = getMitoClusterDensity(rhoWvals,alpj);
    %%
    stM = getSubtreeMitoMass(NT,trunkedge,M1vals.*NT.edgelens');
    stL = []; stV = []; stD = [];
    [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii);
    [asymmetry,nontermjunc] = getAsymmetry(NT,stM, stV);
    ASY2N=asymmetry.*asymmetry;
    countN=0;
    for j=1:NT.nnode
       if (NT.degrees(j)==1)
          countN= countN + 1;
       end
   end
   asyg(i) =sqrt(sum(ASY2N)/(NT.nnode-countN));
end    
plot(gammalist,asyg)
ASYGG=[gammalist' asyg'];
%save("../results/MeanFieldPlots/XALP2ASYGT10.dat","ASYGG","-ascii")
%save("../results/MeanFieldPlots/XALP32ASYGT10.dat","ASYGG","-ascii")
%save("../results/MeanFieldPlots/XALP3ASYGT10.dat","ASYGG","-ascii")
%save("../results/MeanFieldPlots/XEQALP2ASYGT10.dat","ASYGG","-ascii")
save("../results/MeanFieldPlots/XBMINASYGT10.dat","ASYGG","-ascii")

%get error bars
clear
hold off
ND=20;
ASY=cell(1,ND);
Nexp=10;
for i=1:Nexp
  %name=sprintf('../results/MeanFieldPlots/XALP2ASYGT%d.dat',i);
  %name=sprintf('../results/MeanFieldPlots/XALP32ASYGT%d.dat',i);
  %name=sprintf('../results/MeanFieldPlots/XALP3ASYGT%d.dat',i);
  %name=sprintf('../results/MeanFieldPlots/XEQALP2ASYGT%d.dat',i);
  name=sprintf('../results/MeanFieldPlots/XBMINASYGT%d.dat',i);
  data=load(name);
  g=data(:,1);
  for j=1:ND
         ASY{1,j}=[ASY{1,j} ; data(j,2)];
  end  
end

stda=zeros(ND,1);
for j=1:ND
      B=ASY{j}(~isnan(ASY{j}));
      C=B(~isinf(B));
      stda(j)=std(C);
      %stda(j)=std(ASY{j});
end


%meand=load('../results/MeanFieldPlots/EXALP2ASYG.dat'); % run XALP2ASYG.c
%meand=load('../results/MeanFieldPlots/EXALP32ASYG.dat');% run XALP32ASYG.c
%meand=load('../results/MeanFieldPlots/EXALP3ASYG.dat');% run XALP3ASYG.c
%meand=load('../results/MeanFieldPlots/EXEQALP2ASYG.dat');% run XEQALP2ASYG.c
meand=load('../results/MeanFieldPlots/EXBMINASYG.dat');% run XBMINASYG.c


M=meand(:,2);
errorbar(g,M,stda,'b')
GME=[g M stda];
%save("../results/MeanFieldPlots/GEXALP2ASYGE.dat","GME","-ascii")
%save("../results/MeanFieldPlots/GEXALP32ASYGE.dat","GME","-ascii")
%save("../results/MeanFieldPlots/GEXALP3ASYGE.dat","GME","-ascii")
%save("../results/MeanFieldPlots/GEXEQALP2ASYGE.dat","GME","-ascii")
save("../results/MeanFieldPlots/GEXBMINASYGE.dat","GME","-ascii")


