%computes the Distal Enrichment/Cluster Size Ratio per Region versus the fusion exponent
clear
hold off
addpath('../../networktools/') %taken from lenafabr
addpath('../../dendriticmito/dendriticTrees_public/') %taken from lenafabr
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


%% Fig4b Distal enrichment and average cluster size ratio
%using networks in mu

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
rm = 0.24;
doexptradii =false;
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
    
    elseif (doscaledradwexptradii)
        rm = 0; % no minimal radius
        % rescale radii so that the trunk matches the experimental tree trunk
        % radius
        radii = exptradii{fc}';
        rtrunk = radii(trunkedge);
        [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'L/D');
        %[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'equal');
         
        % set radii with power law scaling
        radii = setRadiiFromSubtreeInfo(NT,trunkedge,alph,rtrunk,muvals);
        % set radii with rm parameter
        %[radii,stV] = setRadiiWithRm(NT,trunkedge,alph,rm,rtrunk,stL./stD);

       
    else 
        %set radial scaling
        %trunkedge = NT.nodeedges(NT.rootnode,1);
        %rtrunk = mean(NT.edgewidth{trunkedge}(:,1))*umperpx(fc);
        rm = 0; % no minimal radius
        rtrunk = 1; % trunk radius set to 1
        [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'L/D');
        %[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'equal');
         
        % set radii with power law scaling
        radii = setRadiiFromSubtreeInfo(NT,trunkedge,alph,rtrunk,muvals);
        % set radii with rm parameter
        %[radii,stV] = setRadiiWithRm(NT,trunkedge,alph,rm,rtrunk,stL./stD);

        
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
%max(radii)        
%% set of parameters in mu and sec
ngamma=40;
gammalist=linspace(0,3,ngamma);
voldensity = 0.2; % volume density (unitless)
clustsize = 2.0;
proxstoprate = 0.025; % stopping rate in primary trunk in s^-1
v = 0.6;% velocity in um/s
mitovol = 0.5;% volume of 1 mito per cubic micron
kb=0.01;
rcs=zeros(1,length(gammalist));
stDE=zeros(1,length(gammalist));
for i=1:length(gammalist)
    % get parameters
        %fprintf("i=%d; gamma=%g \n",i,gammalist(i));
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

    rcs(i) = getRatioAverageCluster4(NT,70,M0vals,M1vals,trunkedge);
    N=100;
    stDE(i) = DistalEnrichmentMeanFieldRT4(NT,kp,kb,gammalist(i),Au,v,trunkedge,70,N);
 end    
hold off
plot(gammalist,rcs,'g')
hold on
plot(gammalist,stDE,'r')
RCSGG=[gammalist' rcs'];
DEGG=[gammalist' stDE'];



save("../results/MeanFieldPlots/XNRCSGT10.dat","RCSGG","-ascii")
save("../results/MeanFieldPlots/XNDEGT10.dat","DEGG","-ascii")

%Matlab script for obtaining errors

clear
hold off
ND=40;
ASY=cell(1,ND);
Nexp=10;
for i=1:Nexp %files obtained from RCSDEG.m
  name=sprintf('../results/MeanFieldPlots/XNRCSGT%d.dat',i);
  %name=sprintf('../results/MeanFieldPlots/XNDEGT%d.dat',i);
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


%ensemble averages obtained from XNRCSG.c and XNDEG.c
meand=load('../results/MeanFieldPlots/EXNRCSG.dat'); 
%meand=load('../results/MeanFieldPlots/EXNDEG.dat');

M=meand(:,2);
errorbar(g,M,stda,'b')
GME=[g M stda];
save("../results/MeanFieldPlots/GEXNRCSGE.dat","GME","-ascii")
%save("../results/MeanFieldPlots/GEXNDEGE.dat","GME","-ascii")

