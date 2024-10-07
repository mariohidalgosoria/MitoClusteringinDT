figure
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
max(radii)
Vtree=pi*sum(NT.edgelens.*radii'.^2)        
MitoUnits=Vtree*0.2/0.5

%gamma =2.1;
u=logspace(-2,0,20);
v = 0.6;% velocity in um/s
stDE=zeros(1,length(u));
kb=0.01;
gamma= 1.5;
%u=linspace(0.01,1,20);
MC=800;
stDE=zeros(1,length(u));
eps=1;
meps=0.0001;
for i=1:length(u)
 % Initialize starting number
    kpp = 2.9;
    %kpp = 17;
    MT=1;
     % Loop until number reaches zero
 
    while(abs(MC-MT)>eps)
        %fprintf("dif=%g\n",abs(MC-MT));
        kpp = kpp - meps; % Decrease the number by 0.0001
        %fprintf("kpp=%g\n",kpp);
        if (kpp<0)
            break;
        end    
        UC=((rtrunk^gamma)*u(i)*kb)/kpp;
        % motile mito density, without 2*kp/v prefactor
        rhoWtrunk1 = 1;
        % set motile mito density
        [rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);
        alpj = v*UC./(2*kb*radii.^gamma);
        [M0vals, M1vals] = getMitoClusterDensity(rhoWvals*2*kpp/v,alpj);
        totmass = sum(M1vals.*(NT.edgelens'));
        %[stMFN,pmfmN,rho1edgeN] = getmitoMassMeanFieldRT3(NT,kpp,kb,gamma,UC,v,trunkedge,N);
        MT=totmass;
        %fprintf("MT=%g\n",MT);
        
    end
    
    if(kpp>0)
       %stDE(i) = DEMF(NT,kpp,kb,gamma,UC,v,trunkedge,radii,70);
       N=100;
       stDE(i) = DistalEnrichmentMeanFieldRT3(NT,kpp,kb,gamma,UC,v,trunkedge,70,N);
    else
        %stDE(i) =1e5;
        stDE(i) =0.1;
    end

    fprintf("MT=%g, u=%g, kp=%g, stDE=%g\n",MT,u(i),kpp,stDE(i));
end
figure
loglog(u,stDE)
%GMT21=[u' stDE'];
%save("../results/MeanFieldPlots/XDEUG21T10.dat","GMT21","-ascii")
GMT15=[u' stDE'];
save("../results/MeanFieldPlots/XDEUG15T10.dat","GMT15","-ascii")

%Create Error bars
clear
hold off
ND=20;
ASY=cell(1,ND);
Nexp=10;
for i=1:Nexp
  %name=sprintf('../results/MeanFieldPlots/XDEUG21T%d.dat',i); %run XDEUG21.c
  name=sprintf('../results/MeanFieldPlots/XDEUG15T%d.dat',i); %run XDEUG15.c
  data=load(name);
  g=data(:,1);
  for j=1:ND
         ASY{1,j}=[ASY{1,j} ; data(j,2)];
  end  
end

stda=zeros(ND,1);
for j=1:ND
      stda(j)=std(ASY{j});
end


%meand=load('../results/MeanFieldPlots/EXDEUG21.dat');
meand=load('../results/MeanFieldPlots/EXDEUG15.dat');
M=meand(:,2);
errorbar(g,M,stda,'b')
GME=[g M stda];
%save("../results/MeanFieldPlots/GEXDEUG21E.dat","GME","-ascii")
save("../results/MeanFieldPlots/GEXDEUG15E.dat","GME","-ascii")


MR=[5.3594 3.5606 6.1055 4.9033 4.9277 5.1563 4.5046 3.5167 4.7407 3.3041];
mean(MR)


