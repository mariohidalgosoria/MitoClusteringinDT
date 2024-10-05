%% load in all dendritic tree networks
%dirname = '../data/HStrees_MCFO/';
% load code dependencies as needed
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

    allradii{fc} = radii*umperpx(fc);
        
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


%% compare cluster size and total mass (aka: dimensionless volume density) 
% in multiple trees, properly scaled
% assuming radii obey Da Vinci Law

nkp = 1000;
gamma = 2.3;
alph=2;
Aulist = 1e-2;%[1e-3,1e-2,1e-1];
kblist = 1e-2*ones(size(Aulist));
nkp = 1000;

% go through all trees
kphat = zeros(nkp,length(Aulist),length(allNetworks));
u = kphat;
v=1.5;
rm = 0.12;
doexptradii =false;
dohoward = false;


clear clustsize totmass kplist Auhat kbhat
cmat = jet(length(allNetworks));
for fc = 1:length(allNetworks)
    
    NT = allNetworks(fc);   
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
    else 
        %set radial scaling
        %trunkedge = NT.nodeedges(NT.rootnode,1);
        %rtrunk = mean(NT.edgewidth{trunkedge}(:,1))*umperpx(fc);
        rtrunk =3;
        [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'L/D');
        %[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alph,'equal');

        % set radii with power law scaling
        radii = setRadiiFromSubtreeInfo(NT,trunkedge,alph,rtrunk,muvals);
        % set radii with rm parameter
        [radii,stV] = setRadiiWithRm(NT,trunkedge,alph,rm,rtrunk,stL./stD);





    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tree depth
    depth(fc) = stD(trunkedge);
    b=sum(NT.edgelens)/stD(trunkedge);
    %bushiness(fc)=sum(NT.edgelens)/stD(trunkedge);
    %% -----
    % Try a few different values of Au, kb, d, v. Plot vs dimensionless
    % parameters        

    % motile mito density, without 2*kp/v prefactor
    rhoWtrunk1 = 1;
    % set motile mito density
    [rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk1,radii);

    for ac = 1:length(Aulist) % run with different au values
        Au = Aulist(ac);
        kb = kblist(ac);

        % alpha parameter on each edge
        alpj = v*Au./(2*kb*radii.^gamma);
        kpmax = 1./(max(alpj.*rhoWvals)*2/v);
        kplist{ac} = linspace(0.1,kpmax,nkp);
        for kc = 1:nkp
            kp =kplist{ac}(kc);

            % get linear density of clusters and mass on each edge
            [M0vals, M1vals] = getMitoClusterDensity(rhoWvals*2*kp/v,alpj);

            totmass(kc,ac,fc) = sum(M1vals.*(NT.edgelens'));
            totclust(kc,ac,fc) = sum(M0vals.*(NT.edgelens'));
            % average cluster size for tree
            clustsize(kc,ac,fc) = totmass(kc,ac,fc)/sum(M0vals.*(NT.edgelens'));

        end

        % get dimensionless parameters
        Auhat(ac,fc) = Aulist(ac)/rtrunk^gamma;
        kbhat(ac,fc) = kblist(ac)*depth(fc)/v;
        kphat(:,ac,fc) = kplist{ac}*depth(fc)/v;
        u(:,ac,fc) = kphat(:,ac,fc)*Auhat(ac,fc)./kbhat(ac,fc);
    end
end

ulist = logspace(-4,0,1000);
% cluster size
subplot(1,3,1)
for fc = 1:length(allNetworks) 
    semilogy(u(:,:,fc),clustsize(:,:,fc),'.-','Color',cmat(fc,:))
    hold all 
    A1TU=clustsize(:,:,fc);
    AGMT1=[u(:,:,fc) A1TU];
    %name=sprintf('G3ATUT%d.dat',fc);
    %name=sprintf('ALP22G2ATUT%d.dat',fc);
    %name=sprintf('ALP3G2ATUT%d.dat',fc);
    %name=sprintf('ALLOG2ATUT%d.dat',fc);
    name=sprintf('RMG2ATUT%d.dat',fc);
    s1 = '../results/MeanFieldPlots/';
    s2 = name;
    dirdest = strcat(s1,s2);
    %save(dirdest,"AGMT1","-ascii");
end
semilogy(ulist,0.5*(1-ulist + (1+ulist)./(1-ulist)),'k--', 'LineWidth',2)
%plot(ulist,0.5*((2-ulist + ulist.^2)./(1-ulist)),'k--', 'LineWidth',2)
hold off
ylim([1,1e1])
xlim([1e-4,1])
xlabel('u = kphat *Auhat/kbhat')
ylabel('avg cluster size')

% total mito mass
subplot(1,3,2)
for fc = 1:length(allNetworks) 
    %loglog(u(:,:,fc),totmass(:,:,fc).*(Auhat(:,fc)./kbhat(:,fc)),'.-','Color',cmat(fc,:))
    semilogy(u(:,:,fc),totmass(:,:,fc)./kphat(:,:,fc),'.-','Color',cmat(fc,:))
    hold all
    %M1TU=totmass(:,:,fc).*(Auhat(:,fc)./kbhat(:,fc))';
    M1TU=totmass(:,:,fc)./kphat(:,:,fc);
    GMT1=[u(:,:,fc) M1TU];
    %name=sprintf('G3MTUT%d.dat',fc);
    %name=sprintf('ALP22G2MTUT%d.dat',fc);
    %name=sprintf('ALP3G2MTUT%d.dat',fc);
    %name=sprintf('ALLOG2MTUT%d.dat',fc);
    name=sprintf('RMG2MTUT%d.dat',fc);
    s1 = '../results/MeanFieldPlots/';
    s2 = name;
    dirdest = strcat(s1,s2);
    %save(dirdest,"GMT1","-ascii");
end
%plot(ulist,ulist.*(1+(1+ulist)./(1-ulist).^2),'k--', 'LineWidth',2)
semilogy(ulist,(1+(1+ulist)./(1-ulist).^2),'k--', 'LineWidth',2)
hold off
ylabel('total mito mass / kphat')
ylim([1e0,1e5])
xlim([1e-4,1])
xlabel('u = kphat *Auhat/kbhat')

%clsuter size versus total mass
subplot(1,3,3)
for fc = 1:length(allNetworks) 
    loglog(totmass(:,:,fc)./kphat(:,:,fc),clustsize(:,:,fc),'.-','Color',cmat(fc,:))
    hold all 
    AM1TU=clustsize(:,:,fc);
    AMGMT1=[totmass(:,:,fc)./kphat(:,:,fc) AM1TU];
    %name=sprintf('G3AMTUT%d.dat',fc);
    %name=sprintf('ALP22G2AMTUT%d.dat',fc);
    %name=sprintf('ALP3G2AMTUT%d.dat',fc);
    %name=sprintf('ALLOG2AMTUT%d.dat',fc);
    name=sprintf('RMG2AMTUT%d.dat',fc);
    s1 = '../results/MeanFieldPlots/';
    s2 = name;
    dirdest = strcat(s1,s2);
    %save(dirdest,"AMGMT1","-ascii");
end
loglog((1+(1+ulist)./(1-ulist).^2),(1+(1+ulist)./(1-ulist).^2)./(1+(1+ulist)./(1-ulist)),'k--', 'LineWidth',2)
hold off
ylim([1,1e2])
xlim([1,1e5])
xlabel('total mito mass / kphat')
ylabel('avg cluster size')


