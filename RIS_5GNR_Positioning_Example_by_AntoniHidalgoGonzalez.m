%% OTDOA simulation using PRS and RIS
% 1. SCENARIO CONFIGURATION
nFrames = 1;                   
fc = 3e9;                              
numgNBs = 5;                   
rng('default');    

% Wavelength
c = physconst("lightspeed");
lambda = c/fc;

% RIS parameters 
ris.Enable = true;
ris.Size = [29 29 1]; %Ny x Nz x Npol

% RIS element size
ris.dx = lambda/2;
ris.dy = lambda/2;
% Amplitude of reflection coefficients of each RIS element
ris.A = 0.8; 

% Positions
ris.Pos=[3715.11,7532.00,12]; % 2500m and 30º in azimut from gNB{1}
gNBPos= getgNBPositions(5);
gNBPos{1} = [2465.11, 5366.94, 25]; % Added myself for the simulation
gNBPos{3} = gNBPos{1}+[5000,5000,0]; % Added for the baseline
gNBPos{4} = gNBPos{1}+[-5000,5000,0]; % Added for the baseline
angleUE=deg2rad(180);
step = 10;
iter=20;

initial = 80;
distance2UE= [sqrt((initial+(iter-1)*step)^2 - 10^2)*cos(angleUE), sqrt((initial+(iter-1)*step)^2 - 10^2)*sin(angleUE), -10]; %
UEPos = ris.Pos + distance2UE;


idxMin = choosegNB(gNBPos, ris.Pos, UEPos);
% idxMin = 1;  If we allways want to use the first gNB

% 2. CARRIER OBJECT CONFIGURATION
cellIds = randperm(1008,numgNBs) - 1;

% Configure carrier properties
carrier = repmat(nrCarrierConfig,1,numgNBs);
for gNBIdx = 1:numgNBs
    carrier(gNBIdx).NCellID = cellIds(gNBIdx);
    
end

RIScarrier=pre6GCarrierConfig("NSizeGrid", carrier(idxMin).NSizeGrid, "SubcarrierSpacing", carrier(idxMin).SubcarrierSpacing);
RIScarrier.NCellID=carrier(idxMin).NCellID;
validateCarriers(carrier);

% 3. PRS CONFIGURATION
% Slot offsets of different PRS signals
prsSlotOffsets = 0:2:(2*numgNBs - 1);
prsIDs = randperm(4096,numgNBs) - 1;

% Configure PRS properties
prs = nrPRSConfig;
prs.PRSResourceSetPeriod = [10 0];
prs.PRSResourceOffset = 0;
prs.PRSResourceRepetition = 1;
prs.PRSResourceTimeGap = 1;
prs.MutingPattern1 = [];
prs.MutingPattern2 = [];
prs.NumRB = 52;
prs.RBOffset = 0;
prs.CombSize = 12;
prs.NumPRSSymbols = 12;
prs.SymbolStart = 0;
prs = repmat(prs,1,numgNBs);
for gNBIdx = 1:numgNBs
    prs(gNBIdx).PRSResourceOffset = prsSlotOffsets(gNBIdx);
    prs(gNBIdx).NPRSID = prsIDs(gNBIdx);
   
end
RISprs=prs(gNBIdx);


% 4. PDSCH CONFIGURATION
pdsch = nrPDSCHConfig;
pdsch.PRBSet = 0:51;
pdsch.SymbolAllocation = [0 14];
pdsch.DMRS.NumCDMGroupsWithoutData = 1;
pdsch = repmat(pdsch,1,numgNBs);
RISpdsch = pdsch(idxMin);
validateNumLayers(pdsch);

% 5. PATH LOSS CONFIGURATION

plCfg = nrPathLossConfig;
plCfg.Scenario = 'Uma';

% 6. RESOURCES GENERATION (PRS y PDSCH)
totSlots = nFrames*carrier(1).SlotsPerFrame;
prsGrid = cell(1,numgNBs);
dataGrid = cell(1,numgNBs);
for slotIdx = 0:totSlots-1
    [carrier(:).NSlot] = deal(slotIdx);
    [prsSym,prsInd] = deal(cell(1,numgNBs));
    for gNBIdx = 1:numgNBs
        % Create an empty resource grid spanning one slot in time domain
        slotGrid = nrResourceGrid(carrier(gNBIdx),1);

        % Generate PRS symbols and indices
        prsSym{gNBIdx} = nrPRS(carrier(gNBIdx),prs(gNBIdx));
        prsInd{gNBIdx} = nrPRSIndices(carrier(gNBIdx),prs(gNBIdx));

        % Map PRS resources to slot grid
        slotGrid(prsInd{gNBIdx}) = prsSym{gNBIdx};
        prsGrid{gNBIdx} = [prsGrid{gNBIdx} slotGrid];
    end
    % Transmit data in slots in which the PRS is not transmitted by any of
    % the gNBs (to control the hearability problem)
    for gNBIdx = 1:numgNBs
        dataSlotGrid = nrResourceGrid(carrier(gNBIdx),1);
        if all(cellfun(@isempty,prsInd))
            % Generate PDSCH indices
            [pdschInd,pdschInfo] = nrPDSCHIndices(carrier(gNBIdx),pdsch(gNBIdx));

            % Generate random data bits for transmission
            data = randi([0 1],pdschInfo.G,1);
            % Generate PDSCH symbols
            pdschSym = nrPDSCH(carrier(gNBIdx),pdsch(gNBIdx),data);

            % Generate demodulation reference signal (DM-RS) indices and symbols
            dmrsInd = nrPDSCHDMRSIndices(carrier(gNBIdx),pdsch(gNBIdx));
            dmrsSym = nrPDSCHDMRS(carrier(gNBIdx),pdsch(gNBIdx));

            % Map PDSCH and its associated DM-RS to slot grid
            dataSlotGrid(pdschInd) = pdschSym;
            dataSlotGrid(dmrsInd) = dmrsSym;
        end
        dataGrid{gNBIdx} = [dataGrid{gNBIdx} dataSlotGrid];
    end
end
RISprsGrid = prsGrid{idxMin};


% 7. OFDM MODULATION OF PRS AND DATA SIGNAL AT EACH GNB


txWaveform = cell(1,numgNBs);
for waveIdx = 1:numgNBs
    carrier(waveIdx).NSlot = 0;
    txWaveform{waveIdx} = nrOFDMModulate(carrier(waveIdx),prsGrid{waveIdx} ...
        + dataGrid{waveIdx});
end
% Compute OFDM information using first carrier, assuming all carriers are
% at same sampling rate
ofdmInfo = nrOFDMInfo(carrier(1));

%RIS channel configuration
risCh = hpre6GRISChannel("SampleRate",ofdmInfo.SampleRate,"RISSize", ...
    ris.Size,"CarrierFrequency",fc);

%Angle calculation
[azA, elA, azD, elD, azRG, elRG, distRG] = computeAoAandAoD(UEPos, ris.Pos, ...
    gNBPos{idxMin});
fprintf("geometric AoD → Az: %.2f°, El: %.2f°\n", rad2deg(azD), ...
    rad2deg(elD));
% DIRECTIVITY ESTIMATION -> GAIN
[th, phi, gain] = calculateRISgainCoeff(ris.Enable, risCh, azD, ...
    elD, ris);

% 8. ADD RIS DELAY AND PATH LOSS
sampleDelay = zeros(1,numgNBs);
radius = cell(1,numgNBs);
for gNBIdx = 1:numgNBs
   radius{gNBIdx} = rangeangle(gNBPos{gNBIdx}',UEPos');
   delay = radius{gNBIdx}/c;                      % Delay in seconds
   sampleDelay(gNBIdx) = round(delay*ofdmInfo.SampleRate);   % Delay in samples
end
rxWaveform = zeros(length(txWaveform{1}) + max(sampleDelay),1);
rx = cell(1,numgNBs);
for gNBIdx = 1:numgNBs
    % Calculate path loss for each gNB and UE pair
    losFlag = true; % Assuming the line of sight (LOS) flag as true, as we are only considering the LOS path delays in this example
    PLdB = nrPathLoss(plCfg,fc,losFlag,gNBPos{gNBIdx}(:),UEPos(:));
    if PLdB < 0 || isnan(PLdB) || isinf(PLdB)
        error('nr5g:invalidPL',"Computed path loss (" + num2str(PLdB) + ...
            ") is invalid. Try changing the UE or gNB positions, or path loss configuration.");
    end
    PL(gNBIdx) = 10^(PLdB/10);

    % Add delay, pad, and attenuate
    rx{gNBIdx} = [zeros(sampleDelay(gNBIdx),1); txWaveform{gNBIdx}; ...
                zeros(max(sampleDelay)-sampleDelay(gNBIdx),1)]/sqrt(PL(gNBIdx));

    % Sum waveforms from all gNBs
    rxWaveform = rxWaveform + rx{gNBIdx};
end


% 1) Distances gNB-to-RIS and RIS-to-UE
d_g2r     = norm(gNBPos{idxMin} - ris.Pos);  % m
d_r2u     = norm(ris.Pos     - UEPos);      % m

% 2) Delay in samples
delay_g2r = round((d_g2r/c) * ofdmInfo.SampleRate);
delay_r2u = round((d_r2u/c) * ofdmInfo.SampleRate);
totalDelay= delay_g2r + delay_r2u ;

% 3) Path-loss constants
PLlin_g2r = 10^( nrPathLoss(plCfg,fc,true, gNBPos{idxMin}(:), ris.Pos(:)) / 10 );
PLlin_r2u = 10^( nrPathLoss(plCfg,fc,true, ris.Pos(:),      UEPos(:)) / 10 );

% Use the gNB{1} tx signal
txRef   = txWaveform{idxMin};      % [nSamples×1]
nSamp   = numel(txRef);
numGphi    = numel(gain(:,1));
numGth = numel(gain(1,:));
outLen  = totalDelay + nSamp;

% Pre-aloccate: outLen×numGphixnumGth
rxRIS = zeros(outLen,numGphi, numGth);

% For each RIS coefficient:
for k = 1:numGphi
    for l = 1:numGth
    % a) Attenuation gNB→RIS
    RISsignal = txRef / sqrt(PLlin_g2r);

    % b) Aplicar ganancia del RIS (vector “gain”)
    RISsignal = RISsignal * gain(k,l)*ris.A;

    % c) Attenuation RIS→UE
    RISsignal = RISsignal / sqrt(PLlin_r2u);
    RISsignal_nodelay(:,k,l) = RISsignal; %To later calculate SNR
    % d) Padding for the total geometric delay 
    rxRIS(:,k,l) = [ zeros(totalDelay,1);    % ceros first
                     RISsignal ];                  % processed signal
    end
end

% Equal the samples in RIS and gNB
lenGNB = size(rxWaveform,1);      
lenRIS = size(rxRIS(:,1),    1);       
diflen = abs(lenRIS - lenGNB);
% If RISrx is larger, we padd
if lenRIS > lenGNB
    rxWaveform_padded = [ rxWaveform;
                          zeros(diflen, 1) ];
else
   rxWaveform_padded = rxWaveform;  
   rxRIS(end+1:end+diflen, :, :) = 0;
end

% ADD noise
ue=[5e3 0 2];
gnB=[0 0 28];
pathLoss5km_dB  = nrPathLoss(plCfg,fc,true,ue(:),gnB(:));
pathLoss5km_lin = 10.^(-pathLoss5km_dB/10);
SNR_dB  = 10; % Wanted SNR at 5km
SNR_lin = 10.^(SNR_dB/10);
N0      = pathLoss5km_lin/(sqrt(2)*double(ofdmInfo.Nfft)*SNR_lin);
noise = N0*complex(randn(size(rxWaveform_padded)),randn(size(rxWaveform_padded)));
noise_notpadded = noise((numel(noise)-numel(txWaveform{1})):end);

rxWaveform_padded= rxWaveform_padded/10; % Extra Attenuation (10)
rxnoisy = rxWaveform_padded + noise;
rxRISnoisy = rxRIS + rxWaveform_padded+ noise;

% 9.TOA ESTIMATION 

% CALCULATE RIS CORRELATION (1024xNgphixNgth)
delayRIS=zeros(numGphi, numGth);
delaygNB=zeros(numGphi, numGth);


Ncorr = ceil(ofdmInfo.Nfft * RIScarrier.SubcarrierSpacing/15);
for GphiIdx = 1:numGphi
    for GthIdx = 1:numGth
    
     [~,RISmag] = nrTimingEstimate(RIScarrier, ...
         rxRISnoisy(:,GphiIdx,GthIdx),RISprsGrid);
     corrRIS{GphiIdx, GthIdx} = RISmag(1:(Ncorr));
      
     [pks, locs] = findpeaks(corrRIS{GphiIdx, GthIdx}, ...
         'SortStr','descend', 'NPeaks',2);
     % Findpeaks gives the highest two peaks, sort them in time of arrival
     if locs(2) < locs(1)
     temp = locs(1);
     locs(1) = locs(2);
     locs(2) = temp;
     pks = [pks(2),pks(1)];
     end
     
     peak1(GphiIdx,GthIdx)  = pks(1); % Estimated gNB peak
     peak2(GphiIdx,GthIdx)  = pks(2); % Estimated RIS peak
    
     delaygNB(GphiIdx,GthIdx) = locs(1)-1; % Estimated gNB sample delay
     delayRIS(GphiIdx,GthIdx) = locs(2)-1; % Estimated RIS sample delay
     
    end
end
% Add all the correlation, the correct phi and th should be NyxNz
sumCorr = cellfun(@sum, corrRIS);
[bestValue, linIdx]      = max(sumCorr(:));
[bestPhiIdx, bestThIdx]  = ind2sub(size(sumCorr), linIdx);


bestRISPeakValue = peak2(bestPhiIdx, bestThIdx);

% Range difference estimation:
difDist =  c*(delayRIS(bestPhiIdx, bestThIdx) - delaygNB(bestPhiIdx, bestThIdx))/ofdmInfo.SampleRate;

% Angle estimation:
estAz = (phi(1,bestPhiIdx));
estEl = th(bestThIdx);

% Graph the correlation result:
Fs   = ofdmInfo.SampleRate;
tMag = (0:numel(corrRIS{1})-1)/Fs;

figure;
plot(tMag, corrRIS{bestPhiIdx, bestThIdx}, 'LineWidth',1.4);
xlabel('Time (s)');
ylabel('Correlation magnitude');
title('Correlation (mag) from nrTimingEstimate');
grid on;

% SNR calculation, it is advisable to do it if the correlations have been
% done separatedly, if not, results may be not correct
S_gNB = mean(abs(txWaveform{idxMin}/sqrt(PL(idxMin))).^2);
S_RIS = mean(abs(RISsignal_nodelay(bestPhiIdx,bestThIdx)).^2);
SNR_lin        = [ S_gNB, S_RIS ] ./ mean(abs(noise_notpadded));
SNRdB     = 10*log10( SNR_lin );    




% 10. UE POSITION ESTIMATION:


% RIS
xyz_est = solvePOS(difDist, distRG, azRG, estAz, gNBPos{idxMin}, estEl, -elRG);

% Error computation

EstErrRIS  =  norm(UEPos - xyz_est);
EstErrRISX = abs( UEPos(1) - xyz_est(1) );     
EstErrRISY = abs( UEPos(2) - xyz_est(2) );     
EstErrRISZ = abs( UEPos(3) - xyz_est(3) );   



%% Functions %%
function [azA, elA, azD, elD, azRG, elRG, distRG] = computeAoAandAoD(UEPos, RISPos, gNBPos)
% computeAoAandAoD  Calculates UE↔RIS and RIS↔gNB angles, always positive
%   [azA, elA, azD, elD, azRG, elRG] = computeAoAandAoD(UEPos, RISPos, gNBPos)
%     azA, elA  : azimuth and elevation of arrival UE→RIS (input angle)
%     azD, elD  : azimuth and elevation of departure RIS→UE (output angle)
%     azRG, elRG: azimuth and elevation of departure RIS→gNB
%     distRG    : distance between RIS and gNB
%   All outputs are in radians in the range [0, 2*pi).

    % ---- UE → RIS
    dx  = RISPos(1) - UEPos(1);
    dy  = RISPos(2) - UEPos(2);
    dz  = RISPos(3) - UEPos(3);
    azA = atan2(dy, dx);
    elA = atan2(dz, hypot(dx, dy));

    % ---- RIS → UE (inverse)
    azD = azA + pi;
    elD = -elA;

    % ---- RIS → gNB
    dx2 = gNBPos(1) - RISPos(1);
    dy2 = gNBPos(2) - RISPos(2);
    dz2 = gNBPos(3) - RISPos(3);
    azRG = atan2(dy2, dx2);
    elRG = atan2(dz2, hypot(dx2, dy2));
    distRG = sqrt(dx2^2 + dy2^2 + dz2^2);

    % --- Ensure positive ranges [0, 2*pi)
    azA  = mod(azA,  2*pi);
    azD  = mod(azD,  2*pi);
    azRG = mod(azRG, 2*pi);

end






function idxMin = choosegNB(gNBPos, RISPos, UEPos)
% choosegNB   Selects the index of the gNB on the same “side” as the UE
%             (relative to the RIS) and, if there are several, the one closest to the RIS.
%
%   gNBPos : 1×N cell array where gNBPos{i} = [x y z] of the i-th gNB
%   RISPos : [x y z] of the RIS
%   UEPos  : [x y z] of the UE
%
%   idxMin : index (1…N) of the selected gNB

    % UE→RIS vector in XY
    ueVec = UEPos(1:2) - RISPos(1:2);

    % Initialize distances to Inf
    N = numel(gNBPos);
    dist2R = inf(N,1);

    % Loop through each gNB
    for i = 1:N
        pos = gNBPos{i}(:);
        relX = pos(1) - RISPos(1);

        % Check if it's on the same “half” (X side):
        if ueVec(1) * relX >= 0
            % Calculate distance to the RIS in the XY plane
            dist2R(i) = norm(pos - RISPos(1:2));
        end
    end

    % Choose the minimum among those that meet the condition
    [minDist, idxMin] = min(dist2R);
    if isinf(minDist)
        error("There is no gNB on the same side as the UE.");
    end
end

function xy_est = solvePOS(drest, dgNBRIS, betaAz, alphaAz, gNBPOS, estTH, psi)

alpha = pi - alphaAz ;
beta = betaAz - pi;


Xg =  gNBPOS(1);
Yg =  gNBPOS(2);
Zg = gNBPOS(3);
num = (dgNBRIS)^2 - (dgNBRIS - drest)^2;
den = 2*( dgNBRIS*(1 - cos(psi)*cos(estTH)*(sin(beta)*sin(alpha) - cos(beta)*cos(alpha)) + sin(psi)*sin(estTH)) - drest);
if abs(den) < eps
    error('very small denominator; revise beta and alpha.');
end
dRU = num / den;


% XYZ 


x = dgNBRIS*cos(beta)*cos(psi) - dRU*cos(alpha)*cos(estTH) + Xg;
y = dgNBRIS*sin(beta)*cos(psi) + dRU*sin(alpha)*cos(estTH) + Yg;
z = dgNBRIS*sin(psi) + dRU*sin(estTH) + Zg;
xy_est = [x, y, z];
end


function [thList, phiList, gain] = calculateRISgainCoeff(enableRIS, risCh, ueAz, ueTh, ris)
% calculateRISgainCoeff   RIS coefficients with controlled beamwidth via tapering and iterative pattern plotting
%   enableRIS: flag to enable RIS optimization
%   risCh    : hpre6GRISChannel object
%   carrier  : carrier config
%   ueAz     : known UE azimuth (rad)

% 1) Constants & RIS geometry
c      = physconst("lightspeed");
lambda = c/risCh.CarrierFrequency;

risSize = risCh.RISSize;  % [M N P]
M = risSize(1); 
N = risSize(2); 
P = risSize(3);
dz = ris.dx/lambda; 
dy = ris.dy/lambda;


anY = ones(1,M);
anZ = ones(1,N);

phiList = deg2rad(180):-deg2rad(5):deg2rad(90);



thList = [
  -0.523598775598303, -0.339836909454126, -0.252680255142078, -0.201357920790330, ...
  -0.167448079219690, -0.143347568905365, -0.125327831168065, -0.111341014340964, ...
  -0.100167421161560, -0.0910347780374153, -0.0834300866106150, -0.0769991406568236, ...
  -0.0714894498855205, -0.0667161484102252, -0.0625407617964915, -0.0588575059470812, ...
  -0.0555841732809176, -0.0526559082615699, -0.0500208568057700, -0.0476370626244031, ...
  -0.0454702124169971, -0.0434919707901158, -0.0416787324225779, -0.0400106743539889, ...
  -0.0384710274073283, -0.0370455098120920, -0.0357218823980788, -0.0344895959616788, ...
  -0.0333395092613021, -0.0322636616682469
];

alphaListY = -2*pi * dy * sin(phiList);
alphaListZ = -2*pi * dz * sin(thList);
K    = numel(phiList);
L    = numel(thList);

if ~enableRIS
    % RIS turned OFF
    gain = zeros(K,L);
    return;
end

    
for k = 1:K
    
    for n = 1:L
       % We extract the gain from the specific phiList(k) and thList(n)
    [gainFA,~,~,~]=array2d(anY,anZ,alphaListY(k),alphaListZ(n),dy,dz,ueAz,ueTh,phiList,thList);

    gain(k,n)      = gainFA;
    end  
end
end

function [gain,FAaz,FAth,Daz]=array2d(anY,anZ,alfaY,alfaZ,dy,dz,angAz,angEl,phi,th)

% Modified function from Josep Parron, I added one dimension extra

%Verificando las entradas
[M,N]=size(anY);
if (M~=1) | (N<1), error('an debe ser un vector fila'); end;
aux=size(alfaY);
if (aux~=1), error('alfa tener dimensiones 1x1');end;
aux=size(dy);
if (aux~=1), error('d tener dimensiones 1x1');end;
if (dy<0), error('d tiene que ser mayor que cero');end;


%Angulo electrico
dpsi = pi/1000;
lenY = 2*pi*(1+2*dy);
lenZ = 2*pi*(1+2*dz);
n_ptosY = lenY/dpsi;
n_ptosZ = lenZ/dpsi;
psi_phi=linspace(-pi-2*pi*dy+alfaY,+pi+2*pi*dy+alfaY,n_ptosY);
psi_th = linspace(-pi-2*pi*dy+alfaZ,+pi+2*pi*dy+alfaZ,n_ptosZ);
%Diagrama de radiacion en el espacio electrico
FA_az=abs(freqz(anY,1,psi_phi));
FA_el=abs(freqz(anY,1,psi_th));
%Angulo real


angfa= linspace(pi,0, 2*360+1);
%psi=k*d*cos(th)+alfa
psi_re_az=2*pi*dy*cos(angfa)+alfaY;
psi_re_th=2*pi*dz*cos(angfa)+alfaZ;
%Diagrama de radiacion en el espacio real
FAaz=abs(freqz(anY,1,psi_re_az));
FAth=abs(freqz(anZ,1,psi_re_th));
%Calculo directividad
raz=find((psi_phi > min(psi_re_az)) & (psi_phi < max(psi_re_az)));
rel=find((psi_th > min(psi_re_th)) & (psi_th < max(psi_re_th)));
[~, phiest] = min(abs(angfa - angAz + pi/2));            % mejor azimut
Daz=4*pi*dy*abs(FAaz(phiest)).^2/dpsi/sum(abs(FA_az(raz).^2));

[~, thest] = min(abs(angfa  - (pi/2- angEl)));              % mejor elevación
Del=4*pi*dz*abs(FAth(thest)).^2/dpsi/sum(abs(FA_el(rel).^2));
gain=Daz*Del;

end
%% From here on if you want to compile the code you have to insert 
% the functions validateCarriers, validateNumLayers, colors and plotGrid if
% you want to plot, the other functions have been omitted or are not
% necessary to compile. Also it is necessary to have the hpre6G functions installed
