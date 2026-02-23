clc; clear; close all;
addpath(genpath('/home/xy079/Xin_proc/AveKspace/Bruker_Raw_Read_PV360'));

% ===================== USER SETTINGS =====================
path_root = '/home/xy079/data/';
AFNIdata='02202025';
cutoff=3;

MEflag=1; %multi echo 1; one echo:0
load([AFNIdata,'/23Na_ave_upload/',AFNIdata,'info.mat']);
dd=dir([AFNIdata,'/23Na_ave_upload/',AFNIdata,'info.mat']);
dateStr   = '/20250220_143453_Mouse_23Naplus1_02202025_1_35';
Exp_Number_all = [27:34, 36:42,44:46,52:54,56:61,63:68,70:76,78:89, 94,95,97:105];
%the NARS-fMRI file, the scan numbers that is not listed indicate the single-pulse scan for basic frequency updates of 23Na, 
% %or EPI and other anatomical scans. No data discarded 

fidall=[];fididx=0;
for Exp_Number = Exp_Number_all
% for Exp_Number=Exp_Number_all
% =========================================================
    path = strcat(path_root, dateStr);   
    
    if ~isfile(fullfile(path,num2str(Exp_Number),'method'))
        continue
    end

    file = fopen(fullfile(path,num2str(Exp_Number),'method'),'r');
    methodTxt = fileread(fullfile(path,num2str(Exp_Number),'method'));
    MD=brukerGetScalarSafe(methodTxt,'Method');
    if strcmp(MD,'User:tnnc_2D_121724_trig')
    %%
    fididx=fididx+1;
    praw = fullfile(path, num2str(Exp_Number), 'rawdata.job0');  
    
    
    file = fopen(praw,'r');
    if file == -1
        arr = -2;
        return
    end

    fid = fread(file,'int32');
    fclose(file);
    fidc = reshape(fid,2,[]);
    fidc = complex(fidc(1,:),fidc(2,:));
    fidall(fididx,:)=fidc;
    
    p    = fullfile(path, num2str(Exp_Number));  
    end
end




acqpTxt   = fileread(fullfile(p,'acqp'));
methodTxt = fileread(fullfile(p,'method'));

% ===================== READ PARAMS (Updated for 3D/EPI) =====================
% PVM_Matrix for 3D is [Read, Phase1, Phase2] -> [20, 15, 10]
matOut      = brukerGetArray(methodTxt,'PVM_Matrix');       
NI          = str2double(brukerGetScalar(acqpTxt,'NI')); % 40 frames
nMovie      = str2double(brukerGetScalar(methodTxt,'PVM_NMovieFrames')); 
nStepsTotal = str2double(brukerGetScalar(methodTxt,'PVM_EncGenTotalSteps')); % 6000

% Steps for Phase Encoding 1 and 2
encSteps1 = brukerGetArray(methodTxt,'PVM_EncGenSteps1'); % Length 6000
encSteps2 = brukerGetArray(methodTxt,'PVM_EncGenSteps2'); % Length 6000

%
off_read  = brukerGetArray(methodTxt,'PVM_EffReadOffset');   % likely ~0
off_p1    = brukerGetArray(methodTxt,'PVM_EffPhase1Offset'); % -2.998...
off_p2    = brukerGetArray(methodTxt,'PVM_EffPhase2Offset'); % 0.499...
off_slice = brukerGetArray(methodTxt,'PVM_EffSliceOffset');  % -0.499...
%
fov = brukerGetArray(methodTxt,'PVM_Fov');

fprintf('Exp %d: 3D Matrix=[%s], Frames=%d, TotalSteps=%d\n', ...
    Exp_Number, num2str(matOut), nMovie, nStepsTotal);


%% AVE results
load([dd(1).folder,'/PhaseSTD',AFNIdata,'.mat'])
index=1:size(fidall,1);
StimTime=[18 138 258 378];
thrmethod='thr';
thr=0;
BLOCK=2;
thr1=quantile(PmaskedSTD(:), 1-0.01*thr);
thr2=quantile(PmaskedSTD(:), 0.01*thr);
index=find(PmaskedSTD>thr2&PmaskedSTD<thr1);
% index=find(PmaskedSTD<thr1);
    

    fidC=mean(fidall(index,:),1);
    
    % ===================== INFER nRead FROM DATA =====================
    % In 3D EPI, the total FID size = nRead * nStepsTotal * (NI or nMovie)
    nRead = numel(fidC) / (nStepsTotal); 
    if abs(nRead - round(nRead)) > 1e-9
        error('FID length mismatch: numel(fidC)=%d not divisible by nStepsTotal=%d', numel(fidC), nStepsTotal);
    end
    nRead = round(nRead);
    fprintf('Inferred nRead (samples per encoding step) = %d\n', nRead);
    
    % ===================== RESHAPE & REORDER =====================
    % Bruker storage for 3D GenSteps is typically [nRead, TotalSteps]
    % Since NI=40 and nMovie=40, each 'frame' is embedded in the encoding steps
    % The method file shows @40*(step), meaning the frame is the fastest inner loop
    k_temp = reshape(fidC, [nRead, nMovie, nStepsTotal/nMovie]); 
    
    % We need to map these to a 4D k-space: [Read, Phase1, Phase2, Frame]
    % From method: PVM_EncGenSteps1 has 15 unique values, Steps2 has 10. 15*10*40 = 6000.
    
    k_4d = reshape(permute(k_temp,[1 3 2]), [nRead, matOut(2), matOut(3),nMovie]); 
    
    % --- offsets: choose effective offsets first ---

    % Decide which axes map to your k_4d dims:
    % You built k_4d as [Read, Phase1, Phase2, Frame], so:
    % sx <- read offset, sy <- phase1 offset, sz <- phase2 (or slice) offset.
    % If your third dimension is truly "phase2" encoding (3D), use off_p2.
    % If it is a slice/partition direction, use off_slice. (Often those match in magnitude.)
    Nx = size(k_4d,1);
    Ny = size(k_4d,2);
    Nz = size(k_4d,3);

    % --- convert to voxel units ---
    % If offsets are mm:
      sx = -(off_read  / (fov(1)/Nx));
      sy = -(off_p1    / (fov(2)/Ny));
      sz = -((off_p2 / (fov(3)/Nz) + off_slice / (fov(3)/Nz)));

    % --- build phase ramp (centered indices) ---
    [x, y, z] = ndgrid( (-floor(Nx/2)):(ceil(Nx/2)-1), ...
                        (-floor(Ny/2)):(ceil(Ny/2)-1), ...
                        (-floor(Nz/2)):(ceil(Nz/2)-1) );
    
    phaseRamp = exp( 1i*2*pi*( sx*x/Nx + sy*y/Ny + sz*z/Nz ) );
    
    % --- apply to every frame ---
    for f = 1:size(k_4d,4)
        k_4d(:,:,:,f) = k_4d(:,:,:,f) .* phaseRamp;
    end
    
    % ZIP
    k_full = zeros(matOut(1), matOut(2), matOut(3),nMovie);
    k_full(1+floor((matOut(1) - nRead)/2):nRead+floor((matOut(1) - nRead)/2), :, :, :) = k_4d; 
    k_4d=k_full;
    
    % ===================== RECONSTRUCTION (3D IFFT) =====================
    fprintf('Performing 3D FFT on 40 movie frames...\n');
    img_4d = zeros(size(k_4d), 'like', k_4d);
    for f = 1:nMovie
        % 3D FFT across Read, Phase1, and Phase2
        img_4d(:,:,:,f) = ifftshift(ifftn(fftshift(k_4d(:,:,:,f))));
    end
    img_4dP=angle(img_4d);
    % img_4dmP=angle(mean(img_4d(:,:,:,cutoff+1:end),4));
    img_4d=abs(img_4d);
    
    % img_4d=((img_4d-min(img_4d(:)))/(max(img_4d(:))-min(img_4d(:))))*255;

    img_4dm=mean(img_4d(:,:,:,cutoff+1:end),4);
    img_4dmP=mean(img_4dP(:,:,:,cutoff+1:end),4);
    
    

%% AVE results AVE2


CutN=floor(length(index)/2);
indexall=index(1:floor(length(index)/CutN)*CutN);
GroupN=floor(length(index)/CutN);
% index=find(PmaskedSTD<thr1);
Kall=[];
    for CN=1:CutN
        index=indexall((CN-1)*GroupN+1:CN*GroupN);
        fidC=mean(fidall(index,:),1);
       
        % ===================== INFER nRead FROM DATA =====================
        % In 3D EPI, the total FID size = nRead * nStepsTotal * (NI or nMovie)
        nRead = numel(fidC) / (nStepsTotal); 
        if abs(nRead - round(nRead)) > 1e-9
            error('FID length mismatch: numel(fidC)=%d not divisible by nStepsTotal=%d', numel(fidC), nStepsTotal);
        end
        nRead = round(nRead);
        fprintf('Inferred nRead (samples per encoding step) = %d\n', nRead);
        
        % ===================== RESHAPE & REORDER =====================
        % Bruker storage for 3D GenSteps is typically [nRead, TotalSteps]
        % Since NI=40 and nMovie=40, each 'frame' is embedded in the encoding steps
        % The method file shows @40*(step), meaning the frame is the fastest inner loop
        k_temp = reshape(fidC, [nRead, nMovie, nStepsTotal/nMovie]); 
        
        % We need to map these to a 4D k-space: [Read, Phase1, Phase2, Frame]
        % From method: PVM_EncGenSteps1 has 15 unique values, Steps2 has 10. 15*10*40 = 6000.
        
        k_4d = reshape(permute(k_temp,[1 3 2]), [nRead, matOut(2), matOut(3),nMovie]); 
        
        % --- offsets: choose effective offsets first ---

    % Decide which axes map to your k_4d dims:
    % You built k_4d as [Read, Phase1, Phase2, Frame], so:
    % sx <- read offset, sy <- phase1 offset, sz <- phase2 (or slice) offset.
    % If your third dimension is truly "phase2" encoding (3D), use off_p2.
    % If it is a slice/partition direction, use off_slice. (Often those match in magnitude.)
    Nx = size(k_4d,1);
    Ny = size(k_4d,2);
    Nz = size(k_4d,3);

    % --- convert to voxel units ---
    % If offsets are mm:
      sx = -(off_read  / (fov(1)/Nx));
      sy = -(off_p1    / (fov(2)/Ny));
      sz = -((off_p2 / (fov(3)/Nz) + off_slice / (fov(3)/Nz)));

    % --- build phase ramp (centered indices) ---
    [x, y, z] = ndgrid( (-floor(Nx/2)):(ceil(Nx/2)-1), ...
                        (-floor(Ny/2)):(ceil(Ny/2)-1), ...
                        (-floor(Nz/2)):(ceil(Nz/2)-1) );
    
    phaseRamp = exp( 1i*2*pi*( sx*x/Nx + sy*y/Ny + sz*z/Nz ) );
    
    % --- apply to every frame ---
    for f = 1:size(k_4d,4)
        k_4d(:,:,:,f) = k_4d(:,:,:,f) .* phaseRamp;
    end

        
        % ZIP
        k_full = zeros(matOut(1), matOut(2), matOut(3),nMovie);
        k_full(1+floor((matOut(1) - nRead)/2):nRead+floor((matOut(1) - nRead)/2), :, :, :) = k_4d; 
        k_4d=k_full(:,:,:,1:end);

        fprintf('Performing 3D FFT on 40 movie frames...\n');
        img_4d = zeros(size(k_4d), 'like', k_4d);
        for f = 1:size(img_4d,4)
            % 3D FFT across Read, Phase1, and Phase2
            img_4d(:,:,:,f) = ifftshift(ifftn(fftshift(k_4d(:,:,:,f))));
        end
        img_4dP=angle(img_4d);
        % img_4dmP=angle(mean(img_4d(:,:,:,cutoff+1:end),4));
        img_4d=abs(img_4d);

        Kall(:,:,:,:,CN)=img_4d;

    end
save([dd(1).folder,'/MeanKSpaceZIP',thrmethod,'_raw',num2str(thr),'_',num2str(StimTime(1)),'_cutAVE',num2str(CutN),'ALLImage.mat'],'Kall');
img_4d=mean(Kall,5);

opt.prefix=[dd(1).folder,'/MeanKSpaceZIP',thrmethod,'_raw',num2str(thr),'_',num2str(StimTime(1)),'_cutAVE',num2str(CutN)];
prefix=['MeanKSpaceZIP',thrmethod,'_raw',num2str(thr),'_',num2str(StimTime(1)),'_cutAVE',num2str(CutN)];

opt.OverWrite='y';
opt.BRICK_TYPES = 3; 
info.BRICK_TYPES=3;
WriteBrik(img_4d,info,opt);
delete([opt.prefix,'_minus',num2str(cutoff),'+orig*'])
system(['3dcalc -a ''',opt.prefix,'+orig[',num2str(cutoff),'..',num2str(size(img_4d,4)-1),']'' -expr a -prefix ',opt.prefix,'_minus',num2str(cutoff)])
delete([dd(1).folder,'/blur_ave_',prefix,'+orig*'])
system(['3dmerge -1blur_fwhm .8 -doall -prefix ',dd(1).folder,'/blur_ave_',prefix,' ',opt.prefix,'_minus',num2str(cutoff),'+orig'])
% delete([dd(1).folder,'/demean_func_',prefix,'+orig*'])
% system(['3dDeconvolve -input ',dd(1).folder,'/blur_ave_',prefix,'+orig ' ...
%     '    -nfirst 0 -polort 2 -num_stimts 1 -stim_times 1 ' ...
%     '    ''1D: ',num2str(StimTime),''' ''BLOCK(',num2str(BLOCK),',1)'' ' ...
%     '    -stim_label 1 L_FP -tout -fout -rout ' ...
%     '    -bucket ',dd(1).folder,'/demean_func_',prefix ...
%     '    -errts ',dd(1).folder,'/errts.',prefix])

if MEflag==1
[Vn,infon]=BrikLoad([dd(1).folder,'/blur_ave_',prefix,'+orig']);

optn.OverWrite='y';
optn.prefix=[dd(1).folder,'/blur_aveALL_',prefix];
Vnn=[];
EN=4;%echo num
for i=1:EN
    Vnn(:,:,:,:,i)=Vn(:,:,:,(i-1)*size(img_4d,4)/EN+1:(i-1)*size(img_4d,4)/EN+size(img_4d,4)/EN-cutoff);
end
Vnn=mean(Vnn,5);
infon.DATASET_RANK(2)=size(Vnn,4);
infon.BRICK_STATS=infon.BRICK_STATS(1:size(Vnn,4)*2);
infon.BBRICK_FLOAT_FACS=infon.BRICK_FLOAT_FACS(1:size(Vnn,4));
infon.TAXIS_NUMS(1)=infon.DATASET_RANK(2);
WriteBrik(Vnn,infon,optn);
delete([dd(1).folder,'/demean_func_',prefix,'AVE+orig*'])
system(['3dDeconvolve -input ',dd(1).folder,'/blur_aveALL_',prefix,'+orig ' ...
    '    -nfirst 0 -polort 2 -num_stimts 1 -stim_times 1 ' ...
    '    ''1D: ',num2str(StimTime),''' ''BLOCK(',num2str(BLOCK),',1)'' ' ...
    '    -stim_label 1 L_FP -tout -fout -rout ' ...
    '    -bucket ',dd(1).folder,'/demean_func_',prefix,'AVE' ...
    '    -errts ',dd(1).folder,'/errts.',prefix,'AVE'])
end

%%
function k = reorder_pe_by_steps(k, encSteps1)
    if isempty(encSteps1), return; end
    nPE = size(k,2);
    if numel(encSteps1) ~= nPE, return; end

    encSteps1 = round(encSteps1(:));
    minS = min(encSteps1);
    maxS = max(encSteps1);
    if (maxS - minS + 1) ~= nPE
        warning('PVM_EncSteps1 not contiguous; skipping PE reorder.');
        return;
    end

    k2 = zeros(size(k), 'like', k);
    for i = 1:nPE
        ky = encSteps1(i) - minS + 1;   % maps [-N/2..N/2-1] -> [1..N]
        k2(:, ky, :) = k(:, i, :);
    end
    k = k2;
end

function k = reorder_obj(k, objOrder0)
    if isempty(objOrder0), return; end
    NI = size(k,3);
    if numel(objOrder0) ~= NI
        warning('ACQ_obj_order length mismatch; skipping obj reorder.');
        return;
    end
    invOrder = zeros(1, NI);
    invOrder(objOrder0 + 1) = 1:NI;
    k = k(:,:,invOrder);
end

function img = ifft2c(k)
    img = ifftshift(ifft(ifftshift(k,1),[],1),1);
    img = ifftshift(ifft(ifftshift(img,2),[],2),2);
end

function y = cropCenter(x, nOut, dim)
    nIn = size(x,dim);
    if nOut > nIn, error('nOut > nIn'); end
    % i0 = floor((nIn - nOut)/2) + 1;
    i0 =  1;
    idx = repmat({':'},1,ndims(x));
    idx{dim} = i0:(i0+nOut-1);
    y = x(idx{:});
end

function v = brukerGetScalar(txt, name)
    expr = ['##\$', name, '=\s*([^\n\r]+)'];
    t = regexp(txt, expr, 'tokens', 'once');
    if isempty(t), error('Param not found: %s', name); end
    v = strtrim(t{1});
    v = regexprep(v, '[<>]', '');
end

function arr = brukerGetArray(txt, name)
    expr = ['##\$', name, '=\(\s*([^\)]*)\)\s*[\r\n]+([^\#]+)'];
    t = regexp(txt, expr, 'tokens', 'once');
    if isempty(t), arr = []; return; end
    raw = strtrim(regexprep(t{2}, '\s+', ' '));
    raw = regexp(raw, '^[^#]*', 'match', 'once'); % stop before next ##
    arr = str2num(raw); %#ok<ST2NM>
end

%%
function v = brukerGetScalarSafe(txt, name, defaultVal)
    expr = ['##\$', name, '=\s*([^\n\r]+)'];
    t = regexp(txt, expr, 'tokens', 'once');
    if isempty(t)
        v = defaultVal;
        return;
    end
    v = strtrim(t{1});
    v = regexprep(v, '[<>]', '');
end