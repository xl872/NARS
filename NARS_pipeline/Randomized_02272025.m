clc; clear; close all;
addpath(genpath('/home/xy079/Xin_proc/AveKspace/Bruker_Raw_Read_PV360'));

% ===================== USER SETTINGS =====================
path_root = '/home/xy079/data/';
AFNIdata='02272025';
cutoff=3;
dd=dir([AFNIdata,'/23Na_ave_upload/',AFNIdata,'info.mat']);
[ROI,info]=BrikLoad([dd(1).folder,'/ROI019+orig']);
ROI_flip=flip(ROI,1);
ROI_flip=circshift(ROI_flip,4,1);
opt.prefix=[dd(1).folder,'/ROI019bothside'];opt.OverWrite='y';
WriteBrik(ROI+ROI_flip,info,opt)

load([AFNIdata,'/23Na_ave_upload/',AFNIdata,'info.mat']);

dateStr = '/20250227_144637_Mouse_23Naplus1_02272025_1_43';


%% RANDOM results
rng(0);
StimTime=[18];
thrmethod='thr';
thr=0;
BLOCK=2;

load([dd(1).folder,'/MeanKSpaceZIPthr_raw0_18_cutAVE55ALLImage.mat']);
RandomN=1000;
n=size(Kall,5);
k=32;
%%
index_matrix = zeros(k, RandomN);
for i = 1:RandomN
    index_matrix(:, i) = sort(randperm(n, k))';
end

for i=1:RandomN

img_4d=mean(Kall(:,:,:,:,index_matrix(:,i)),5);

opt.prefix=[dd(1).folder,'/random/MeanKSpaceZIP',num2str(i),'_cutAVE',num2str(k)];
prefix=['MeanKSpaceZIP',num2str(i),'_cutAVE',num2str(k)];

opt.OverWrite='y';
opt.BRICK_TYPES = 3; 
info.BRICK_TYPES=3;
WriteBrik(img_4d,info,opt);
delete([opt.prefix,'_minus',num2str(cutoff),'+orig*'])
system(['3dcalc -a ''',opt.prefix,'+orig[',num2str(cutoff),'..',num2str(size(img_4d,4)-1),']'' -expr a -prefix ',opt.prefix,'_minus',num2str(cutoff)])
delete([dd(1).folder,'/random/blur_ave_',prefix,'+orig*'])
system(['3dmerge -1blur_fwhm .8 -doall -prefix ',dd(1).folder,'/random/blur_ave_',prefix,' ',opt.prefix,'_minus',num2str(cutoff),'+orig'])
delete([dd(1).folder,'/random/demean_func_',prefix,'+orig*'])
system(['3dDeconvolve -input ',dd(1).folder,'/random/blur_ave_',prefix,'+orig ' ...
    '    -nfirst 0 -polort 2 -num_stimts 1 -stim_times 1 ' ...
    '    ''1D: ',num2str(StimTime),''' ''BLOCK(',num2str(BLOCK),',1)'' ' ...
    '    -stim_label 1 L_FP -tout -fout -rout ' ...
    '    -bucket ',dd(1).folder,'/random/demean_func_',prefix])

end
%%
p_ROI=[];p_ROI_flip=[];
for i=1:RandomN
    prefix=['MeanKSpaceZIP',num2str(i),'_cutAVE',num2str(k)];
    [StatV,info]=BrikLoad([dd(1).folder,'/random/demean_func_',prefix,'+orig']);
    dof = info.BRICK_STATAUX(10); 
    Pt = 2 * (1 - tcdf(abs(StatV(:,:,:,4)), dof));

    p_ROI(i)=mean(Pt(ROI~=0));
    p_ROI_flip(i)=mean(Pt(ROI_flip~=0));
end
%%
figure(Position=[100 100 500 300])

h1 = histogram(p_ROI, 10,'Normalization', 'probability'); 
h1.FaceColor = [0.2 0.4 0.8];    
h1.EdgeColor = 'none';           
h1.FaceAlpha = 0.5;              
hold on
h2 = histogram(p_ROI_flip, 10,'Normalization', 'probability');
h2.FaceColor = [0.8 0.2 0.2];   
h2.EdgeColor = 'none';
h2.FaceAlpha = 0.5;              
xlabel('p value (T-stat)');
ylabel('Probability Density (%)');

legend({'ROI', 'Control'}, 'Box','off',Location='bestoutside');
fontsize(gcf,16,"points")
grid off;
print(gcf,[AFNIdata,'RANDOM019.png'],'-dpng', '-r300')
exportgraphics(gcf, [AFNIdata,'RANDOM019.pdf'], 'ContentType', 'vector');
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