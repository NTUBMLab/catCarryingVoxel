
function [meanvalue, voxelvalue, voxmni, voxcor] = catSPVoxel(sp_center, sp_radius, data)
% catSPVoxel ver 0.2
% Create spherical roi and extract voxel value from your data.
%
% Usage: [meanvalue, voxelvalue, voxmni, voxcor] = catSPVoxel(sp_center, sp_radius, data)
%
%   Output:
%       meanvalue: the mean voxel value for each data. (mask x data matirx)
%       voxelvalue: give you every voxel value for each data. (mask x data)
%       voxmni: the MNI coordinates for each voxel in "voxelvalue".
%       voxcor: the index of 3D image matrix for each voxel.
%
%   Input:
%       sp_center: center of sphereical ROI. (n x 3 numeric matrix)
%
%       sp_radius: radius of sphereical ROI. (n x 1 or 1 x 1 matrix)
%
%       data: string or cell array for your data image. (*.nii or *.img)
%
% Example:
%   [meanvalue, voxelvalue, voxmni, voxcor] = catCarryingVoxel(...
%       spm_select(Inf, 'image', 'Select mask file'), spm_select(Inf, ...
%       'image', 'Select data file'), 1, 0.5);
%
% Dependicies: SPM (tested on SPM12 revision 6906)
%
% Version log:
% Yu-Shiang Su Mar 01 2017: ver 0.1 created this script.
% Yu-Shiang Su Mar 01 2017: ver 0.2 minor bug fixed.



%% Get image header
% if the data is cell array, transform to string array. It will be more
% easy to get the index from each image volume.
if iscell(data)
    data_V = spm_vol(char(data));
else
    data_V = spm_vol(data);
end

sp_n = size(sp_center,1);
if size(sp_radius, 1) ~= size(sp_center, 1)
    if size(sp_radius, 1) == 1
        sp_radius = repmat(sp_radius, size(sp_center, 1), 1);
    else
        fprintf('%s\n', 'the length of Sphere center is not equal to Sphere raidus')
        return
    end
end

byeachdata = 0; % default
if any(ismember(cat(1,data_V.dim), data_V(1).dim, 'rows') == 0)
    fprintf('%s\n', 'Your data are not in the same space (V.dim). ROI space will calculate by each data.');
    byeachdata = 1;
end
if any(ismember(reshape(cat(1,data_V.mat)', 16, [])', reshape(data_V(1).mat', 16, [])', 'rows') == 0)
    fprintf('%s\n', 'Your data have different transformation matrix (V.mat), ROI space will calculate by each data.');
    byeachdata = 1;
end
if byeachdata == 1;
    fprintf('%s\n', 'If you have big data here. This will be dramaticaaly slower than regular processing.');
    voxmni = cell(sp_n, length(data_V));
    voxcor = cell(sp_n, length(data_V));
else
    voxmni = cell(sp_n, 1);
    voxcor = cell(sp_n, 1);
end

meanvalue = nan(sp_n, length(data_V));
voxelvalue = cell(sp_n, length(data_V));

%% get coordinate of sphere
for sp_counter = 1:sp_n
    sp_mni = [];
    for sp_x = [0:1:sp_radius(sp_counter), -1:-1:-sp_radius(sp_counter)]
        for sp_y = [0:1:sp_radius(sp_counter), -1:-1:-sp_radius(sp_counter)]
            for sp_z = [0:1:sp_radius(sp_counter), -1:-1:-sp_radius(sp_counter)]
                if sqrt(sp_x^2 + sp_y^2 + sp_z^2) <= sp_radius(sp_counter)
                    sp_mni = [sp_mni; sp_x, sp_y, sp_z];
                end
            end
        end
    end
    sp_mni = [sp_mni + repmat(sp_center(sp_counter, :), size(sp_mni, 1),1) ...
        ones(size(sp_mni,1), 1)];
    
    %% extract voxel value
    if byeachdata == 1
        for data_counter = 1:length(data_V)
            data_voxxyz = sp_mni * inv(data_V(data_counter).mat)';
            rawdata = spm_get_data(data_V(data_counter), data_voxxyz');
            invalidvox = (sum(rawdata, 1) == 0) | isnan(sum(rawdata,1));
            meanvalue(sp_counter, data_counter) = mean(rawdata(:,~invalidvox),2);
            voxelvalue{sp_counter, data_counter} = rawdata(~invalidvox);
            rawdata_mnixyz = (data_V((data_counter)).mat*[data_voxxyz(~invalidvox,1) data_voxxyz(~invalidvox,2) data_voxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
            voxmni{mask_counter, data_counter} = rawdata_mnixyz(:,1:3);
            voxcor{mask_counter, data_counter} = unique(round(data_voxxyz(~invalidvox, 1:3)), 'rows');
            
            fprintf('%s%d%s%d%s%d%s%d%s%d%s\n', 'ROI ', sp_counter, ...
                ', Data ', data_counter, ': Total ', ...
                length(data_voxxyz), ...
                ' voxels (1x1x1mm) defined in ROI, get ', ...
                size(voxcor{sp_counter, 1}, 1), ' valid voxels and ', ...
                size(unique(round(data_voxxyz(invalidvox, 1:3)), 'rows'), 1), ...
                ' invalid voxels from data space.');
        end
    else
        data_voxxyz = sp_mni * inv(data_V(1).mat)';
        rawdata = spm_get_data(data_V, data_voxxyz');
        invalidvox = (sum(rawdata, 1) == 0) | isnan(sum(rawdata,1));
        meanvalue(sp_counter, :) = mean(rawdata(:,~invalidvox),2);
        voxelvalue(sp_counter, :) = mat2cell(rawdata(:,~invalidvox), ones(1,size(rawdata,1)), sum(~invalidvox));
        rawdata_mnixyz = (data_V(1).mat*[data_voxxyz(~invalidvox,1) data_voxxyz(~invalidvox,2) data_voxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
        voxmni{sp_counter, 1} = rawdata_mnixyz(:,1:3);
        voxcor{sp_counter, 1} = unique(round(data_voxxyz(~invalidvox, 1:3)), 'rows');
        
        fprintf('%s%d%s%d%s%d%s%d%s\n', 'ROI ', sp_counter, ...
            ': Total ', ...
            length(data_voxxyz), ...
            ' voxels (1x1x1mm) defined in ROI, get ', ...
            size(voxcor{sp_counter, 1}, 1), ' valid voxels and ', ...
            size(unique(round(data_voxxyz(invalidvox, 1:3)), 'rows'), 1), ...
            ' invalid voxels from data space.');
    end
end