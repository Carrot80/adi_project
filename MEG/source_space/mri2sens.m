function outcoor=mri2sens(cfg)
%Transforms MRI-coordinates gained from volumenormalized data into Sensor
%coordinates using a provided transformation matrix.
%
%Attention:
%   - When validating your outcoor with outunit cm on MRIs be sure to
%   multiply output by 10 (i.e. convert to mm) as the MRI-units are usually
%   in mm.
%   - If you have -in prior steps- have subtracted the sphere origin of the
%   voxel-points (ususally visible with a z-origin ~0), you need to have to 
%   add it back: outcoor + orig
%
%INPUT:
%   - cfg.source:    volumenormalized sourceanalysis data
%   - cfg.mricoord:  MRI Matrix as 2D (Nx3); rows are locations of interest
%                    e.g. [72 1 67; -48 -29 50]; in mm
%OPTIONAL:
%   - cfg.outunit: 'cm' (default), 'mm'
%
%Nathan August, 2008
 
transMat=cfg.source.cfg.final;
coordMat=cfg.mricoord;
 
if ~isfield(cfg,'outunit')
    d=10; %assumes cm
elseif strcmpi(cfg.outunit,'cm')
    d=10;
elseif strcmpi(cfg.outunit,'mm')
    d=1;
end
 
coordMat(:,4)=ones(1,size(coordMat,1));
 
outcoor=inv(transMat)*coordMat';
 
outcoor=outcoor(1:3,:)'/d;


end

