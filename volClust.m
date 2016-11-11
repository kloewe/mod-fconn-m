function out = volClust(mask,vol,varargin)
%VOLCLUST Volume clustering and cluster properties.
%
%   C = VOLCLUST(B,V) finds clusters and computes cluster properties such
%   as cluster size and mass based on the binary volume B and the data
%   volume V. The results are returned in the structure C with the following
%   fields and subfields:
%
%     Field            Content
%     -----            -------
%     clust            Found clusters.
%
%                      Subfield  Content
%                      --------  -------
%                      n         Number of found clusters.
%                      volume    Label volume.
%                                The value of each voxel is the number of
%                                its containing cluster (values 1 - N,
%                                where N is the number of clusters).
%
%     size             Property: cluster size.
%
%                      Subfield  Content
%                      --------  -------
%                      values    Cluster size values.
%                                Contains for each cluster the number of
%                                comprised voxels. The i-th value corresponds
%                                to the i-th cluster, where 1 <= i <= N.
%                      volume    Cluster size volume.
%                                The value of each voxel is the size of its
%                                containing cluster or NaN if the voxel
%                                doesn't belong to any cluster.
%
%     mass             Property: cluster mass.
%
%                      Subfield  Content
%                      --------  -------
%                      values    Cluster mass values.
%                                Contains for each cluster the sum of all
%                                values within that cluster. The i-th value
%                                corresponds to the i-th cluster, where
%                                1 <= i <= N.
%                      volume    Cluster mass volume.
%                                The value of each voxel is the mass of its
%                                containing cluster or NaN if the voxel
%                                doesn't belong to any cluster.
%
%   C = VOLCLUST(B,V,'PARAM1',VAL1,'PARAM2',VAL2,...) can be used to specify
%   additional parameters and their values:
%
%     Parameter       Value
%     ---------       -----
%     'Connectivity'  Voxel connectivity.
%                     6  ->  6-connected neighborhood
%                     18 -> 18-connected neighborhood  (default)
%                     26 -> 26-connected neighborhood
%
%     'Properties'    Cluster properties.
%                     This parameter specifies which properties are to be
%                     computed using a cell array of property names. By
%                     default all supported properties are computed.
%
%     'Output'        'All' (default) or 'NoVolumes'
%
%   Author: Kristian Loewe

assert(nargin >= 2 && nargin <= 6 && mod(nargin,2) == 0, ...
  'Unexpected number of input arguments.');
assert(isequal(size(mask), size(vol)));

opts.Connectivity = 18;
opts.Output = 'all';
opts.Properties = {'size','mass'};

for i = 1:2:numel(varargin)
  pn = varargin{i};   % parameter name
  assert(ischar(pn), 'Parameter names must be of type char.');
  pv = varargin{i+1}; % parameter value
  switch pn
    case 'Connectivity'
      assert(isscalar(pv) && ismember(pv, [6 18 26]));
    case 'Output'
      assert(ischar(pv) && ismember(lower(pv), {'all', 'novolumes'}));
    case 'Properties'
      if ischar(pv)
        assert(ismember(pv, opts.Properties));
      elseif iscell(pv) && cellfun(@ischar, pv)
        assert(all(ismember(pv, opts.Properties)));
      else
        error('Unexpected parameter value.');
      end
    otherwise
      error('Unexpected parameter name.');
  end
  opts.(pn) = pv;
end

% find connected components in 3D
if verLessThan('matlab', '7.8')                 % variant 1
  [lbl, nC] = bwlabeln(mask, opts.Connectivity);
  out.clust.n = nC;
  out.clust.volume = lbl;
  clear lbl;
  out.clust.volume = cast(out.clust.volume, class(vol));
  
  if ismember('size', opts.Properties)
    out.size.values = zeros(nC, 1, class(vol));
    for iC = 1:nC
      out.size.values(iC) = sum(out.clust.volume(:) == iC);
    end
  end
  
  if ismember('mass', opts.Properties)
    out.mass.values = zeros(nC, 1, class(vol));
    for iC = 1:nC
      out.mass.values(iC) = sum(vol(out.clust.volume(:) == iC));
    end
  end
  
  if ~strcmp(opts.Output, 'all')
    out.clust = rmfield(out.clust, 'volume');
  else
    out.clust.volume(out.clust.volume == 0) = NaN;
  end
  
else                                            % variant 2 (faster)
  cc = bwconncomp(mask, opts.Connectivity);
  out.clust.n = cc.NumObjects;
  
  if strcmp(opts.Output, 'all')
    out.clust.volume = NaN(cc.ImageSize);
    for iC = 1:cc.NumObjects
      out.clust.volume(cc.PixelIdxList{iC}) = iC;
    end
  end
  
  if ismember('size', opts.Properties)
    out.size.values = flat(cellfun(@numel, cc.PixelIdxList));
  end
  
  if ismember('mass', opts.Properties)
    out.mass.values = zeros(cc.NumObjects, 1, class(vol));
    for iC = 1:cc.NumObjects
      out.mass.values(iC) = sum(vol(cc.PixelIdxList{iC}));
    end
  end
end

if ismember('size', opts.Properties) && strcmpi(opts.Output, 'all')
  out.size.volume = NaN(size(mask), class(vol));
  for iC = 1:numel(out.size.values)
    out.size.volume(out.clust.volume == iC) = out.size.values(iC);
  end
end

if ismember('mass',  opts.Properties) && strcmpi(opts.Output, 'all')
  out.mass.volume = NaN(size(mask), class(vol));
  for iC = 1:numel(out.size.values)
    out.mass.volume(out.clust.volume == iC) = out.mass.values(iC);
  end
end

end
