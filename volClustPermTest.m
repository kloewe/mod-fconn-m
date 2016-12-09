function c = volClustPermTest(vol,volperm,thr,cstat,alpha,varargin)
%VOLCLUSTPERMTEST FWE-corrected cluster-level p values based on permutations.
%
%   C = VOLCLUSTPERMTEST(V,VP,THR,CSTAT,ALPHA) assesses the statistical
%   significance of the clusters in V with respect to the specified cluster
%   statistic CSTAT using a suprathreshold cluster test based on
%   permutations (see Nichols and Holmes, 2001). V is an X-by-Y-by-Z volume
%   of statistics, VP is an X-by-Y-by-Z-by-NP array containing NP
%   permutation-derived volumes of statistics. Based on the threshold THR,
%   clusters are formed in the volume V and in each of the
%   permutation-derived volumes in VP, and the specified cluster statistic
%   CSTAT is computed for each cluster, where CSTAT can be 'size' or
%   'mass'. An FWE-corrected p-value is computed for each cluster based on
%   the permutation distribution of the maximal statistic. The results are
%   returned in the structure C with the following fields and subfields:
%
%     Field           Content
%     -----           -------
%     n               Number of found clusters.
%
%     label           Cluster labels.
%          .volume    Label volume.
%                     The value of each voxel is the number of its
%                     containing cluster (values 1 - N, where N is the
%                     number of clusters) or NaN if the voxel doesn't
%                     belong to any cluster.
%
%     CSTAT           Cluster statistic.
%         .values     Array of statistic values.
%                     Contains the specified property for each of the found
%                     clusters. The i-th value corresponds to the i-th
%                     cluster, where 1 <= i <= N.
%         .volume     Volume of statistics.
%                     The value of each voxel is the statistic value of its
%                     containing cluster or NaN if the voxel doesn't belong
%                     to any cluster.
%         .maxstats   Maximal statistics.
%                     Statistic value from the maximal clusters in V and in
%                     each of the volumes in VP. The first value
%                     corresponds to the maximal cluster in V, and the
%                     (i+1)-th value corresponds to maximal cluster in the
%                     i-th permutation-derived volume VP(:,:,:,i).
%         .pvalues    Array of p values.
%                     Contains for each of the found clusters the
%                     FWE-corrected p value (corrected based on the
%                     permutation distribution of the maximal statistic).
%         .pvolume    Volume of p values.
%                     The value of each voxel is the p value of its
%                     containing cluster or NaN if the voxel doesn't belong
%                     to any cluster.
%         .hvalues    Array of test results.
%                     Contains for each of the found clusters the test
%                     result. (CSTAT).hvalues(i) = 1 indicates a
%                     rejection of the null hypothesis at the specified
%                     significance level ALPHA for the i-th cluster,
%                     and (CSTAT).hvalues(i) = 0 indicates a failure to
%                     reject the null hypothesis for the i-th cluster.
%         .hvolume    Volume of test results.
%                     The value of each voxel is the test result of its
%                     containing cluster or NaN if the voxel doesn't belong
%                     to any cluster.
%
%   C = VOLCLUSTPERMTEST(V,VP,THR,CSTAT,ALPHA,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   can be used to specify additional parameters and their values:
%
%     Parameter       Value
%     ---------       -----
%     'Connectivity'  Voxel connectivity.
%                     6  ->  6-connected neighborhood
%                     18 -> 18-connected neighborhood  (default)
%                     26 -> 26-connected neighborhood
%
%   Author: Kristian Loewe

assert(nargin >= 5 && mod(numel(varargin),2) == 0, ...
  'Unexpected number of input arguments.');
assert(ndims(vol) == 3);
assert(ndims(volperm) == 4);
assert(isequal(size(vol), size(volperm(:,:,:,1))));
assert(ismember(class(vol), {'single', 'double'}));

nP = size(volperm, 4);                  % number of permutations

% defaults
debug = 0;
opts.Connectivity = 18;

% optional parameter-value pairs
for i = 1:2:numel(varargin)
  pn = varargin{i};   % parameter name
  assert(ischar(pn), 'Parameter names must be of type char.');
  pv = varargin{i+1}; % parameter value
  switch pn
    case 'Connectivity'
      assert(isscalar(pv) && ismember(pv, [6 18 26]));
    otherwise
      error('Unexpected parameter name.');
  end
  opts.(pn) = pv;
end

% nans
nanmask = isnan(vol);
if any(nanmask(:))
  for iP = 1:nP
    assert(isequal(nanmask, isnan(volperm(:,:,:,iP))));
  end
end

% apply cluster-forming threshold -> binary masks
mask = vol >= thr;
maskperm = volperm >= thr;

% find clusters and compute cluster statistics
c = volClust(mask, vol, ...
  'Connectivity', opts.Connectivity, 'Properties', {cstat});
c.(cstat).maxstats = zeros(nP+1, 1);
if c.n == 0
  c.(cstat).maxstats(1) = 0;
else
  c.(cstat).maxstats(1) = max(c.(cstat).values);
end
for i = 1:nP
  if any(flat(maskperm(:,:,:,i)))
    cperm = volClust(maskperm(:,:,:,i), volperm(:,:,:,i), ...
      'Connectivity', opts.Connectivity, 'Properties', {cstat}, ...
      'Output', 'NoVolumes');
    c.(cstat).maxstats(i+1) = max(cperm.(cstat).values);
  end
end
  
% determine cluster-level p values
c.(cstat).pvalues = pvalClustHlpr(c.(cstat).values, c.(cstat).maxstats);
if 0
  c.(cstat).pvolume = pcorClustHlpr(c.(cstat).volume, c.(cstat).maxstats);
else
  c.(cstat).pvolume = NaN(size(c.(cstat).volume));
  for i = 1:c.n
    c.(cstat).pvolume(c.label.volume == i) = c.(cstat).pvalues(i);
  end
end

% h values
c.(cstat).hvalues = cast(c.(cstat).pvalues < alpha, class(vol));
c.(cstat).hvolume = cast(c.(cstat).pvolume < alpha, class(vol));
c.(cstat).hvolume(isnan(c.label.volume)) = NaN;

% some additional assertions for debugging
if debug
  for i = 1:c.n
    vali = c.(cstat).values(i);
    maxstats = c.(cstat).maxstats;
    pvali = c.(cstat).pvalues(i);
    assert(isequal(pvali, unique(c.(cstat).pvolume(c.label.volume == i))));
    assert(isequal(pvali, sum(maxstats >= vali)/numel(maxstats)));
    assert(~any(isnan(c.(cstat).pvolume(c.label.volume == i))));
  end
  if any(nanmask(:))
    assert(all(isnan(c.(cstat).pvolume(nanmask))));
    assert(all(isnan(c.(cstat).pvolume(isnan(c.label.volume)))));
  end
end

end

function p = pvalClustHlpr(stats,maxstats)
% helper function: compute p values for the observed clusters based on the
% maximal statistics

assert(~isempty(maxstats));
N = numel(maxstats);
p = ones(size(stats));
p(isnan(stats)) = NaN;
maxstats = sort(maxstats);
for i = 1:N
  idx = stats > maxstats(i);
  if ~isempty(find(idx,1))
    p(idx) = (N-i)/N;
  else
    break;
  end
end
end
