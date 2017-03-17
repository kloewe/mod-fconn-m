function k = fcmNodeMetrics(nm,data,thr,varargin)
%FCMNODE Nodal metrics of functional connectomes.
%
%   K = FCMNODEMETRICS(NM,FD,THR) computes the nodal metric specified by NM
%   for the nodes in the functional connectivity graph G implicitly
%   constructed based on the input functional data FD and the binarization
%   threshold THR. FD is a T-by-N array of functional data, where T is the
%   number of observations (corresponding to the number of point in time or
%   scans), and N is the number of variables (nodes, definend, for example,
%   as voxels or grayordinates). That is, there are T observations for each
%   of the N nodes. Valid values for NM are the following:
%
%     Metric    Description
%     ------    -----------
%     'degree'  node degree aka degree centrality
%
%   A = FCMNODEMETRICS(NM,FD,THR,'PARAM1',VAL1,'PARAM2',VAL2,...) can be
%   used to specify the following additional parameters and their values:
%
%     Parameter     Value
%     ---------     -----
%     'MaxThreads'  Max. number of threads used for multi-threading.
%                   'auto' (default)     auto-determine
%                   0                    single-threaded
%                   1-n                  multi-threaded with P threads
%
%     'MaxMemory'   Max. amount of memory that should be used.
%                   'auto' (default)     auto-determine
%                   0                    use as little memory as possible
%                   1-n                  try to use at most n GiB of memory
%
%     'ConMeasure'  Functional connectivity measure.
%                   'Pearson' (default)  Pearson sample corr. coefficient
%                   'tetrachoric'        tetrachoric corr. coefficient
%
%   References
%   ----------
%   If you use this program, please cite:
%
%   Memory-efficient analysis of dense functional connectomes
%   K. Loewe, S.E. Donohue, M.A. Schoenfeld, R. Kruse, C. Borgelt
%   Frontiers in Neuroinformatics 10:50 (2016)
%
%   Fast construction of voxel-level functional connectivity graphs
%   K. Loewe, M. Grueschow, C. Stoppel, R. Kruse, and C. Borgelt
%   BMC Neuroscience 15:78 (2014)
%
%   Author: Kristian Loewe

assert(nargin >= 3 && nargin <= 9, 'Unexpected number of input arguments.');
assert(ischar(nm), 'Unexpected metric name.');

% input data
dtype = class(data);
assert(isnumeric(data) && isreal(data) && ~issparse(data) ...
  && ischar(dtype) && ismember(dtype, {'single', 'double'}), ...
  'Unsupported input data type.');
assert(ndims(data) == 2 && all(size(data) > 1), ...
  'Unexpected input data dimensions.');

% precision
if strcmp(dtype, 'single')
  mxNodeDeg = @mxNodeDegFlt;
else
  mxNodeDeg = @mxNodeDegDbl;
end

% threshold
assert(isfloat(thr) && thr >= -1.0 && thr <= 1.0);

% nodal metric
switch nm
  case 'degree'
    %
  otherwise
    error('Unexpected function name.');
end

% optional parameter-value pairs
opts = checkCommonParams(size(data,2), varargin{:});

args = {data, cast(thr, dtype), ...
  opts.MaxThreads, opts.CacheParam, opts.MaxMemory, opts.ConMeasure};

% call mex function
k = mxNodeDeg(args{:});

k = cast(k, dtype);
end
