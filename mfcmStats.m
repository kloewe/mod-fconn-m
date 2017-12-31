function a = mfcmStats(fun,data,varargin)
%MFCMSTATS Element-wise statistics across multiple functional connectomes.
%
%   A = MFCMSTATS(FUN,F,...) computes a N-by-N matrix of connection-level
%   statistics across multiple functional connectomes and returns its upper
%   triangular elements in A, where A is a vector of length N*(N-1)/2 and F
%   is a T-by-N-by-S array containing S functional data sets (for example,
%   from multiple sessions or subjects) from which S functional connectomes
%   are derived. T is the number of observations (corresponding to the number
%   of points in time or scans), and N is the number of variables (nodes,
%   defined, for example, as voxels or grayordinates). That is, in each of
%   the S data sets, there are T observations for each of the N variables.
%   FUN is a string that specifies which statistic is computed. Valid values
%   for FUN are the following:
%
%     Function  Description
%     --------  -----------
%     'mean'    Sample means.
%               A = MFCMSTATS('mean',F) computes the matrix of sample
%               means across multiple functional connectivity matrices
%               derived from the input functional data sets in F.
%
%     'std'     Sample standard deviations.
%               A = MFCMSTATS('std',F) computes the matrix of sample
%               standard deviations computed across multiple functional
%               connectivity matrices derived from the input functional
%               data sets in F.
%
%     'var'     Sample variances.
%               A = MFCMSTATS('var',F) computes the matrix of sample
%               variances computed across multiple functional connectivity
%               matrices derived from the input functional data sets in F.
%
%     'tstat2'  T statistics (for independent two-sample t-tests).
%               A = MFCMSTATS('tstat2',F,G) computes the matrix of
%               t statistics based on the functional connectivity matrices
%               derived from the input functional data sets in F and the
%               membership vector G. F is an T-by-N-by-S array containing S
%               functional data sets, where S = S1 + S2, and G is a vector
%               of length S indicating sample membership for the data sets
%               in F.
%               Using A = MFCMSTATS('tstat2',F1,F2), the data sets can
%               be passed separately. In this case, F1 (F2) is a T-by-N-by-S1
%               (T-by-N-by-S2)) array containing S1 (S2) functional data
%               sets.
%
%     'corr'    Linear correlation between FC and an additional variable.
%               A = MFCMSTATS('corr',F,V) returns the upper triangular
%               elements of the matrix of correlations between the FC
%               values and the variable V. For each element, the functional
%               connectivity values that are correlated with V come from
%               the corresponding elements of the S functional connectivity
%               matrices derived from the S input functional data sets in F,
%               where F is an T-by-N-by-S array and V is a vector of length S.
%
%   A = MFCMSTATS(...,'PARAM1',VAL1,'PARAM2',VAL2,...) can be used to
%   specify the following additional parameters and their values:
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

assert(nargin >= 2, 'Unexpected number of input arguments.');
assert(ischar(fun), 'Unexpected function name.');

% input data
dtype = class(data);
assert(isnumeric(data) && isreal(data) && ~issparse(data) ...
  && ischar(dtype) && ismember(dtype, {'single', 'double'}), ...
  'Unsupported input data type.');
assert(ndims(data) == 3 && all(size(data) > 1), ...
  'Unexpected input data dimensions.');

% precision
if strcmp(dtype, 'single')
  mxFcm = @mxFcmFlt;
else
  mxFcm = @mxFcmDbl;
end

% get function/statistic and optional parameter-value pairs
[tf,fno] = ismember(fun, {'mean','std','var','tstat2','corr'});
assert(tf && fno > 0, 'Unexpected function name.');
fno = int32(fno);

if fno <= 3  % mean, std, var
  opts = checkCommonParams(size(data,2), varargin{:});
else         % tstat2, corr
  opts = checkCommonParams(size(data,2), varargin{2:end});
end

% collect args
args = {data, ...
  opts.MaxThreads, opts.CacheParam, opts.MaxMemory, opts.ConMeasure, fno};

if fno == 4 ...              % tstat2: mfcmStats('tstat2',F1,F2)
    && size(data,1) == size(varargin{1},1) ...
    && size(data,2) == size(varargin{1},2)
  v = int32([zeros(size(data,3),1); ones(size(varargin{1},3),1)]);
  args{1} = cat(3,data,varargin{1});
  args{end+1} = v;
elseif fno == 4 || fno == 5  % tstat2: mfcmStats('tstat2',F1,F2) or corr
  v = varargin{1};
  assert(isvector(v) && size(data,3) == numel(v));
  args{end+1} = v;
end

% compute statistics (call mxFcm or use an alternative variant)
if 0  % alternative variants
  if opts.ConMeasure == 1
    fpCorr = @pcc;     % fpCorr = @(d) pcc(d,'avx',1,48);
  elseif opts.ConMeasure == 2
    fpCorr = @tetracc; % fpCorr = @(d) tetracc(d,'m128i',0,48);
  end
  switch fun
    case 'mean'
      fno = [];
      n = size(data, 2);
      s = size(data, 3);
      a = zeros(n*(n-1)/2, 1, dtype);
      for i = 1:s
        fcmi = atanh(fpCorr(data(:,:,i)));
        a = a + fcmi;
      end
      a = a./s;
    case 'std'
    case 'var'
    case 'tstat2'
    case 'corr'
  end
end
if ~isempty(fno) % call mxFcm
  a = mxFcm(args{:});
end

end
