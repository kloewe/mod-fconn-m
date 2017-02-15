function a = mfcmStats(fun,data,varargin)
%MFCMEDGESTATS Element-wise statistics across multiple functional connectomes
%
%   A = MFCMEDGESTATS(FUN,F,...) computes a N-by-N matrix of connection-level
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
%               A = MFCMEDGESTATS('mean',F) computes the matrix of sample
%               means across multiple functional connectivity matrices
%               derived from the input functional data sets in F.
%
%     'std'     Sample standard deviations.
%               A = MFCMEDGESTATS('std',F) computes the matrix of sample
%               standard deviations computed across multiple functional
%               connectivity matrices derived from the input functional
%               data sets in F.
%
%     'var'     Sample variances.
%               A = MFCMEDGESTATS('var',F) computes the matrix of sample
%               variances computed across multiple functional connectivity
%               matrices derived from the input functional data sets in F.
%
%   A = MFCMEDGESTATS(...,'PARAM1',VAL1,'PARAM2',VAL2,...) can be used to
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

% optional parameter-value pairs
opts = checkCommonParams(size(data,2), varargin{:});

% function/statistic
fno = [];
switch fun
  case 'mean'
    if 1
      fno = 1;
    else
      if opts.ConMeasure == 1
        fpCorr = @pcc;     % fpCorr = @(d) pcc(d,'avx',1,48);
      elseif opts.ConMeasure == 2
        fpCorr = @tetracc; % fpCorr = @(d) tetracc(d,'m128i',0,48);
      end
      n = size(data, 2);
      s = size(data, 3);
      a = zeros(n*(n-1)/2, 1, dtype);
      for i = 1:s
        fcmi = atanh(fpCorr(data(:,:,i)));
        a = a + fcmi;
      end
      a = a./s;
    end
  case 'std'
    fno = 2;
  case 'var'
    fno = 3;
  otherwise
    error('Unexpected function name.');
end

if ~isempty(fno)
  args = {data, ...
    opts.MaxThreads, opts.CacheParam, opts.MaxMemory, opts.ConMeasure, ...
    int32(fno)};

  % call mex function
  a = mxFcm(args{:});
end

end
