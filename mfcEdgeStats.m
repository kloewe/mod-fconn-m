function varargout = mfcEdgeStats(fun,data,varargin)
%MFCEDGESTATS Edge-level statistics across multiple func. connectomes.
%
%   A = MFCEDGESTATS(FUN1,F,...) computes an N-by-N matrix of connection-level
%   statistics across multiple functional connectomes and returns its upper
%   triangular elements in A, where A is a vector of length N*(N-1)/2 and F is
%   a T-by-N-by-S array containing S functional data sets (for example, from
%   multiple sessions or subjects) from which S functional connectomes are
%   derived. T is the number of observations (corresponding to the number of
%   points in time or scans), and N is the number of variables (nodes,
%   defined, for example, as voxels or grayordinates). That is, in each of the
%   S data sets, there are T observations for each of the N variables.
%
%   FUN1 is a string that specifies which statistic is computed. Valid values
%   for FUN1 are the following:
%
%     FUN1            Description
%     ----            -----------
%     'sum'           sum                    A = MFCEDGESTATS('sum',F)
%
%     'mean'          sample mean            A = MFCEDGESTATS('mean',F)
%
%     'var'           sample variance        A = MFCEDGESTATS('var',F)
%
%     'std'           sample standard dev.   A = MFCEDGESTATS('std',F)
%
%     'tstat'         t statistic            A = MFCEDGESTATS('tstat',F)
%                     (one-sample)
%
%                     [A,P] = MFCEDGESTATS('tstat',F) also returns the
%                     corresponding p values in P.
%
%   A = MFCEDGESTATS(FUN2,F,G,...) can be used to compute statistics based on
%   two samples. F is an T-by-N-by-S array containing S functional data sets,
%   where S = S1 + S2, and G is a vector of length S indicating sample
%   membership for the data sets in F.
%
%   A = MFCEDGESTATS(FUN2,F1,F2) can be used to pass the data sets forming
%   the two samples separately. In this case, F1 (F2) is a T-by-N-by-S1
%   (T-by-N-by-S2) array containing S1 (S2) functional data sets.
%
%   FUN2 is a string that specifies which statistic is computed. Valid values
%   for FUN2 are the following:
%
%     FUN2            Description
%     ----            -----------
%     'mdiff'         mean differences       A = MFCEDGESTATS('mdiff',F,G)
%                                            A = MFCEDGESTATS('mdiff',F1,F2)
%
%     'tstat2'        t statistic            A = MFCEDGESTATS('tstat2',F,G)
%                     (two indep. samples,   A = MFCEDGESTATS('tstat2',F1,F2)
%                      equal variances)
%
%                     The t statistics are computed using
%                       t = (M1-M2) / (S12 * sqrt(1/S1 + 1/S2)), where
%                       S12 = sqrt(((S1-1)*VAR1 + (S2-1)*VAR2) / (S1+S2-2)).
%
%                     [A,P] = MFCEDGESTATS('tstat2',...) also returns the
%                     corresponding p values in P.
%
%     'pairedt'       t statistics           A = MFCEDGESTATS('pairedt',F,G)
%                     (two dependent/paired  A = MFCEDGESTATS('pairedt',F1,F2)
%                      samples)
%
%   A = MFCEDGESTATS('corrv',F,V) can be used to compute for each connection
%   the Pearson correlation between the connectivity values from the
%   individual connectomes and the specified variable V. F is an T-by-N-by-S
%   array containing S functional data sets. V is a numeric vector of
%   length S.
%
%   A = MFCEDGESTATS('didt',FA1,FA2,FB1,FB2,...) computes a difference-in-
%   differences (DiD) t statistic, where FA1 and FA2 (FB1 and FB2) are the
%   first and second measurements from group A (B).
%
%   By default, the functional connectomes are being derived from the
%   functional data using the Pearson sample correlation coefficient. Note
%   that the correlation values are Fisher-z transformed prior to computing
%   any statistic across the connectomes.
%
%   A = MFCEDGESTATS(...,'PARAM1',VAL1,'PARAM2',VAL2,...) can be used to
%   specify the following additional parameters and their values:
%
%     Parameter       Value
%     ---------       -----
%     'MaxThreads'    Max. number of threads used for multi-threading.
%                     'auto' (default)     auto-determine
%                     0                    single-threaded
%                     1-n                  multi-threaded with P threads
%
%     'MaxMemory'     Max. amount of memory that should be used.
%                     'auto' (default)     auto-determine
%                     0                    use as little memory as possible
%                     1-n                  try to use at most n GiB of memory
%
%     'ConMeasure'    Functional connectivity measure.
%                     'Pearson' (default)  Pearson sample corr. coefficient
%                     'tetrachoric'        tetrachoric corr. coefficient
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
assert(ischar(fun), 'The function name must be of type ''char''');

% get function name and optional parameter-value pairs
assert(~isempty(fun), 'Unexpected function name.');
[tf,fno] = ismember(fun, {...
  'sum', ...     //  1
  'mean', ...    //  2
  'var', ...     //  3
  'std', ...     //  4
  'tstat', ...   //  5
  'mdiff', ...   //  6
  'tstat2', ...  //  7
  'pairedt', ... //  8
  'didt', ...    //  9
  'corrv', ...   // 10
  });
assert(tf && fno > 0, 'Unexpected function name: ''%s''.', fun);
fno = int32(fno);

% input data
assert(isnumeric(data) && isreal(data) && ~issparse(data) ...
  && ismember(class(data), {'single', 'double'}), ...
  'Unsupported input data type.');
assert(ndims(data) == 3 && all(size(data) > 1), ...
  'Unexpected input data dimensions.');

% precision
if isa(data, 'single')
  mxFcm = @mxFcmFlt;
else
  mxFcm = @mxFcmDbl;
end

switch fun
  case {'sum', 'mean', 'var', 'std', 'tstat'}
    optargin = varargin;
    n = size(data,3);

    args = {fno, data, int32(n)};

  case {'mdiff', 'tstat2', 'pairedt'}
    assert(isnumeric(varargin{1}) && isreal(varargin{1}) ...
      && ~issparse(varargin{1}), 'Unsupported input data type (arg #3).');

    if isvector(varargin{1})
      g = varargin{1};
      assert(isnumeric(g) && isreal(g) && isvector(g) && ~issparse(g), ...
        'G must be a numeric vector.');
      assert(isequal(class(data), class(g)), ...
        'F and G must have the same data type.');
      data1 = data(g == 0);
      data2 = data(g == 1);
    else
      data1 = data;
      data2 = varargin{1};
    end

    optargin = varargin(2:end);
    n1 = size(data1,3);
    n2 = size(data2,3);

    assert(isequal(class(data1), class(data2)), ...
      'F1 and F2 must have the same data type.');
    assert(size(data1,1) == size(data2,1), ...
      'F1 and F2 must have the same #observations per variable.');
    assert(size(data1,2) == size(data2,2), ...
      'F1 and F2 must have the same #variables per data set.');
    assert(~strcmp(fun, 'pairedt') || n1 == n2, ...
      'If FUN2 is ''pairedt'', S1 must be equal to S2.');

    args = {fno, cat(3,data1,data2), int32([n1 n2])};

  case 'corrv'
    v = varargin{1};
    assert(isnumeric(v) && isreal(v) && isvector(v) && ~issparse(v), ...
      'V must be a numeric vector.');
    assert(isequal(class(data), class(varargin{1})), ...
      'F and G must have the same data type.');

    optargin = varargin(2:end);
    n = size(data,3);

    assert(n == numel(varargin{1}), 'V must have S elements.');

    args = {fno, data, int32(n), varargin{1}};
    
  case 'didt'
    na1 = size(data,3);
    na2 = size(varargin{1},3);
    nb1 = size(varargin{2},3);
    nb2 = size(varargin{3},3);
    assert(na1 == na2 && nb1 == nb2);
    optargin = varargin(4:end);
    args = {fno, cat(3, data, varargin{1:3}), int32([na1 nb1])}; 

  otherwise
    error('Unexpected function name: ''%s.''', fun);
end

opts = checkCommonParams(size(data,2), optargin{:});
opts.MaxMemory = cast(opts.MaxMemory, class(data));
args(end+1:end+4) = {...
  opts.MaxThreads, opts.CacheParam, opts.MaxMemory, opts.ConMeasure};

a = mxFcm(args{:});

if nargout == 2
  if     strcmp(fun, 'tstat')
    p = 2 * cast(1 - tcdf(double(abs(a)), n-1), class(data));
  elseif strcmp(fun, 'tstat2')
    p = 2 * cast(1 - tcdf(double(abs(a)), n1+n2-2), class(data1));
  elseif strcmp(fun, 'pairedt')
    p = 2 * cast(1 - tcdf(double(abs(a)), n1-1), class(data));
  end
  varargout = {a,p};

else
  varargout = {a};
end

end
