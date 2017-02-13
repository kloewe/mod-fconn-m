function opts = checkCommonParams(n,varargin)
%CHECKCOMMONPARAMS Check common parameter-value pairs
%
%   Author: Kristian Loewe

assert(mod(numel(varargin),2) == 0, 'Unexpected number of input arguments.');

% defaults
opts.MaxThreads = int32(-1);  % auto
opts.MaxMemory  = -1;         % auto
opts.CacheParam = int32(-1);  % auto
opts.ConMeasure = int32( 1);  % Pearson

for i = 1:2:numel(varargin)
  pn = varargin{i};   % parameter name
  assert(ischar(pn), 'Parameter names must be of type char.');
  pv = varargin{i+1}; % parameter value
  switch pn
    case {'MaxThreads','MaxMemory'}
      if ischar(pv) && strcmp(pv, 'auto')
        pv = -1;
      else
        assert(isscalar(pv) && isnumeric(pv) && isreal(pv) ...
          && ~issparse(pv) && isfinite(pv), ...
          'Unexpected parameter value specified for %s.\n', pn);
      end
      if strcmp(pn, 'MaxThreads')
        pv = int32(min(pv, n));
        assert(pv == -1 || (pv >= 0 && pv <= 1024), ...
          'Unexpected parameter value specified for MaxThreads.');
      elseif strcmp(pn, 'MaxMemory')
        assert(pv <= 2048, ...
          'Unexpected parameter value specified for MaxMemory.');
      end
    case 'ConMeasure'
      assert(ischar(pv), ...
        'ConMeasure needs to be specified as a string.');
      switch pv
        case {'Pearson','pcc'}
          pv = int32(1);
        case {'tetrachoric','tcc', 'tetracc'}
          pv = int32(2);
        otherwise
          error('Unknown functional connectivity measure specified.');
      end
    case 'Paired'
      assert(isscalar(pv) && ~issparse(pv) ...
        && (islogical(v) || v == 0 || v == 1));
    otherwise
      error('Unexpected parameter name %s.',pn);
  end
  opts.(pn) = pv;
end

end
