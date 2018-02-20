%TEST_MFCEDGESTATS
%
%   Requires the following modules:
%
%     cpuinfo-m
%     corr-m
%     fconn-m
%     stats-m
%     util-m
%
%   Author: Kristian Loewe

%%
doTest = [1 1 1 1];

dtypes = {'single', 'double'};
fmt = {'%15.9f','%25.17f'};

tVals = [100 300];
nVals = [100 500];

sVals  = [5 10 20 30];
sVals1 = [5 10 20 30];
sVals2 = [5 10 20 30];

R = 2;

mtVals = [1 corecnt() proccnt()];
mmVals = [0 1 10]; % 100];


%% sum, mean, var, std, tstat
if doTest(1)

fprintf('Test 1: sum, mean, var, std, tstat ...\n');

stats = {'sum','mean','var','std','tstat'};
refFuncs = {@sum @mean @var @std @(varargin) refwrap('tstat', varargin{:})};

for iDT = 1:numel(dtypes)
  fprintf('\n--- %s ---\n', dtypes{iDT});

  for iStat = 1:numel(stats)
    prefix = sprintf('mfcEdgeStats(''%s''', stats{iStat});

    fprintf('\n%s, ...)\n', prefix);
    fprintf('----------------------\n');

    for T = tVals           % n. of timepoints
      for N = nVals         % n. of nodes
        for S = sVals       % n. of data sets

          mad = (N*S)*5*eps(dtypes{iDT});
          mad = min(mad, 0.001);

          t = zeros(1,8);   % timings
          str = cell(1,8);  % descriptive strings

          for iR = 1:R      % repetitions

            data = 10 + rand(T, N, S, dtypes{iDT});

            str{1} = sprintf('%s, data)', prefix);
            tic;
            s = mfcEdgeStats(stats{iStat}, data);
            t(1) = t(1) + toc;

            if strcmp(stats{iStat}, 'tstat')
              [s,p] = mfcEdgeStats(stats{iStat}, data);
            end

            for iMT = 1:numel(mtVals)
              mt = mtVals(iMT);
              str{1+iMT} ...
                = sprintf('%s, data, ''MaxThreads'', %d)', prefix, mt);
              tic;
              s2 = mfcEdgeStats(stats{iStat}, data, 'MaxThreads', mt);
              t(1+iMT) = t(1+iMT) + toc;
              assert(isequal(s, s2));
              clear s2;
            end

            for iMM = 1:numel(mmVals)
              mm = mmVals(iMM);
              str{1+numel(mtVals)+iMM} ...
                = sprintf('%s, data, ''MaxMemory'', %d)', prefix, mm);
              tic;
              s2 = mfcEdgeStats(stats{iStat}, data, 'MaxMemory', mm);
              t(1+numel(mtVals)+iMM) = t(1+numel(mtVals)+iMM) + toc;
              d = s2 - s;
              assert(isequal(s, s2) | max(abs(d)) < mad | max(abs(d)./s) < mad);
              clear s2;
            end

            r = NaN(N*(N-1)/2, S, dtypes{iDT});
            for iS = 1:S;
              r(:,iS) = bstats('fr2z', pcc(data(:,:,iS)));
            end
            assert(~any(isnan(r(:))));

            for iE = 1:N*(N-1)/2
              if strcmp(stats{iStat}, 'tstat')
                [sRef,pRef] = refFuncs{iStat}(r(iE,:));
                d = pRef - p(iE);
                if ~isequal(p(iE), pRef) && abs(d) > mad && abs(d)/pRef > mad
                  fprintf([str{1} '  N: %6d  R: %3d' ...
                    '  pRef: ' fmt{iDT} '  p: ' fmt{iDT} ...
                    '  d: %15.9f  \n'], ...
                    N, iR, pRef, p(iE), d);
                end
              else
                sRef = refFuncs{iStat}(r(iE,:));
              end
              d = sRef - s(iE);
              if ~isequal(s(iE), sRef) && abs(d) > mad && abs(d)/sRef > mad
                fprintf([str{1} '  N: %6d  R: %3d' ...
                  '  sRef: ' fmt{iDT} '  s: ' fmt{iDT} ...
                  '  d: %15.9f  \n'], ...
                  N, iR, sRef, s(iE), d);
              end
            end
          end

          fprintf('T = %d; N = %d; S = %d\n', T, N, S);
          for i = 1:numel(str)
            fprintf('%-55s    %f\n', str{i}, t(i)/R);
          end
        end
      end
    end
  end
end

clear data r s;

fprintf('\n[PASSED]\n\n');

end


%% mdiff, tstat2, pairedt
if doTest(2)

fprintf('Test 2: mdiff, tstat2, pairedt, ...\n');

stats = {'mdiff','tstat2','pairedt'};
refFuncs = {...
  @(a,b) refwrap('mdiff', a,b) ...
  @(varargin) refwrap('tstat2', varargin{:}) ...
  @(varargin) refwrap('pairedt', varargin{:})};

for iDT = 1:numel(dtypes)
  fprintf('\n--- %s ---\n', dtypes{iDT});

  for iStat = 1:numel(stats)
    prefix = sprintf('mfcEdgeStats(''%s''', stats{iStat});

    fprintf('\n%s, ...)\n', prefix);
    fprintf('----------------------\n');

    for T = tVals             % n. of timepoints
      for N = nVals           % n. of nodes
        for S1 = sVals1       % n. of data sets
          for S2 = sVals2     % n. of data sets

            if strcmp(stats{iStat}, 'pairedt') && S1 ~= S2
              continue;
            end

            mad = (N*S1*S2)*5*eps(dtypes{iDT});
            mad = min(mad, 0.001);

            t = zeros(1,8);   % timings
            str = cell(1,8);  % descriptive strings

            for iR = 1:R      % repetitions

              data1 = 10 + rand(T, N, S1, dtypes{iDT});
              data2 = 10 + rand(T, N, S2, dtypes{iDT});

              str{1} = sprintf('%s, data1, data2)', prefix);
              tic;
              s = mfcEdgeStats(stats{iStat}, data1, data2);
              t(1) = t(1) + toc;

              if ismember(stats{iStat}, {'tstat2','pairedt'})
                [s,p] = mfcEdgeStats(stats{iStat}, data1, data2);
              end

              for iMT = 1:numel(mtVals)
                mt = mtVals(iMT);
                str{1+iMT} ...
                  = sprintf('%s, data1, data2, ''MaxThreads'', %d)', prefix, mt);
                tic;
                s2 = mfcEdgeStats(stats{iStat}, data1, data2, 'MaxThreads', mt);
                t(1+iMT) = t(1+iMT) + toc;
                assert(isequal(s, s2));
                clear s2;
              end

              for iMM = 1:numel(mmVals)
                mm = mmVals(iMM);
                str{1+numel(mtVals)+iMM} ...
                  = sprintf('%s, data1, data2, ''MaxMemory'', %d)', prefix, mm);
                tic;
                s2 = mfcEdgeStats(stats{iStat}, data1, data2, 'MaxMemory', mm);
                t(1+numel(mtVals)+iMM) = t(1+numel(mtVals)+iMM) + toc;
                d = s2 - s;
                assert(isequal(s, s2) | max(abs(d)) < mad | max(abs(d)./s) < mad);
                clear s2;
              end

              r1 = NaN(N*(N-1)/2, S1, dtypes{iDT});
              for iS1 = 1:S1;
                r1(:,iS1) = bstats('fr2z', pcc(data1(:,:,iS1)));
              end
              assert(~any(isnan(r1(:))));
              r2 = NaN(N*(N-1)/2, S2, dtypes{iDT});
              for iS2 = 1:S2;
                r2(:,iS2) = bstats('fr2z', pcc(data2(:,:,iS2)));
              end
              assert(~any(isnan(r2(:))));

              for iE = 1:N*(N-1)/2
                if ismember(stats{iStat}, {'tstat2','pairedt'})
                  [sRef,pRef] = refFuncs{iStat}(r1(iE,:), r2(iE,:));
                  d = pRef - p(iE);
                  if ~isequal(p(iE), pRef) && abs(d) > mad && abs(d)/pRef > mad
                    fprintf([str{1} '  N: %6d  R: %3d' ...
                      '  pRef: ' fmt{iDT} '  p: ' fmt{iDT} ...
                      '  d: %15.9f  \n'], ...
                      N, iR, pRef, p(iE), d);
                  end
                else
                  sRef = refFuncs{iStat}(r1(iE,:), r2(iE,:));
                end
                d = sRef - s(iE);
                if ~isequal(s(iE), sRef) && abs(d) > mad && abs(d)/sRef > mad
                  fprintf([str{1} '  N: %6d  R: %3d' ...
                    '  sRef: ' fmt{iDT} '  s: ' fmt{iDT} ...
                    '  d: %15.9f  \n'], ...
                    N, iR, sRef, s(iE), d);
                end
              end
            end

            fprintf('T = %d; N = %d; S1 = %d; S2 = %d\n', T, N, S1, S2);
            for i = 1:numel(str)
              fprintf('%-55s    %f\n', str{i}, t(i)/R);
            end
          end
        end
      end
    end
  end
end

clear data1 data2 s r1 r2;

fprintf('[PASSED]\n\n');

end


%% didt

if doTest(3)

fprintf('Test 3: didt \n');

stats = {'didt'};

for iDT = 1:numel(dtypes)
  fprintf('\n--- %s ---\n', dtypes{iDT});

  for iStat = 1:numel(stats)
    prefix = sprintf('mfcEdgeStats(''%s''', stats{iStat});

    fprintf('\n%s, ...)\n', prefix);
    fprintf('----------------------\n');

    for T = tVals             % n. of timepoints
      for N = nVals           % n. of nodes
        for S1 = sVals1       % n. of data sets group 1
          for S2 = sVals2     % n. of data sets group 2

            mad = (N*S1*S2)*5*eps(dtypes{iDT});
            mad = min(mad, 0.001);

            t = zeros(1,8);   % timings
            str = cell(1,8);  % descriptive strings

            for iR = 1:R      % repetitions

              dataA1 = 10 + rand(T, N, S1, dtypes{iDT});
              dataA2 = 10 + rand(T, N, S1, dtypes{iDT});
              dataB1 = 10 + rand(T, N, S2, dtypes{iDT});
              dataB2 = 10 + rand(T, N, S2, dtypes{iDT});

              str{1} = sprintf('%s, data)', prefix);
              tic;
              s = mfcEdgeStats(stats{iStat}, dataA1, dataA2, dataB1, dataB2);
              t(1) = t(1) + toc;
            end
          end
        end
      end
    end
  end
end

end


%% corrv
if doTest(4)

fprintf('Test 4: corrv ...\n');

stats = {'corrv'};

refFuncs = {@(varargin) refwrap('corrv', varargin{:})};

for iDT = 1:numel(dtypes)
  fprintf('\n--- %s ---\n', dtypes{iDT});

  for iStat = 1:numel(stats)
    prefix = sprintf('mfcEdgeStats(''%s''', stats{iStat});

    fprintf('\n%s, ...)\n', prefix);
    fprintf('----------------------\n');

    for T = tVals           % n. of timepoints
      for N = tVals         % n. of nodes
        for S = tVals       % n. of data sets

          mad = (N*S)*5*eps(dtypes{iDT});
          mad = min(mad, 0.001);

          t = zeros(1,8);   % timings
          str = cell(1,8);  % descriptive strings

          for iR = 1:R      % repetitions

            data = 10 + rand(T, N, S, dtypes{iDT});
            v = rand(S, 1, dtypes{iDT});

            str{1} = sprintf('%s, data, v)', prefix);
            tic;
            s = mfcEdgeStats(stats{iStat}, data, v);
            t(1) = t(1) + toc;

            for iMT = 1:numel(mtVals)
              mt = mtVals(iMT);
              str{1+iMT} ...
                = sprintf('%s, data, v, ''MaxThreads'', %d)', prefix, mt);
              tic;
              s2 = mfcEdgeStats(stats{iStat}, data, v, 'MaxThreads', mt);
              t(1+iMT) = t(1+iMT) + toc;
              assert(isequal(s, s2));
              clear s2;
            end

            for iMM = 1:numel(mmVals)
              mm = mmVals(iMM);
              str{1+numel(mtVals)+iMM} ...
                = sprintf('%s, data, v, ''MaxMemory'', %d)', prefix, mm);
              tic;
              s2 = mfcEdgeStats(stats{iStat}, data, v, 'MaxMemory', mm);
              t(1+numel(mtVals)+iMM) = t(1+numel(mtVals)+iMM) + toc;
              assert(isequal(s, s2) | max(abs(s2 - s)) < mad);
              clear s2;
            end

            r = NaN(N*(N-1)/2, S, dtypes{iDT});
            for iS = 1:S;
              r(:,iS) = bstats('fr2z', pcc(data(:,:,iS)));
            end
            assert(~any(isnan(r(:))));

            for iE = 1:N*(N-1)/2
              sRef = refFuncs{iStat}(r(iE,:), v(:));
              d = sRef - s(iE);
              if ~isequal(s(iE), sRef) && abs(d) > mad && abs(d)/sRef > mad
                fprintf([str{1} '  N: %6d  R: %3d' ...
                  '  sRef: ' fmt{iDT} '  s: ' fmt{iDT} ...
                  '  d: %15.9f  \n'], ...
                  N, iR, sRef, s(iE), d);
              end
            end

            fprintf('T = %d; N = %d; S = %d\n', T, N, S);
            for i = 1:numel(str)
              fprintf('%-55s    %f\n', str{i}, t(i)/R);
            end
          end
        end
      end
    end
  end
end

clear data v r s;

fprintf('[PASSED]\n\n');

end
