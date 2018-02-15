%TEST_FCNODEMETRICS
%
%   Requires the following modules:
%
%     cpuinfo-m
%     corr-m
%     fconn-m
%
%   Author: Kristian Loewe

%% test 1
fprintf('Test 1 ... \n');

thr = 0.1;

for dtypes = {'single', 'double'}
  dtype = dtypes{1};
  
  for moas = {'pcc', 'tetracc'}
    moa = moas{1};
    func = str2func(moa);
    
    for nT = 20:20:200
      for nV = 50:50:2000
        dat = rand(nT, nV, dtype);
        utm = func(dat);
        mat = toSymMat(utm,0);
        
        a = flat(sum(mat > thr));
        
        if strcmp(moa, 'pcc')
          b = fcNodeMetrics('degree', dat, thr);
          if ~isequal(a,b)
            sod = sum(abs(a - b));
            if sod < nV/100
              fprintf('WARNING: sod: %d   N: %d\n', sum(abs(a - b)), nV);
            else
              error('ERROR: sod: %d   N: %d\n',   sum(abs(a - b)), nV);
            end
          end
          c = fcNodeMetrics('degree', dat, thr);
          assert(isequal(b,c));
        end
        
        c = fcNodeMetrics('degree', dat, thr, ...
          'ConMeasure', moa);
        if strcmp(moa, 'pcc')
          assert(isequal(b,c), '%d\n', sum(abs(b - c)), nV);
        end
        
        for mt = [1 corecnt() proccnt()]
          d = fcNodeMetrics('degree', dat, thr, ...
            'ConMeasure', moa, 'MaxThreads', mt);
          assert(isequal(c,d), '%d\n', sum(abs(c - d)));
        end
        
        for mm = [1 10 100]
          e = fcNodeMetrics('degree', dat, thr, ...
            'ConMeasure', moa, 'MaxMemory', mm);
          assert(isequal(c,e), '%d\n', sum(abs(c - e)));
        end
      end
    end
  end
end
fprintf('[PASSED]\n');
