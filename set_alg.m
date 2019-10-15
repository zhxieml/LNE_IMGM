% Graph matching module used
% Currently, RRWM (Cho et al's, ECCV10) is adopted
% other algorithms can be also used instead.

nAlg = 1;
if 1
    alg(nAlg).fhandle = @RRWM;
    alg(nAlg).variable = {'affinityMatrix', 'group1', 'group2'};
    alg(nAlg).strName = 'RRWM';
    alg(nAlg).param = {};
    alg(nAlg).postDiscretization = 1;
    alg(nAlg).postProcess = @postGreedy;
    alg(nAlg).color = 'r';
    alg(nAlg).lineStyle = '-';
    alg(nAlg).marker = 'p';
end

% if 1    
%     alg(nAlg).fhandle = @SM;
%     alg(nAlg).variable = {'affinityMatrix'};
%     alg(nAlg).strName = 'SM';
%     alg(nAlg).param = {};
%     alg(nAlg).postDiscretization = 1;
%     alg(nAlg).postProcess = @postGreedy;
%     alg(nAlg).color = 'k';
%     alg(nAlg).lineStyle = '-';
%     alg(nAlg).marker = 'x';
% end

% if 1
%     alg(nAlg).fhandle = @ipfp_gm;
%     alg(nAlg).variable = {'affinityMatrix', 'L12'};
%     alg(nAlg).strName = 'IPFP';
%     alg(nAlg).param = {};
%     alg(nAlg).postDiscretization = 1;
%     alg(nAlg).postProcess = @postGreedy;
%     alg(nAlg).color = 'b';
%     alg(nAlg).lineStyle = '--';
%     alg(nAlg).marker = 'x';
% end