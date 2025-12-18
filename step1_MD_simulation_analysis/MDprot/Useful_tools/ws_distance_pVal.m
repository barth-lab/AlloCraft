function [wsd,pVal,wsdH0] = ws_distance_pVal(u_samples, v_samples, options)
% ws_distance_pVal Calculates  Wasserstein distance between two discrete
% probability measures as well as a p value based on statistical
% bootstrapping
%
% Usage: 
%   [wsd,pVal] = ws_distance_pVal(u_samples, v_samples)
%   [wsd,pVal] = ws_distance_pVal(u_samples, v_samples, options)
%
% * nBlocks: is the number of bootstrapping iterations 
%
% * combineUV: Specifies mode of calculation, if true, u and v are combined
% and significance is tested against randomly permutated data from both u
% and v nBlocks times. If false, u is used as reference and is divided
% into nBlocks and significance is tested using the bootstrapping method in 
% McClendon et al., J. Chem. Theory Comput. (2012). In this mode, "excess"
% wsd is subtracted and reported to wsdH0. Defaults to false
%
% * p: corresponds to the p-Wasserstein distance, p must be 1 or 2
%
% see also ws_distance

    arguments
        u_samples
        v_samples
        options.nBlocks = 10
        options.combineUV = false
        options.p = 1
    end
    p = options.p;
    nBlocks = options.nBlocks;

    wsd = ws_distance(u_samples, v_samples, p);
    nu = length(u_samples);
    

    if options.combineUV % Combine u and v and randomize the data series
        uv_samples = [u_samples ; v_samples];
        isSig = 0;
        nv = length(v_samples);
        nTot = nu + nv;
        for i = 1:nBlocks
            rpu = randperm(nTot,nu);
            rpv = randperm(nTot,nv);
            wsdTest = ws_distance(uv_samples(rpu), uv_samples(rpv), p);

            if wsdTest >= wsd
                isSig = 1 + isSig;
            end
        end
        pVal = isSig/nBlocks;

    else % Use u_sample as reference for intra-sample variability
        % Make temp variable divided by blocks
        interface = floor(linspace(0, (size(u_samples,1) - mod(size(u_samples,1),nBlocks)), nBlocks+1));
    
        uTemp = zeros((size(u_samples,1) - mod(size(u_samples,1),nBlocks))/nBlocks, nBlocks);
    
        for iblock = 1:nBlocks
            istart = interface(iblock)+1;
            iend = interface(iblock+1);
            uTemp(:,iblock) = u_samples(istart:iend,:);
        end
        
        combinations =  nchoosek(1:nBlocks,nBlocks/2); 
        % Do we repeat 1 2 with 3 4 and then 3 4 with 1 2? No! Since WSDis 
        % is symmetric
        combinations = [combinations flipud(combinations)];
        nCombs = nchoosek(nBlocks,nBlocks/2)/2;
        % Remove first half of combinations to avoid repeats:
        combinations = combinations(1:(nCombs),:);

        wsdH0 = zeros(nCombs,1);
        for iblock = 1:nCombs % Combinations of blocks
        
            wsdH0(iblock) = ws_distance(reshape(uTemp(:,combinations(iblock,1:nBlocks/2)),[],1,1) ...
                , reshape(uTemp(:,combinations(iblock,(nBlocks/2+1):end)),[],1,1), p);
        end
        wsdH0Mean = sum(wsdH0)/nCombs; % Average intra-ensemble variability
        
        % Now get the p-value for every d.o.f., zero out non-significant WS distances,
        % and subtract "excess" distance from significant ones
        pVal = sum((wsd - wsdH0)<0) / nCombs;
        
        wsd = wsd - wsdH0Mean;

    end

end