function constls = genConstls(moduTypes)
% constls = genConstls(moduTypes) output the constellation info of a cell
% list of moduTypes.
% % example
% moduTypes = {'qpsk','qpsk','16QAM'};
% constls = genConstls(moduTypes)

    constls = cell(3,length(moduTypes));
    for i = 1: length(moduTypes)
        moduType = moduTypes{i};
        switch lower(moduType)
        case 'bpsk'
            M =2;
        case 'qpsk'
            M = 4;
        case '16qam'
            M = 16;
        case  '64qam'
            M = 64;
        case '256qam'
            M = 256;
        end
    
        K = log2(M);
        constls{1,i} = K;
    
        bitLabels = de2bi(0:M-1,'left-msb');
        constls{3,i} = bitLabels;
    
        symbBitsIn = bitLabels.';
        constSymbs = nrModuMapper(symbBitsIn(:),lower(moduType));
    
        constls{2,i} = constSymbs;
    
    end
end