
clear all;
addpath('Functions');

Nt = 4;
Nr = 8;
moduTypes = cell(1,Nt);
for i = 1:ceil(Nt/2)
    moduTypes{1,i} = 'qpsk';
end
for i = ceil(Nt/2)+ 1: Nt
    moduTypes{1,i} = '64qam';
end

constls = genConstls(moduTypes); 

Ks = zeros(1, Nt);
offsets = zeros(1, Nt+1);
for i = 1: Nt
    Ks(i) = constls{1,i};
    offsets(i+1) = offsets(i)+Ks(i);
end
offsets(end) = [];
Ms = 2.^Ks;

% number bits to transform
totBits = sum(Ks);

%generate msg bits
msg = randi([0,1],totBits,1);

txSymbs = zeros(Nt,1);
% modulation
for i = 1: Nt
    txSymbs(i) =  nrModuMapper(msg(offsets(i) + (1: Ks(i))),moduTypes{i});
end


%flat channel
H = sqrt(1/2)*(randn(Nr,Nt) + 1j* randn(Nr,Nt));
H = H/sqrt(Nt);
rxSymbs = H*txSymbs;

% noise
N0 = 0.01;
rxSymbs = rxSymbs + sqrt(N0/2)*(randn(Nr,1)+ 1j* randn(Nr,1));

%------------------------------------------------------------------------- 
outType = 'hard';
[out,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType);

nErrs = sum(out ~= msg);

disp('Sim Results')
disp('------------------')
fprintf("Tx Antennas:%d,Rx Antennas: %d \n", Nt,Nr);
fprintf('Modulations:')
disp(moduTypes)
fprintf('(SD) Error bits: %d out of %d\n',nErrs,length(msg) );
fprintf('--(Hard) %d of %d leaf nodes are visited.', nVistedNodes,prod(Ms) );
%------------------------------------------------------------------------- 
outType = 'soft';
[outSoft,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType);

nErrs = sum( (outSoft < 0) ~= msg);
%fprintf('(Soft) Error bits: %d out of %d\n',nErrs,length(msg) );
fprintf(' (Soft) %d of %d leaf nodes are visited\n', nVistedNodes,prod(Ms) );

% ======================================================================
% ZF
eqSymbs = pinv(H)*rxSymbs;
for i = 1 : Nt
 out(offsets(i) + (1: Ks(i))) = nrSoftModuDemapper(eqSymbs(i),moduTypes{i},N0,'max-log-map');
end
nErrs = sum((out<0) ~= msg);
fprintf('(ZF) Error bits: %d out of %d\n',nErrs,length(msg) );

% ======================================================================
% MMSE
eqSymbs = inv(H' * H + N0*eye(Nt))*H'*rxSymbs;
for i = 1 : Nt
 out(offsets(i) + (1: Ks(i))) = nrSoftModuDemapper(eqSymbs(i),moduTypes{i},N0,'max-log-map');
end
nErrs = sum((out<0) ~= msg);
fprintf('(MMSE) Error bits: %d out of %d\n',nErrs,length(msg) );














