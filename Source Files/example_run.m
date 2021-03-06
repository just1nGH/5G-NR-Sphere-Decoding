clear all;
addpath('Functions');


Nt = 3;
Nr = 4;
moduTypes = {'qpsk','16QAM','qpsk'};
%constls = genConstls(moduTypes);

% bits per symbol for each antenna
Ks = [2,4,2];
offsets = [0,2,6];
Ms = [4,16,4];

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
N0 = 0.5;
rxSymbs = rxSymbs + sqrt(N0/2)*(randn(Nr,1)+ 1j* randn(Nr,1));

% ==================================================================
outType = 'hard';
[out,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType);
nErrs = sum(out ~= msg);

disp('Sim Results')
disp('------------------')
fprintf("Tx Antennas:%d,Rx Antennas: %d \n", Nt,Nr);
fprintf('Modulations:')
disp(moduTypes)

fprintf('(SD) Error bits: %d out of %d\n',nErrs,length(msg) );
fprintf(' - (Hard SD) %d of %d leaf nodes are visited.', nVistedNodes,prod(Ms) );

% ==================================================================
outType = 'soft';
[outSoft,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType);
nErrs = sum((outSoft < 0) ~= msg);

%fprintf('(Soft SD)Error bits: %d out of %d\n',nErrs,length(msg) );
fprintf(' - (Soft SD)%d of %d leaf nodes are visited\n', nVistedNodes,prod(Ms) );

% ======================================================================
% ZF
eqSymbs = pinv(H)*rxSymbs;
out(:) = 0.0;
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






