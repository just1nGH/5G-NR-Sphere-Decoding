# 5G-NR-Sphere-Decoding
Matlab Implenmentation of 5G NR MIMO Sphere Decoder
- Tree traverse stratigies: single tree search 
- Support different modulation schemes for each Tx antenna

The implementation sphere decoding algorithms(single tree seach)for seeking the maximum-likelihood solution 
for a set of symbols transmitted over the MIMO channel. 
```
function [out,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType)
% [out,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType) uses sphere decoding algorithms
% (single tree seach)for seeking the maximum-likelihood solution 
%  for a set of symbols transmitted over the MIMO channel. 
% Input:
%   H - channel matrix Nr X Nt
%   rxSymbs - received symbols: Nr X 1
%   moduTypes - nodulation for each tx antenna in a cell array of 1 X Nt
%   outType - 'soft' or 'hard'
% Output:
%   out - soft bits or hard bits output, soft bits (not sacled by 1/N0)
%   nVistedNodes - number of leaf nodes traversed
% source paper: 
%   "Soft-Output Sphere Decoding: Performance and Implementation Aspects"
% Author: Dr J Mao
% Email: juquan.justin.mao@gmail.com
% 2021 Nov
%--------------------------------------------------------------------------
```

### Example run simlation results
Start by try `example_run.m`

```
Sim Results
------------------
Tx Antennas:4,Rx Antennas: 8 
Modulations:    {'qpsk'}    {'qpsk'}    {'64qam'}    {'64qam'}

(Hard SD) Error bits: 0 out of 16
--(Hard) 96 of 65536 leaf nodes are visited. (Soft) 15448 of 65536 leaf nodes are visited
(ZF) Error bits: 9 out of 16
(MMSE) Error bits: 9 out of 16
```
