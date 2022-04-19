% minFunc
fprintf('Compiling minFunc files...\n');
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/mcholC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c

