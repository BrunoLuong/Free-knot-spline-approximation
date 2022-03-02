function minmaxfilter_install
% function minmaxfilter_install
% Installation by building the C-mex files for minmax filter package
%
% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 20-Sep-2009

here = fileparts(mfilename('fullpath'));
oldpath = cd(here);

arch=computer('arch');
mexopts = {'-O' '-v' ['-' arch]};
if ~verLessThan('MATLAB','9.4')
    R2018a_mexopts = {'-R2018a'};
else
    R2018a_mexopts = {};    
    % 64-bit platform
    if ~isempty(strfind(computer(),'64'))
        mexopts(end+1) = {'-largeArrayDims'};
    end
end

opmopt = {};
% Look for C compiler
cc = mex.getCompilerConfigurations('C');
[ccfound,loc] = ismember('C',{cc.Language});
if ccfound
    cc = cc(loc);
    if contains(cc.ShortName,'MSVC')
        opmopt = {'COMPFLAGS="$COMPFLAGS /openmp"'};
        
    elseif contains(cc.ShortName,'INTEL')
        opmopt = {'COMPFLAGS="$COMPFLAGS /MD /Qopenmp"'};
        
    elseif contains(cc.ShortName,{'gcc' 'mingw64'}) % not tested
        opmopt = {'CFLAGS="$CFLAGS -fopenmp"'
            'LDFLAGS="$LDFLAGS -fopenmp"'};
    else
        fprintf('Error: Not known compiler\n')
        fprintf('You might want to edit %s if you know how to setup openmp options\n', me)
        return
    end
else
    fprintf('Please download a supported C compiler, run:\n');
    fprintf('>> mex -setup C\n');
    fprintf('then rerun %s again\n', me);
end

% invoke MEX compilation tool
mex(mexopts{:},opmopt{:},'lemire_nd_minengine.c');
mex(mexopts{:},opmopt{:},'lemire_nd_maxengine.c');

% Those mex files are no longer used, we keep for compatibility
mex(mexopts{:},R2018a_mexopts{:},'lemire_nd_engine.c');
mex(mexopts{:},R2018a_mexopts{:},'lemire_engine.c');

cd(oldpath);
