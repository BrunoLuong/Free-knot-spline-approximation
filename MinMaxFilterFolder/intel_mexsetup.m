function intel_mexsetup()
% intel_mexsetup()
% Setup mex for intel oneAPI C compiler

intelCompilerPath = 'C:\Program Files (x86)\Intel\oneAPI\compiler\2021.2.0\windows\';
mexoptspath = fullfile(matlabroot,'bin','win64','mexopts');
oldpath = cd(mexoptspath);
setenv('ICPP_COMPILER21', intelCompilerPath);
mex -setup:intel_c_21_vs2019.xml -v c

cd(oldpath)

end % intel_mexsetup
