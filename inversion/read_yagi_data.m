function [H,d]=read_yagi_data(kernelfile,datafile)

% function [H,d]=read_yagi_data(kernelfile,datafile)
%
% takes yagi-format kernel and data files and converts it into a format
% suitable for inversions
%
% tres facile
%
% gjf, 14-feb-2007 (kissy-kissy)
% geniusinaction.com



findata=fopen(datafile);
tmpin=fscanf(findata,'%lf',[inf]);
fclose(findata);

d=tmpin(2:tmpin(1)+1);

finkernel=fopen(kernelfile);
tmpin=fscanf(finkernel,'%lf',[inf]);
fclose(finkernel);

nrows=tmpin(1);
ncols=tmpin(2);

H=reshape(tmpin(3:(nrows*ncols)+2),ncols,nrows)';


