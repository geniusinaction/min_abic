function [staname,stasigw,stasamp,statstp,stalocn]=read_stations(infile)

% function [staname,stasigw,stasamp,statstp,stalocn]=read_stations(infile)
%
% read in information about stations. useful for plotting, not much else
%
% gjf, 14-feb-2007 (mwah)
% geniusinaction.com

fin=fopen(infile,'r');
input=textscan(fin,'%s%s%f%d%f%f%f');
fclose(fin);

staname(:,1)=input{1};  % station name code
staname(:,2)=input{2};  % station component
stasigw=input{3};       % sigw, whatever the hell that is
stasamp=input{4};       % number of time samples
statstp=input{5};       % time step (sample interval)
stalocn(:,1)=input{7};  % station distance
stalocn(:,2)=input{6};  % station azimuth


