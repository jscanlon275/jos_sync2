function [lat] = geteventlat(EEG,type)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
lat = sort([EEG.event(strcmp({EEG.event.type}, type)).latency]);

end

