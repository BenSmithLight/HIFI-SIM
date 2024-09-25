function [ wFilter] = WienerFilterW2_3D(sp)
    wFilter.sp=sp;
    wFilter.wDenom=updateCache(sp);
end

function [wDenom] = updateCache(sp)
w=sp.imgSize;
h=sp.imgSize;
wDenom=zeros(2*h,2*w);
for d=1:sp.nrDirs
    for b=1:sp.nrBands
        wd=wDenom;
        [wDenom]=addWienerDenominator_3D(wd,sp,d,b);
    end
end
end


