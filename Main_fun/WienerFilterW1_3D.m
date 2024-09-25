function [ wFilter ] = WienerFilterW1_3D( sp )
    wFilter.sp=sp;
    wFilter.wDenom=updateCache(sp);
end

function [ wDenom ] = updateCache( sp )
    wDenom=zeros(2*sp.imgSize,2*sp.imgSize);
    for d=1:sp.nrDirs
            wDenom=wDenom+WienerDenominatorNoAtt_3D(sp,d,1)+WienerDenominatorUseAtt_3D(sp,d,2)+WienerDenominatorUseAtt_3D(sp,d,3);
    end
end


