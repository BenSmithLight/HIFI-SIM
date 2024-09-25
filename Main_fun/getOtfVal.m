function [ val ]= getOtfVal(ret,cycl)
    mask=cycl>ret.cutoff;
    cycl(mask)=0; 

    pos=cycl./ret.cyclesPerMicron;
    cpos=pos+1;
    lpos=floor(cpos);
    hpos=ceil(cpos);
    f=(cpos-lpos); 

    retl=ret.vals(lpos).*(1-f);
    reth=ret.vals(hpos).*f;
    val=retl+reth;
    val(mask)=0;
end

