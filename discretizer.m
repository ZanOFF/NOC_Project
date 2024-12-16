function [xc_out, yc_out, xin_out, yin_out, xout_out, yout_out] = discretizer(xc_in, yc_in, xin_in, yin_in, xout_in, yout_in, Ss_new, Ss_old)
    
    j = 1;
    for i = 1:Ss_new/Ss_old:length(xc_in)
        xc_out(j) = xc_in(i);
        j = j+1;
    end
    j = 1;
    for i = 1:Ss_new/Ss_old:length(yc_in)
        yc_out(j) = yc_in(i);
        j = j+1;
    end
    j = 1;
    for i = 1:Ss_new/Ss_old:length(xin_in)
        xin_out(j) = xin_in(i);
        j = j+1;
    end
    j = 1;
    for i = 1:Ss_new/Ss_old:length(yin_in)
        yin_out(j) = yin_in(i);
        j = j+1;
    end
    j = 1;
    for i = 1:Ss_new/Ss_old:length(xout_in)
        xout_out(j) = xout_in(i);
        j = j+1;
    end
    j = 1;
    for i = 1:Ss_new/Ss_old:length(yout_in)
        yout_out(j) = yout_in(i);
        j = j+1;
    end
end
