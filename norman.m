function [image_out] = norman(image_in, TH_lo,TH_up)

ave_int = nanmean(image_in(:));
image_in(find(image_in >= TH_up*ave_int)) = TH_up*ave_int;
image_in(find(image_in <= TH_lo*ave_int)) = TH_lo*ave_int;
image_out = mat2gray(image_in);

end
