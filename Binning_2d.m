function [ out_matrix ] = Binning_2d( in_matrix,x_bin, y_bin )
%BINNING_2D Summary of this function goes here
%   Detailed explanation goes here
    out_matrix=zeros(floor(size(in_matrix,1)/y_bin),floor(size(in_matrix,2)/x_bin));
    x_block=[ones(1,floor(size(in_matrix,1)/y_bin)-1)*y_bin,y_bin+mod(size(in_matrix,1),y_bin)];
    y_block=[ones(1,floor(size(in_matrix,2)/x_bin)-1)*x_bin,x_bin+mod(size(in_matrix,2),x_bin)];
    binned_matrix=mat2cell(in_matrix,x_block,y_block);

    for iii=1:size(out_matrix,1)
        for jjj=1:size(out_matrix,2)
            out_matrix(iii,jjj)=sum(sum(binned_matrix{iii,jjj}));
        end
    end

end

