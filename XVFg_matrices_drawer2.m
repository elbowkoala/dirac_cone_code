function XVFg_matrices = XVFg_matrices_drawer2(VFs,k_gamma,e_gamma,matsize)

%matsize = size(cone_mod);
%gamma = 2;

matsize1 = matsize(1);
matsize2 = matsize(2);
XVFg_matrices = zeros(matsize1,matsize2,length(VFs));

x = 1:matsize1;
x_0 = round(matsize1/2);
y_0 = round(matsize2/2);

for v_i = 1:length(VFs)
    
    VF = VFs(v_i);
    y1 = -VF*(x - x_0) + y_0;
    y2 = +VF*(x - x_0) + y_0;
    XVF_matrix1 = zeros(matsize1,matsize2);
    XVF_matrix2 = zeros(matsize1,matsize2);
    for y = 1:matsize2
        y1_srch = abs(y1 - y);
        y2_srch = abs(y2 - y);
        x1_i = y1_srch==min(y1_srch);
        x2_i = y2_srch==min(y2_srch);
        x01 = x(x1_i);
        x02 = x(x2_i);
        XVF_matrix1(x01,y) = XVF_matrix1(x01,y) + 1;
        XVF_matrix2(x02,y) = XVF_matrix2(x02,y) + 1;
    end
    
    for xn = 1:matsize1
        yn1 = round(y1(xn));
        if yn1 >= 1 && yn1 <= matsize2
            XVF_matrix1(xn,yn1) = XVF_matrix1(xn,yn1)+1;
        end
        yn2 = round(y2(xn));
        if yn2 >= 1 && yn2 <= matsize2
            XVF_matrix2(xn,yn2) = XVF_matrix2(xn,yn2)+1;
        end
    end
    XVF_matrix1(XVF_matrix1>0) = 1;
    XVF_matrix2(XVF_matrix2>0) = 1;
    
    XVF_matrix = XVF_matrix1 + XVF_matrix2;
     

    kernel_k = zeros(41,1);
    kermidrow = length(kernel_k)/2;
    for kerk = 1:length(kernel_k)
        kernel_k(kerk) = 1 / ((kerk - kermidrow)^2+k_gamma^2);
    end
    kernel_kspread = (k_gamma/(2*pi))*repmat(kernel_k,[1,41]);
    
    kernel_e = zeros(1,41);
    kermidcol = length(kernel_e)/2;
    for kere = 1:length(kernel_e)
        kernel_e(kere) = 1 / ((kere - kermidcol)^2+e_gamma^2);
    end
    kernel_espread = (e_gamma/(2*pi))*repmat(kernel_e,[41,1]);
    
    kernel_bothspread = kernel_kspread .* kernel_espread;

    XVFg_matrix = imfilter(XVF_matrix,kernel_bothspread);    
 
    
    XVFg_matrices(:,:,v_i) = XVFg_matrix;
end

% figure, subplot(121), imagesc(XVFg_matrices{1}+XVFg_matrices{end}), axis xy
% subplot(122), imagesc(XVFg_matrices{end}), axis xy

end
