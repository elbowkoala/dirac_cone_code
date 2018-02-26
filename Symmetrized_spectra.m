function [ spectra_sym,kaxis_sym ] = Symmetrized_spectra( spectra, Kaxis )
%SYMMETRIZED_SPECTRA Summary of this function goes here
%   Detailed explanation goes here
sum_list=[];
spectra_list={};

for kpl=round(size(spectra,1)/2)-5:round(size(spectra,1)/2)+5
    asyspectra=zeros(min([abs(kpl),abs(size(spectra,1)-kpl)]),size(spectra,1));
    asyspectra=spectra(kpl+1:kpl+size(asyspectra,1),:)-flipud(spectra(kpl-size(asyspectra,1)+1:kpl,:));
    sym_spectra=(spectra(kpl+1:kpl+size(asyspectra,1),:)+flipud(spectra(kpl-size(asyspectra,1)+1:kpl,:)))/2;
    
    sum_list(end+1,:)=[kpl,sum(sum(asyspectra(:).^2))];
    spectra_list(end+1)={sym_spectra};
end

[~,m_kpl] = min(sum_list(:,2));
k_zero = sum_list(m_kpl,1);
spectra_sym = [flipud(spectra_list{m_kpl}(1:end-1,:));spectra_list{m_kpl}];
k_length = (size(spectra_sym,1)-1)/2;

kaxis_sym = [Kaxis(k_zero-k_length:k_zero-1)-Kaxis(k_zero);-flipud(Kaxis(k_zero-k_length:k_zero)-Kaxis(k_zero))];

end

