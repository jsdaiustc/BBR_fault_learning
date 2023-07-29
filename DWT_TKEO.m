function Hia=DWT_TKEO(ia)
[ia_DWT, L] = wavedec(ia, 9, 'db45');
[d3]=wrcoef('d',ia_DWT, L,'db45',3);
[d4]=wrcoef('d',ia_DWT, L,'db45',4);
ia = d3+d4;
Hia= ia(2:end-1).^2-ia(1:end-2).*ia(3:end);
Hia=(Hia-mean(Hia))/mean(Hia);