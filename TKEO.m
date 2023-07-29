function Hia=TKEO(ia)

Hia= ia(2:end-1).^2-ia(1:end-2).*ia(3:end);



Hia=(Hia-mean(Hia))/mean(Hia);
% figure
% plot(Hia)