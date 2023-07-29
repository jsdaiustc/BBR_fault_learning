function Hia=TOEO(ia)
Hia=ia(2:end-2).*ia(3:end-1)-ia(1:end-3).*ia(4:end);
Hia=(Hia-mean(Hia))/mean(Hia);