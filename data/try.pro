ha =  [17.8995  ,    18.6627   ,   4.27664  ,    3.25920]   
hea = [14.9968  ,    18.8913   ,   4.68834  ,    4.32324]
;ha =  [ 12.9788   ,   16.0474   ,   3.47310  ,    3.34094 ]
;hea = [ 22.2881   ,   19.5259   ,   4.73076  ,    4.33493 ]
lia = [ 11.5225   ,   18.4206   ,   4.51023  ,    4.95009 ]
bea = [ 7.90157   ,   18.2806   ,   4.34940  ,    5.44576 ]
ba =  [ 8.58667   ,   18.0113   ,   4.24663  ,    5.71902 ]
ca =  [ 9.31262   ,   17.8741   ,   4.17401  ,    5.95372 ]
na =  [ 7.62292   ,   17.7093   ,   4.07259  ,    6.212   ]
oa =  [ 30.5089   ,   22.6961   ,   5.08210  ,    6.38208 ]	
                                                                              
;hb =  [ -3.00245  ,    16.1024  ,    3.53766 ,   -0.813012]
;heb = [ -4.09955  ,    18.6471  ,    4.51866 ,   -0.725887]
hb =  [-3.67605   ,   18.0410   ,  4.17601  ,  -0.839795]
heb = [-3.51482   ,   18.9156   , 4.67875  ,  -0.758665]
lib = [ -2.58087  ,    18.4229  ,    4.50744 ,   -0.715188]
beb = [ -1.70145  ,    18.2918  ,    4.34880 ,   -0.699983]
bb =  [ -1.83272  ,    18.0244  ,    4.24983 ,   -0.664671]
cb =  [ -1.97643  ,    17.8867  ,    4.17913 ,   -0.639631]
nb =  [ -1.58207  ,    17.7210  ,    4.07723 ,   -0.632859]
ob =  [ -6.33045  ,    22.7251  ,    5.09112 ,   -0.615300]

temps0 = [0.05, 0.10, 0.30, 0.50, 0.75, 1.0, 1.25, 1.50, 2.0, 3.0, 5.0, 10.0]
temps0 = temps0*10000.

recs = [  2.201124,   2.011993    ,    1.705094    ,    1.552790    ,    1.425534    ,    1.331022    ,    1.255273    ,    1.191451    ,    1.087781    ,   0.9339932    ,   0.7279477    ,   0.4271614]

recs = 10.^recs

temps=temps0
log10Ne = alog10(100.)
zfit = (log10Ne-HA(1))/HA(2)
Afit = HA(0)*exp(-zfit*zfit/2.)+HA(3)
zfit = (log10Ne-Hb(1))/Hb(2)
Bfit = HB(0)*exp(-zfit*zfit/2.)+HB(3)
Hbeta = fltarr(12)
for i = 0, 11 do begin
    Hbeta(i) = 10.^(Afit+Bfit*alog10(temps(i)))
endfor


Log10Ne = alog10(100.)
zfit = (log10Ne-HEA(1))/HEA(2)
Afit = HEA(0)*exp(-zfit*zfit/2.)+HEA(3)
zfit = (log10Ne-HEb(1))/HEb(2)
Bfit = HEB(0)*exp(-zfit*zfit/2.)+HEB(3)
Hebeta = fltarr(12)
for i = 0, 11 do begin
    Hebeta(i) = 10.^(Afit+Bfit*alog10(temps(i)))
endfor

stop

temps = temps0/9.
log10Ne = alog10(100.)
zfit = (log10Ne-LIA(1))/LIA(2)
Afit = LIA(0)*exp(-zfit*zfit/2.)+LIA(3)
zfit = (log10Ne-LIb(1))/LIb(2)
Bfit = LIB(0)*exp(-zfit*zfit/2.)+LIB(3)
LIbeta = 10.^(Afit+Bfit*alog10(temps))

temps = temps0/16.
log10Ne = alog10(100.)
zfit = (log10Ne-BEA(1))/BEA(2)
Afit = BEA(0)*exp(-zfit*zfit/2.)+BEA(3)
zfit = (log10Ne-BEb(1))/BEb(2)
Bfit = BEB(0)*exp(-zfit*zfit/2.)+BEB(3)
BEbeta = 10.^(Afit+Bfit*alog10(temps))

temps = temps0/25.
log10Ne = alog10(100.)
zfit = (log10Ne-BA(1))/BA(2)
Afit = BA(0)*exp(-zfit*zfit/2.)+BA(3)
zfit = (log10Ne-Bb(1))/Bb(2)
Bfit = BB(0)*exp(-zfit*zfit/2.)+BB(3)
Bbeta = 10.^(Afit+Bfit*alog10(temps))

temps = temps0/36.
log10Ne = alog10(100.)
zfit = (log10Ne-CA(1))/CA(2)
Afit = CA(0)*exp(-zfit*zfit/2.)+CA(3)
zfit = (log10Ne-Cb(1))/Cb(2)
Bfit = CB(0)*exp(-zfit*zfit/2.)+CB(3)
Cbeta = 10.^(Afit+Bfit*alog10(temps))

temps = temps0/49.
log10Ne = alog10(100.)
zfit = (log10Ne-NA(1))/NA(2)
Afit = NA(0)*exp(-zfit*zfit/2.)+NA(3)
zfit = (log10Ne-Nb(1))/Nb(2)
Bfit = NB(0)*exp(-zfit*zfit/2.)+NB(3)
Nbeta = 10.^(Afit+Bfit*alog10(temps))

temps = temps0/64.
log10Ne = alog10(100.)
zfit = (log10Ne-OA(1))/OA(2)
Afit = OA(0)*exp(-zfit*zfit/2.)+OA(3)
zfit = (log10Ne-Ob(1))/Ob(2)
Bfit = OB(0)*exp(-zfit*zfit/2.)+OB(3)
Obeta = 10.^(Afit+Bfit*alog10(temps))

@colors
;temps=temps0
plot,  temps, hbeta    , yrange=[0.,500.]
oplot, temps, hebeta  , color = red
stop
oplot, temps, libeta  , color = green
oplot,  temps, bebeta , color = blue
oplot, temps, bbeta   , color = cyan
oplot, temps, cbeta   , color = magenta
oplot,  temps, nbeta  , color = yellow
oplot, temps, obeta  , color = brown


stop 
end
