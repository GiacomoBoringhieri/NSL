set terminal qt persist
binwidth = 0.1
bin(x,width) = width*floor(x/width) + width/2.0
plot "pofv_iniz.dat" using (bin($1,binwidth)):(1.0) smooth freq with boxes