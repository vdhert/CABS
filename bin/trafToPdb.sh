#! /bin/bash

awk -v fn=$2 '

BEGIN{
 while(getline<fn) seq[++n]=$0;
 s=1;
}

/\./{
 if(s==1) printf "MODEL%9i\n",++m;
 l=$2;
 i=0;
}

!/\./{
 for(j=1;j<NF;j+=3){
  if(i>0 && i<l-1){
   printf "ATOM%7i  CA %6s%5s   ",s,substr(seq[s],8,6),substr(seq[s],2,5);
   for(k=0;k<3;++k) printf "%8.3f",$(j+k)*0.61;
   printf "%6.2f%6.2f\n",substr(seq[s],14,3),substr(seq[s++],17,6);
  }
  ++i;
  if(s>n){
   printf "ENDMDL\n";
   s=1;
  }
 }
}

' $1
