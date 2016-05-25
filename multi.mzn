include "adj.dzn";
int: nAminos = npos*nProteins;



array[1..npos,1..nAminos] of var 0..1: s;


forall(i in 1..npos)(
forall (k in 1..nProteins)(
    sum(j in (k-1)*npos+1..k*npos)(
     s[i,j] 
   )<=1
)
);


constraint 
forall (j in 1..nAminos)(
   sum(i in 1..npos)(
     s[i,j] 
   )<=1
);



solve maximize sum(i in 1..npos,j in 1..nAminos,k in 1..nAminos)(
    adj[j,k] * s[i,j] * s[i,k]
);




