include "adj.dzn";
int: nAminos = npos*nProteins;



array[1..npos,1..nAminos] of var 0..1: s;
array[1..npos,1..nAminos,1..nAminos] of var 0..1:y;


constraint 
forall (i in 1..nProteins)(
    sum(j in (i-1)*npos+1..i*npos)(
     s[i,j] 
   )<=1
);

constraint 
forall (j in 1..nAminos)(
   sum(i in 1..npos)(
     s[i,j] 
   )<=1
);

constraint
forall (i in 1..npos)(
  forall(j in 1..nAminos)(
   forall(k in 1..nAminos)(
       y[i,j,k] <= s[i,j]
     )
   )
);
constraint
forall (i in 1..npos)(
  forall(j in 1..nAminos)(
   forall(k in 1..nAminos)(
       y[i,j,k] <= s[i,k]
     )
   )
);

solve maximize sum(i in 1..npos,j in 1..nAminos,k in 1..nAminos)(
    adj[j,k] * y[i,j,k]
);




