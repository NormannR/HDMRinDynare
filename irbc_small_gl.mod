@#include "irbc_small_inc"

@#for j in 1:N
  limits!("k_@{j}", min = 0.7, max = 1.3);
  limits!("a_@{j}", min = -0.22, max = 0.22);
@#endfor