@#include "rbc_inc.mod"

limits!("k", min = 0.05, max = 0.08);
%limits!("k", min = 0.01, max = 0.13);
limits!("z", min = -4*sig/sqrt(1-rho^2), max = 4*sig/sqrt(1-rho^2));