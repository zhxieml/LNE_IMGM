 function p = GRPdur(n)
% -------------------------------------------------------------
% Generate a random permutation p(1:n) using Durstenfeld's 
% O(n) Shuffle Algorithm, CACM, 1964. 
% See Knuth, Section 3.4.2, TAOCP, Vol 2, 3rd Ed.
% Complexity: O(n)
% USE: p = GRPdur(10^7);
% Derek O'Connor, 8 Dec 2010.  derekroconnor@eircom.net
% -------------------------------------------------------------
      p = 1:n;                  % Start with Identity permutation
  for k = n:-1:2    
      r = 1+floor(rand*k);      % random integer between 1 and k
      t    = p(k);
      p(k) = p(r);               % Swap(p(r),p(k)).
      p(r) = t;                  
  end
  return % GRPdur