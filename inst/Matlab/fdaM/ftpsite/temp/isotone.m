function mony = isotone(y)
%  ISOTONE computes the isotonic regression function,
%  the monotonic polygonal line that minimizes error
%  sum of squares.

% Last modified 20 July 2006

  n = length(y);
  mony = y;
  eb   = 0;
  while (eb < n)  
    bb = eb + 1;
    eb = bb;
    while eb < n && mony(bb) == mony(eb+1) 
      eb = eb + 1;
    end
    poolflg = -1;
    while poolflg ~= 0
      if eb >= n || mony(eb) <= mony(eb+1)
        poolflg = 1;
      end
      if poolflg == -1   
        br = eb+1;
        er = br;
        while er < n && mony(er+1) == mony(br) 
          er = er + 1;
        end
        pmn = (mony(bb)*(eb-bb+1) + mony(br)*(er-br+1))/(er-bb+1);
        eb = er;
        mony(bb:eb) = pmn;
        poolflg = 1;
      end
      if poolflg == 1  
        if bb <= 1 || mony(bb-1) <= mony(bb)
          poolflg = 0;
        else
          bl = bb-1;
          el = bl;
          while bl > 1 && mony(bl-1) == mony(el) 
            bl = bl - 1;
          end
          pmn = (mony(bb)*(eb-bb+1) + mony(bl)*(el-bl+1))/(eb-bl+1);
          bb = bl;
          mony(bb:eb) = pmn;
          poolflg = -1;
        end
      end
    end
  end
 
      
     