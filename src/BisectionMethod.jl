module BisectionMethod 
export bisection_method

function bisection_method(f::Function, a::Number, b::Number;
                          tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]") 
    i = 0 
    local c
    while b-a > tol
       i += 1
       i != maxiter || error("Max iteration exceeded")
       c = (a+b)/2
       fc = f(c)
       if fc == 0
          break
       elseif fa*fc > 0
          a = c # root is in the right half of [a,b].
          fa = fc
       else
          b = c # root is in the left side of [a,b].
       end
    end
    return c

end

end