module BrentsMethod
export brents_method

function brents_method(f::Function, x0::Number, x1::Number, args::Tuple=();
         xtol::AbstractFloat=1e-7, ytol=2eps(Float64),
         maxiter::Integer=50)
    EPS = eps(Float64)
    y0 = f(x0, args...)
    y1 = f(x1, args...)
    if abs(y0) < abs(y1) 
       # swap lower and upper bounds
       x0, x1 = x1, x0
       y0, y1 = y1, y0
    end
    x2 = x0
    y2 = y0
    x3 = x2
    bisection = true
    for _ in 1:maxiter
       # x-tolerance
       if abs(x1 - x0) < xtol
          return x1
       end

       # Use inverse quadratic interpolation if f(x0) != f(x1) != f(x2)
       # and linear interpolation (secant method) otherwise.
       if abs(y0 - y2) > ytol && abs(y1 - y2) > ytol
           x = x0*y1*y2/((y0-y1)*(y0-y2)) +
                x1*y0*y2/((y1-y0)*(y1-y2)) +
                x2*y0*y1/((y2-y0)*(y2-y1))
       else
            x = x1 - y1 * (x1-x0)/(y1-y0)
       end

       # Use bisection method if satisfies the conditions.
       delta = abs(2EPS*abs(x1))
       min1 = abs(x - x1)
       min2 = abs(x1 - x2)
       min3 = abs(x2 - x3)
       if (x < (3x0 + x1)/4 && x > x1) ||
          (bisection && min1 >= min2/2) ||
          (!bisection && min1 >= min3/2) ||
          (bisection && min2 < delta) ||
          (!bisection && min3 < delta)
          x = (x0 + x1)/2
          bisection = true
       else
          bisection = false
       end

       y = f(x, args...)
       # y-tolerance.
       if abs(y) < ytol
          return x
       end
       x3 = x2
       x2 = x1
       if sign(y0) != sign(y)
          x1 = x
          y1 = y 
       else
          x0 = x
          y0 = y
       end
       if abs(y0)  < abs(y1)
          #swap lower and upper bounds
          x0, x1 = x1, x0
          y0, y1 = y1, y0
       end 
    end
    error("Max iteration exceeded")
       
end

end 