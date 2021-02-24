module SecantMethod

function secant_method(f::Function, x0::Number, x1::Number, args::Tuple=();
                       tol::AbstractFloat=1e-5, maxiter::Integer=50)
    for _ in 1:maxiter
       y1 = f(x1, args...)
       y0 = f(x0, args...)
       x = x1 - y1* (x1 - x0)/(y1 - y0)
       if abs(x - x1) < tol
          return x
       end
       x0 = x1
       x1 = x
    end
    error("Max iteration exceeded")
end

end