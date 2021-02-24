module NewtonRaphson

export newton_raphson

function newton_raphson(f::Function, x0::Number, fprime::Function, args::Tuple=();
                tol::AbstractFloat=1e-8, maxiter::Integer=50, eps::AbstractFloat=1e-10)
    for _ in 1:maxiter
        yprime = fprime(x0, args...)
        if abs(yprime) < eps
           warn("First derivative is zero")
           return x0
        end
        y = f(x0, args...)
        x1 = x0 - y/yprime
        if abs(x1-x0) < tol
          return x1
        end
        x0 = x1
    end
    error("Max iteration exceeded")
end

end