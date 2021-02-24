module InverseQuadraticInterpolation 
export invquadinterp

function invquadinterp(f::Function, x0::Number, x1::Number, x2::Number,
                        args::Tuple=(); xtol::AbstractFloat=1e-5,
                        ytol=2eps(Float64), maxiter::Integer=50)
    y0 = f(x0, args...)
    y1 = f(x1, args...)
    y2 = f(x2, args...)
    for _ in 1:maxiter
       x = x0*y1*y2/((y0-y1)*(y0-y2)) +
           x1*y0*y2/((y1-y0)*(y1-y2)) +
           x2*y0*y1/((y2-y0)*(y2-y1))
        # x-tolerance
        if min(abs(x-x0), abs(x-x1), abs(x-x2)) < xtol 
           return x
        end
        y = f(x, args...)
        # y-tolerance.
        if abs(y) < ytol
            return x
        end
        x0 = x1
        y0 = y1
        x1 = x2
        y1 = y2
        x2 = x
        y2 = y
    end
    error("Max iteration exceeded")
end 

end