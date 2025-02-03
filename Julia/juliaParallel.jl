using NLsolve

function f(x)
    [(x[1]+3)*(x[2]^3-7)+18,
    sin(x[2]*exp(x[1])-1)]
end

sol = nlsolve(f, [ 0.1, 1.2])
sol.zero