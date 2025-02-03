using Plots

# Define the Runge-Kutta 4th Order Method
function runge_kutta_4(f, y0, t0, tf, h)
    n = Int((tf - t0) / h)
    t = t0
    y = y0

    for i in 1:n
        k1 = h * f(t, y)
        k2 = h * f(t + 0.5h, y + 0.5k1)
        k3 = h * f(t + 0.5h, y + 0.5k2)
        k4 = h * f(t + h, y + k3)
        y += (k1 + 2k2 + 2k3 + k4) / 6
        t += h
    end

    return y
end

#function euler_basic(f, y0, t0, tf, h)
#end

# reduced mass of the spaceship
mu = 1

# first of four hamiltonian equations for the
# spaceship's motion in the Earth-Moon cm system
#r(t, pr,phi,pphi) = pr/mu
r(t, pr) = pr/mu
#phi(t,r,pr,pphi) = pphi/(mu*r^2)
#pr()

r0 = 0.5
t0 = 0.0
tf = 12.0
h = 0.1

timeinterval = t0:h:tf
fx = log.(timeinterval)
fy = cos.(timeinterval)

plot(fx, fy, label="Parametric Plot")

display(plot(fx, fy, label="Parametric Plot"))

solution = runge_kutta_4(r, r0, t0, tf, h)
println("The solution at t = $tf is y = $solution")