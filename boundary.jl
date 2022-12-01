C1 = 0.8/365
C2 = 0.2/365
E = 0.01
G = 0.03

function calc_g(t, T = 40)
    (t * (C1 + C2) - T * C2)/T
    end

t = [0:2:40;]
g = calc_g.(t)

plot(t, g)

function calc_t(L, W0 = 0.8, W1 = 0.1)
    (W0 - W1) / (E * L)
    end

function calc_g_from_L(L, T = 40)
    calc_g(calc_t(L))
end

function calc_g_from_tau(tau)
    calc_g(calc_t(calc_L(tau)))
end

function calc_L(tau)
    L = 1 - tau^3 * G^3
    if L <= 0
        return 0
    else return L
    end
end


tau = [1:0.1:100;]
a = calc_L.(tau)
t = calc_t.(a)
g = calc_g.(t)
using QuadGK
function int_tau(input)
    quadgk(calc_g_from_tau, 0.1, input, rtol=1e-3)
end

bigL = int_tau.(tau[a.>0])
plot(tau[a.>0], first.(bigL))

function int_L(L)
    quadgk(calc_g_from_L, 0.1, L, rtol=1e-3)
end

bigL = int_L.(a)
plot(a,first.(bigL))
