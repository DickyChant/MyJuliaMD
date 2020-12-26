using LinearAlgebra

const Length = 15.0
const particles = 2
const dt = 0.00001
const sigma = 1.0
const eps = 1.0
const mass = 1.0
const temp = 1.0
const delta = 0.001

function minDist(x1,x2)
    if ( x1 - x2 ) > (Length / 2)
        return x1 - x2 - Length
    elseif ( x1 - x2 ) < (Length / 2)
        return x1 - x2 + Length
    else
        return x1 - x2
    end
end 

function calcForce(xi,xj,r)
    xDiff = minDist(xi,xj)
    sig6 = (sigma^6)/(r^6)
    return 6*sig6*eps*(2*sig6-1)*(xDiff)/(r^2)
end

function boudaryBounce(x1)
    x = x1
    for i in 1:particles
        x[i] = boudarybounce(x[i])
    end
end

function boudarybounce(x1)
    if (x1>Length)
        return x1 - Length * (x1 ÷ Length)
    elseif x1 < 0
        return x1 + Length * ((-x1) ÷ Length + 1)
    else 
        return x1
    end
end

function calcKE(vx,vy,vz)
    return 0.5 * mass * (vx.^2+vy.^2+vz.^2)
end

function KE_all(vxs,vys,vzs)
    KE = calcKE(vxs,vys,vzs)
    return sum(KE)
end

function PE_lj(r)
    sig6 = sigma^6 .* r.^(-6)
    PE = eps.*sig6.*(sig6-1)
end

function PE_all(xs,ys,zs)
    p_tot = 0.
    for i in 1:length(xs) - 1
        for j in i:length(xs)
            xmin = minDist(xs[i],xs[j])
            ymin = minDist(ys[i],ys[j])
            zmin = minDist(zs[i],zs[j])
            r = sqrt(xmin^2+ymin^2+zmin^2)
            p_tot += PE_lj(r)
        end
    end
    return p_tot
end

function updatePos(xs,vxs,fxs)
        boudaryBounce(xs .+ vxs .* dt + fxs .* dt^2 ./ 2 ./ mass)
end

function updateVelosity(vxs,fxs)
    vxs .+ fxs ./ 2 ./ mass .* dt
end

function updateForce(x,y,z)
    fxs = zeros(particles)
    fys = zeros(particles)
    fzs = zeros(particles)
    for i = 1:particles-1
        fx[i] = 0
        for j = i:particles
            xDiff = minDist(x[i],x[j])
            yDiff = minDist(y[i],y[j])
            zDiff = minDist(z[i],z[j])

            r = sqrt(xDiff^2+yDiff^2+zDiff^2)

            fxs[i] += calcForce(x[i],x[j],r)
            fxs[j] -= calcForce(x[i],x[j],r)
            fys[i] += calcForce(y[i],y[j],r)
            fys[j] -= calcForce(y[i],y[j],r)
            fzs[i] += calcForce(z[i],z[j],r)
            fzs[j] -= calcForce(z[i],z[j],r)


        end
    end
    return fxs,fys,fzs

end

vx = [1 2 3]
r = sqrt(sum(vx.^2))
println(PE_lj(r))