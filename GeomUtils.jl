module GeomUtils
using Base.Math
    export Vec2D
    struct Vec2D
        x::Float64
        y::Float64
    end

    Base.:+(a::Vec2D, b::Vec2D) = Vec2D(a.x + b.x, a.y + b.y)
    Base.:-(a::Vec2D, b::Vec2D) = Vec2D(a.x - b.x, a.y - b.y)
    Base.:*(a::Vec2D, b) = b isa Float64 ? Vec2D(a.x * b, a.y * b) : a.x * b.x + a.y * b.y
    Base.:-(a::Vec2D) = Vec2D(-a.x, -a.y)
    Base.:/(a::Vec2D, b) = Vec2D(a.x/b, a.y/b)
    Base.:*(a::Matrix{Float64},b::Vec2D) = Vec2D(b.x*a[1]+b.y*a[3] , b.x*a[2]+b.y*a[4])

    export magnitude
    function magnitude(a::Vec2D)::Float64
        return hypot(a.x, a.y)
    end

    export distance
    function distance(a::Vec2D, b::Vec2D)::Float64
        return magnitude(b - a)
    end

    export norm
    function norm(a::Vec2D)::Vec2D
        return a / magnitude(a)
    end

    export rotate
    function rotate(a::Vec2D, angle::Float64)::Vec2D
        m = magnitude(a)
        n = norm(a)
        C = cos(angle)
        S = sin(angle)
        R = [C -S; S C]

        return (R*n)*m
    end

    export turn_right
    function turn_right(a::Vec2D)::Vec2D
        X = a.y
        Y = -a.x
        return Vec2D(X, Y)
    end

    export turn_left
    function turn_left(a::Vec2D)::Vec2D
        return -turn_right(a)
    end

    export angle_between
    function angle_between(a::Vec2D, b::Vec2D)::Float64
        crossprod = a.x * b.y - b.x * a.y
        dotprod = a.x * b.x + a.y * b.y
        return atan(crossprod,dotprod)
    end

    export rotate_towards
    function rotate_towards(a::Vec2D, b::Vec2D)::Vec2D
        diff = norm(b - a)
        return diff*magnitude(a)
    end

    Base.show(io::IO, a::Vec2D) = print(io, "($(a.x), $(a.y))")

    export lerp
    function lerp(v1::Vec2D, v2::Vec2D, t::Float64)::Vec2D
        return v1*(1.0-t)+v2*t
    end

    export interpolate
    function interpolate(vs::Vector{Vec2D}, t::Float64)::Vec2D
        if length(vs) == 1
            return vs[1]
        end
        
        vls = fill(Vec2D(0.0,0.0), length(vs) - 1)
        
        for i in 2:length(vs)
            vls[i - 1] = lerp(vs[i], vs[i - 1], t)
        end
        
        return interpolate(vls, t)
    end

    export barycentric
    function barycentric(P1::Vec2D, P2::Vec2D, P3::Vec2D, P::Vec2D)::Tuple{Float64,Float64,Float64}

        P1x, P1y = P1.x, P1.y
        P2x, P2y = P2.x, P2.y
        P3x, P3y = P3.x, P3.y
        Px, Py = P.x, P.y

        dg = ((P1y - P2y) * P3x + (P2x - P1x) * P3y+ P1x * P2y - P2x * P1y)
        db = ((P1y - P3y) * P2x + (P3x - P1x) * P2y+ P1x * P3y- P3x * P1y)

        if db == 0 || dg == 0
            return (-1,-1,-1)
        end

        gamma = ((P1y- P2y) * Px + (P2x - P1x) * Py+ P1x * P2y- P2x * P1y) / dg

        beta = ((P1y- P3y) * Px + (P3x - P1x) * Py+ P1x * P3y- P3x * P1y) / db

        alpha = 1 - beta - gamma

        return alpha, beta, gamma
    end

    export solveEq2d
    function solveEq2d(a::Float64,b::Float64,c::Float64)::Tuple{Float64,Float64}
        Delta = b^2 - 4*a*c
        if Delta < 0 
            throw(DomainError(Delta,"Delta < 0 implies complex solutions."))
        end
        # Angular coefficient
        x1 = (-b + sqrt(Delta))/(2*a)
        x2 = (-b - sqrt(Delta))/(2*a)
        return x1,x2
    end

    export lineEq
    function lineEq(m::Float64,q::Float64)::var"#f#1"{Float64, Float64}
        return f(x) = m*x + q
    end

    export tangents_to_circle
    function tangents_to_circle(xp::Float64,yp::Float64,xo::Float64,yo::Float64,r::Float64)::Tuple{var"#f#1"{Float64, Float64}}
        # Computing tangent lines to circle passing through the point d1.pos
        dx = xo-xp
        dy = yo-yp
        a = dx^2 - r^2
        b = 2*dx*dy
        c = dy^2 - r^2
                    
        # Angular coefficient
        m1,m2 = solveEq2d(a,b,c)
        
        # Intersection with y axis
        q1 = yp-m1*xp
        q2 = yp-m2*xp
        return lineEq(m1,q1),lineEq(m2,q2)
    end

    export leastSq
    function leastSq(points::Vector{Vec2D})
        
        Xavg = sum([p.x for p in points])/length(points)

        Yavg = sum([p.y for p in points])/length(points)

        Sxy = sum([(p.x-Xavg)*(p.y-Yavg) for p in points])
        Sx2 = sum([(p.x-Xavg)^2 for p in points])

        m = Sxy/Sx2
        q = Yavg - m*Xavg

        return lineEq(m,q)
    end

end