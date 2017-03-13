function interior_point_method(A, b, c; tol = 1e-5,max_time = 60, max_iter = 50, verbose = true)
    (m,n) = size(A)
    e = ones(n)
    τ = 0.9
    iter = 0
    ef = 0
    el_time = 0.0
    st_time = time()
    x = e; s = e; λ = ones(m)
    #x, s, λ = start_point(A, b, c, tol, iter, el_time)
    μ = dot(x,s)/n
    L_xλ = A' * λ + s - c
    cx = A * x - b
    S = diagm(s);X = diagm(x)
    if verbose
        @printf("%4s  %9s  %9s  %13s  %4s  %8s\n", "iter", "||c(x)||", "||L_xλ||", "||c'x - b'λ||", "λ", "μ")
    end

    while norm(L_xλ) > tol || norm(cx) > tol || norm(c'*x - b'*λ) > tol
        cx_ant = cx
        S = diagm(s);X = diagm(x); M = [-L_xλ;-cx;-X * S * e]
        B = [zeros(n,n) A' eye(n); A zeros(m,m) zeros(m,n);S zeros(n,m) X]
        d = B\M
        #d,stats = cg(B,M)
        dx = d[1:n]
        dλ = d[n+1:n+m]
        ds = d[n+m+1:m+2*n]
        M = [zeros(n);zeros(m);-X * S * e]
        d = B\M
        #d,stats = cg(B,M)
        dx += d[1:n]
        dλ += d[n+1:n+m]
        ds += d[n+m+1:m+2*n]
        I = find(dx.<0)
        J = find(ds.<0)
        if length(I) == 0
            α1 = 0.5
        else
            α1 = min(1,minimum(-x[I] ./ dx[I]))
        end
        if length(J) == 0
            α2 = 0.5
        else
            α2 = min(1,minimum(-s[J] ./ ds[J]))
        end
        μ_aff = dot(x + α1 * dx,s + α2 * ds)/n
        σ = (μ_aff/μ)^3
        M = [-L_xλ;-cx;-X * S * e - diagm(dx) * diagm(ds) * e + (σ * μ) * e]
        d = B\M
        #d,stats = cg(B,M)
        dx = d[1:n]
        dλ = d[n+1:n+m]
        ds = d[n+m+1:m+2*n]
        I = find(dx.<0)
        J = find(ds.<0)
        if length(I) == 0
            α1 = 0.5
        else
            α1 = minimum(-x[I] ./ dx[I])
        end
        if length(J) == 0
            α2 = 0.5
        else
            α2 = minimum(-s[J] ./ ds[J])
        end
        α1 = min(1,τ*α1)
        α2 = min(1,τ*α2)
        x += α1 * dx
        s += α2 * ds
        λ += α2 * dλ
        μ = dot(x,s)/n
        L_xλ = A' * λ + s - c
        cx = A * x - b
        if norm(cx) <= norm(cx_ant)
            τ = min(1,1.001*τ)
        end
        iter = iter + 1
        verbose && @printf("%2d  %9.1e  %9.1e  %12.1e  %9.1e  %9.1e\n", iter, norm(cx), norm(L_xλ), norm(c'*x - b'*λ), norm(λ), μ)
        if iter >= max_iter
            ef = 1
            break
        end
        el_time = time() - st_time
        if el_time >= max_time
            ef = 2
            break
        end
    end

    return x, λ, norm(cx), norm(L_xλ), iter, ef, el_time
end

function start_point(A, b, c, tol, iter, eltime)
    (m,n) = size(A)
    e = ones(n)
    AA_t = inv(A * A')
    x = A' * AA_t * b
    x = vec(x)
    λ = AA_t * A * c
    s = c - A' * λ
    δ_x = max((-3/2) * minimum(x),0)
    δ_s = max((-3/2) * minimum(s),0)
    x = x + δ_x * e
    s = s + δ_s * e
    μ = dot(x,s)/n
    L_xλ = A' * λ + s - c
    cx = A * x - b
    if norm(cx) < tol && norm(L_xλ) < tol
       return x, s, λ
    end
    δ_x = 0.5 * dot(x,s)/dot(e,s)
    δ_s = 0.5 * dot(x,s)/dot(e,x)
    x = x + δ_x * e
    s = s + δ_s * e
    return x, s, λ
end
