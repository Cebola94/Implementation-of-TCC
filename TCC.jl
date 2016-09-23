##Funçã para calcular o passo "d" e o tamanho alpha
function metodo_newton_PC(B, L_xλ, cx, λ, y;tol = 1e-6, tau = 0.5)
    (m,n) = size(A)
    d = B\[L_xλ;cx;-diagm(λ)*diagm(y)*ones(m)]
    dx = d[1:n]
    dy = d[n+1:n+m]
    dλ = d[n+m+1:n+2*m]
    μ = dot(y,λ)/m
    α = 1.0
    for i = 1:m
        while y[i] +α*dy[i] <= 0 || λ[i] +α*dλ[i] <= 0
            α = α*0.9
        end
    end
    μ_aff = dot(y + α*dy,λ + α*dλ)/m
    σ = (μ_aff/μ)^3
    d = B\[L_xλ;cx;-diagm(λ)*diagm(y)*ones(m)-diagm(dλ)*diagm(dy)*ones(m)+(σ*μ)*ones(m)]
    dx = d[1:n]
    dy = d[n+1:n+m]
    dλ = d[n+m+1:n+2*m]
    I = find(dy.<0)
    J = find(dλ.<0)
    if length(I) == 0
        α1 = 0.5
    else
        α1 = minimum(-y[I]./dy[I])
    end
    if length(J) == 0
        α2 = 0.5
    else
        α2 = minimum(-λ[J]./dλ[J])
    end
    α = tau*min(α1,α2)
    return d, α
end

##Funçã principal
function pre_corr_QP(f, x0, A, b; tol = 1e-6, max_iter = 1000, max_time = 60)
    exit_flag = 0
    ∇f(x) = ForwardDiff.gradient(f, x)
    H(x) = ForwardDiff.hessian(f, x)
    (m,n) = size(A)

    x = copy(x0) # Cópia de x0
    y = A*x - b
    for i = 1:m
        if y[i] <= 0
            y[i] = 1
        end
    end
    iter = 0
    start_time = time()
    elapsed_time = 0.0
    fx = f(x)
    ∇fx = ∇f(x)
    G = H(x)
    λ = ones(m)
    L_xλ = ∇fx - A'*λ
    while norm(L_xλ) > tol && norm(A*x-b, 1) > tol
        cx = A*x - y - b
        B = [G zeros(n,m) -A'; A -eye(m) zeros(m,m);zeros(m,n) diagm(λ) diagm(y)]
        d, α = metodo_newton_PC(B ,-L_xλ ,-cx ,λ ,y)
        x += α*d[1:n]
        y += α*d[n+1:n+m]
        for i = 1:m
            if y[i] <= 0
                y[i] = 1
            end
        end
        λ += α*d[n+m+1:n+2*m]
        for i = 1:m
            if λ[i] <= 0
                λ[i] = 1
            end
        end
        fx = f(x)
        ∇fx = ∇f(x)
        G = H(x)
        L_xλ = ∇fx - A'*λ
        iter = iter + 1
        if iter >= max_iter
            exit_flag = 1
            break
        end
        elapsed_time = time() - start_time
        if elapsed_time >= max_time
            exit_flag = 2
            break
        end
    end
    return x, fx, ∇fx, exit_flag, iter, elapsed_time
end
