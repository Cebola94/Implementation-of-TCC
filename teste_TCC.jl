using FactCheck, Krylov
include("interior_point.jl")

facts("min x1 + 2x2 s.a x1 + x2 = 1, x1,x2 >= 0") do
    A = [1 2]
    b = 1
    c = [1;1]
    x, λ, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = false)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("min x1 + 2x2 s.a 2x1 + x2 = 0.5, x1,x2 >= 0") do
    A = [2 1]
    b = 0.5
    c = [1;2]
    x, λ, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = false)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("min x1 + x2 s.a x1 + 2x2 = 3, 3x1 + x2 = -5, x1,x2 >= 0") do
    A = [1 2;3 1]
    b = [3,3]
    c = [1;0]
    x, λ, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = false)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("min x1 - x2 s.a 2x1 + x2 = 4;x1 + 3x2 = 2, x1,x2 >= 0") do
    A = [2 1;1 3]
    b = [4;2]
    c = [1;-1]
    x, λ, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = false)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("min x1 + 2x2 s.a x1 + x2 = 3;x1 + 2x2 = 2, x1,x2 >= 0") do
    A = [2 1;1 2]
    b = [3;2]
    c = [1;2]
    x, λ, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = false)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end
