### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 6cc00c6e-ddaa-11ec-0ada-4302207df1e9
begin
	import Pkg;
	Pkg.activate();

	using Plots;
	using OrderedCollections;
end;

# ╔═╡ 0ec70f3c-9dc6-480c-8690-3156c8f40a37
md"""
# Лабораторная работа #1. Метод итераций

Выполнил: Боровский Илья
Вариант: 6
"""

# ╔═╡ d135bf80-3b9a-4435-b59a-4b64f1549d78
md"""
$x^3 - 1.5x^2 + 0.58x - 0.057 = 0$
"""

# ╔═╡ e11af52e-e9a4-4804-a157-37747dae9a70
f(x) = x^3 - 1.5 * x^2 + 0.58 * x - 0.057

# ╔═╡ 04f3f6b5-88ad-4971-8ea2-74a6ff8b9209
plot(
	-1:.1:2,
	f.(-1:.1:2);
	label=:none
)

# ╔═╡ adeb96ec-9346-4805-bec7-3d5388639035
md"""
Возьмём отрезок $[-1, 2]$ в качестве $[a, b]$.

Приведём $x^3 - 1.5x^2 + 0.58x - 0.057 = 0$ к канонической форме, выразив $x$.

$- x^3 + 1.5x^2 + 0.057 = 0.58x$

$\frac{- x^3 + 1.5x^2 + 0.057}{0.58} = x$

имеем

$φ(x) = \frac{- x^3 + 1.5x^2 + 0.057}{0.58}$

$φ'(x) = \frac{- 3*x^2 + 3x}{0.58}$
"""

# ╔═╡ bebb9185-b601-4f64-ae31-9331bf227db8
φ(x) = (-x^3 + 1.5x^2 + 0.057) / 0.58

# ╔═╡ 3ab565f7-9ddb-4c1f-97ed-88a554ae8288
dφ(x) = (-3x^2 + 3x) / 0.58

# ╔═╡ 23ebaf18-0795-47c5-a2f8-aefe3e80cdb6
plot(
	-.2:.01:.4,
	abs.(dφ.(-.2:.01:.4)); 
	label=:none
)

# ╔═╡ 2d0be0ea-2c73-40f3-b788-c4d7e0eab2ae
hline!([1]; label="1")

# ╔═╡ 9bf4fec6-a92d-4af1-9452-78abd30acea4
vline!([-.085, .105]; label=["x₀ - δ", "x₀ + δ"])

# ╔═╡ 380daf34-f705-47ae-b967-0e73c574a09b
hline!([.5]; label="q")

# ╔═╡ 4d137535-d929-4751-9f89-c86eb781350e
vline!([.01]; label="x₀")

# ╔═╡ 3350461e-37bc-4aa3-a0df-2220d2bea4af
md"""
В качестве $|x - x₀| ≤ δ$ отрезок $[-0.065, 0.081]$, тогда $x₀ = 0.01, δ = 0.095, q = 0.5$.

Вычислим $m$:

$m ≤ δ * (1 - q)$

$m ≤ 0.095 * (1 - 0.5) = 0.0475$
"""

# ╔═╡ b7b2bc19-b5da-421f-aef7-7173246e6bed
md"""
## Вычисление корня
"""

# ╔═╡ da7ee5ae-677a-441c-829f-a45f4b00d167
x = [.01]

# ╔═╡ 7611737c-360b-4f96-8220-f7294b674462
push!(x, φ(x[end]))

# ╔═╡ 4ff99569-3803-45c3-bc31-a39478aaf60f
while abs(x[end] - x[end-1]) > 0.00001
	push!(x, φ(x[end]))
end

# ╔═╡ e01cdef2-b31b-4a48-a39b-5b4915491d12
OrderedDict(x .=> f.(x))

# ╔═╡ 813d79cd-96cd-4917-81bc-56564c1f54f4
f(x[end])

# ╔═╡ Cell order:
# ╠═6cc00c6e-ddaa-11ec-0ada-4302207df1e9
# ╟─0ec70f3c-9dc6-480c-8690-3156c8f40a37
# ╟─d135bf80-3b9a-4435-b59a-4b64f1549d78
# ╠═e11af52e-e9a4-4804-a157-37747dae9a70
# ╠═04f3f6b5-88ad-4971-8ea2-74a6ff8b9209
# ╟─adeb96ec-9346-4805-bec7-3d5388639035
# ╠═bebb9185-b601-4f64-ae31-9331bf227db8
# ╠═3ab565f7-9ddb-4c1f-97ed-88a554ae8288
# ╠═23ebaf18-0795-47c5-a2f8-aefe3e80cdb6
# ╠═2d0be0ea-2c73-40f3-b788-c4d7e0eab2ae
# ╠═9bf4fec6-a92d-4af1-9452-78abd30acea4
# ╠═380daf34-f705-47ae-b967-0e73c574a09b
# ╠═4d137535-d929-4751-9f89-c86eb781350e
# ╟─3350461e-37bc-4aa3-a0df-2220d2bea4af
# ╟─b7b2bc19-b5da-421f-aef7-7173246e6bed
# ╠═da7ee5ae-677a-441c-829f-a45f4b00d167
# ╠═7611737c-360b-4f96-8220-f7294b674462
# ╠═4ff99569-3803-45c3-bc31-a39478aaf60f
# ╠═e01cdef2-b31b-4a48-a39b-5b4915491d12
# ╠═813d79cd-96cd-4917-81bc-56564c1f54f4
