### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ e979dd52-76fa-11ed-3ecd-d1f9941a18d1
begin
	import Pkg; Pkg.activate()

	using LinearAlgebra
	using OrderedCollections
	using SymPy
	using Plots
end

# ╔═╡ 15637745-5a70-4370-9058-553550a4c463
md"""
# Метод замены ядра на вырожденное

$φ(x) - 0.2 * \int\limits_{0}^{1} \frac{1}{10 - x*y} φ(y) 𝕕 y = 1 + x^2$
"""

# ╔═╡ 6c9f99a5-59c8-4fed-a10c-1ea6f49aa7ee
A, B = 0, 1

# ╔═╡ b58ac615-6f2b-49a2-9a41-cdb36406f165
λ = 0.2

# ╔═╡ c93395e8-b244-4a10-9a35-c3797426e3e0
K(x, y) = 1 / (10 - x*y)

# ╔═╡ f0f89c25-0424-4115-b597-645138579c19
f(x) = 1 + x^2

# ╔═╡ 5d9526a8-3194-4e1b-a76b-aa8cda69058c
md"""
## Разложение в ряд Тейлора в (1/2, 1/2)
"""

# ╔═╡ 367911ea-f525-4b4a-bc3c-6a9fa3e41f56
r = 4

# ╔═╡ 80516d38-acc0-408f-8bd3-bd97dbc532c2
@syms 𝕩, 𝕪

# ╔═╡ 292cf495-f860-4615-907d-7b06b5de46be
𝕂 = 1 / (10 - 𝕩 * 𝕪)

# ╔═╡ 97854f9f-9ba8-4e00-b993-98e38c78e07e
𝕗 = 1 + 𝕩^2

# ╔═╡ 10eb2d2c-e91e-496e-a565-7529193ae5c0
diff(𝕂, 𝕩)

# ╔═╡ 730ea0a5-47ad-4350-95eb-734ffa74647c
!(n::Int) = n == 0 ? 1 : n * !(n-1)

# ╔═╡ 4528c937-1aad-4c4a-8c70-4572f5e3acf1
𝕏 = sum(
	diff(diff(𝕂, 𝕩, p), 𝕪, r-p).subs(Dict(
		𝕩 => .5, 𝕪 => .5
	)) * (𝕩-.5)^p * (𝕪-.5)^(r-p) / !(p) / !(r-p)
	for p in 0:r
) |> expand |> simplify

# ╔═╡ ab58b8f5-1da5-4879-bd86-0b478e0fe6f3
𝕏.coeff(𝕩)

# ╔═╡ 516d7c40-63ca-4522-a0f9-d69008eabaac
prods = [
	(𝕩^p * 𝕪^q, p, q)
	for p in 0:r, q in 0:r if p+q ≤ 4
] |> reverse

# ╔═╡ 47304f39-654b-400a-9147-44d407bbbc45
begin
	buf = 𝕏
	coeffs = []
	for p in prods[1:end-1] .|> first
		push!(coeffs, buf.coeff(p) |> N)
		global buf -= buf.coeff(p) * p
	end
	push!(coeffs, buf |> N |> BigFloat)
	coeffs
end

# ╔═╡ e9cf6bff-fff7-47e0-9c2e-c2e07d94e78f
terms = OrderedDict(prods .=> coeffs)

# ╔═╡ 4f639682-365c-4875-84b7-caefca06069d
md"""
Каждый член суммы представим в виде $c_{p+q} * x^{p-1} * y^{q-1}$, тогда $c_{p+q} = a_{p+q} * p * q = (a_{p+q} * p^2 / q) * (a_{p+q} * q^2 / p)$. Функции $α(x), β(y)$ примут вид:

$\begin{cases}
	α(x) = a_{p+q} * q * x^{p-1} \\
	β(y) = p * y^{q-1}
\end{cases}$
"""

# ╔═╡ fcfc41de-8a30-4bd1-a395-b788d5b5dbdb
begin
	α, β = [], []
	term_keys = collect(keys(terms))
	for key in term_keys[1:end-1]
		push!(
			α, 
			terms[key] * (key[3]+1) * 𝕩^key[2] / sqrt((key[2]+1) * (key[3]+1))
		)
		push!(
			β, 
			(key[2]+1) * 𝕪^key[3] / sqrt((key[2]+1) * (key[3]+1))
		)
	end
	α[end] += terms[term_keys[end]] / β[end]
	α, β
end

# ╔═╡ 2ce03790-ed56-46d5-b6ce-5299af1d1882
# проверка на маленьковость ошибки
sum(
	α[i] * β[i]
	for i in 1:14
) - 𝕏

# ╔═╡ 5c517018-92d2-4ed0-a68a-82f3ba5f394e
md"""
## Построение системы и вычисление $Αᵢ$
"""

# ╔═╡ 333de524-2fad-4615-8163-58f8c200c45e
Β = -λ * [
	integrate(
		α[i] * β[j].subs(Dict(𝕪=>𝕩)),
		(𝕩, A, B)
	) |> N |> Float64
	for i in 1:length(α), j in 1:length(β)
] + I

# ╔═╡ 48bba8a2-1d1f-4ea0-8023-97574cf77b8b
F = [
	integrate(
		𝕗 * β[j].subs(Dict(𝕪=>𝕩)),
		(𝕩, A, B)
	) |> N
	for j in 1:length(β)
]

# ╔═╡ 54daf8c2-d927-45b2-a0f5-71538f157a17
Α = Β \ F

# ╔═╡ 60e8bb87-6557-46ee-8e42-f846462c27bc
𝕗 + λ * sum(Α[i] * α[i] for i in 1:14)

# ╔═╡ 9051c8d1-d27b-440f-8263-b14acac76652
φ(x::Number) = (𝕗 + λ * sum(Α[i] * α[i] for i in 1:14)).subs(Dict(𝕩 => x)) |> N

# ╔═╡ ebe5875b-91d4-4ca8-835e-3d1e057d1156
φ.(0:.01:1) |> plot

# ╔═╡ Cell order:
# ╟─e979dd52-76fa-11ed-3ecd-d1f9941a18d1
# ╟─15637745-5a70-4370-9058-553550a4c463
# ╠═6c9f99a5-59c8-4fed-a10c-1ea6f49aa7ee
# ╠═b58ac615-6f2b-49a2-9a41-cdb36406f165
# ╠═c93395e8-b244-4a10-9a35-c3797426e3e0
# ╠═f0f89c25-0424-4115-b597-645138579c19
# ╟─5d9526a8-3194-4e1b-a76b-aa8cda69058c
# ╠═367911ea-f525-4b4a-bc3c-6a9fa3e41f56
# ╠═80516d38-acc0-408f-8bd3-bd97dbc532c2
# ╠═292cf495-f860-4615-907d-7b06b5de46be
# ╠═97854f9f-9ba8-4e00-b993-98e38c78e07e
# ╠═10eb2d2c-e91e-496e-a565-7529193ae5c0
# ╠═730ea0a5-47ad-4350-95eb-734ffa74647c
# ╠═4528c937-1aad-4c4a-8c70-4572f5e3acf1
# ╠═ab58b8f5-1da5-4879-bd86-0b478e0fe6f3
# ╠═516d7c40-63ca-4522-a0f9-d69008eabaac
# ╠═47304f39-654b-400a-9147-44d407bbbc45
# ╠═e9cf6bff-fff7-47e0-9c2e-c2e07d94e78f
# ╟─4f639682-365c-4875-84b7-caefca06069d
# ╠═fcfc41de-8a30-4bd1-a395-b788d5b5dbdb
# ╠═2ce03790-ed56-46d5-b6ce-5299af1d1882
# ╟─5c517018-92d2-4ed0-a68a-82f3ba5f394e
# ╠═333de524-2fad-4615-8163-58f8c200c45e
# ╠═48bba8a2-1d1f-4ea0-8023-97574cf77b8b
# ╠═54daf8c2-d927-45b2-a0f5-71538f157a17
# ╠═60e8bb87-6557-46ee-8e42-f846462c27bc
# ╠═9051c8d1-d27b-440f-8263-b14acac76652
# ╠═ebe5875b-91d4-4ca8-835e-3d1e057d1156
