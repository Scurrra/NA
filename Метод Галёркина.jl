### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 4c0586b8-4577-11ed-1644-933e14cee4c8
begin
	import Pkg; 	Pkg.activate()

	using SymPy
end

# ╔═╡ dbe335ed-e652-437f-baf8-a962ad945684
md"""
# Метод Галёркина

$\begin{cases}
	y''(x) + \frac{0.9}{1.8 x + 1}y'(x) = 0\\
	y(0) = 2 \\
	y(1) = 2 \sqrt{2.8}
\end{cases}$

"""

# ╔═╡ dc297401-0be2-4a55-8a6b-c6b030f6b0c2
A, B = 0., 1.

# ╔═╡ bd070091-3d90-4b13-909f-c8796c47ab4b
α₀, α₁ = 1., 1.

# ╔═╡ 939ada0b-ca09-4ff7-8d3e-02ffb9c581e4
γ₀, γ₁ = 2., 2*sqrt(2.8)

# ╔═╡ 70345809-aa46-4019-8b70-5a7cb01463de
@syms a[0:10], x

# ╔═╡ 0933c652-63b8-45c3-ba11-02165e33ba05
p = 0.9 / (1.8x + 1)

# ╔═╡ ce78d978-ee8b-49cb-b3db-9f9744d03fc1
φ = Dict(
	i-1 => a[1] + sum(a[j+1]*x^j for j in 1:i)
	for i in 1:5
)

# ╔═╡ 349fbf01-e76c-40ba-96e3-fafe47a4416c
md"""
## φ₀
"""

# ╔═╡ 36abfde9-1ba0-4082-a0fb-2059d6fc01d1
coeffs₀ = solve([
	Eq(α₀*φ[0].subs(Dict(x=>A)), γ₀),
	Eq(α₁*φ[0].subs(Dict(x=>B)), γ₁)
], a[1:2])

# ╔═╡ b3d6efce-bd8a-4b53-b2fb-751566a4603f
φ[0] = φ[0].subs(coeffs₀)

# ╔═╡ e9a082e9-a65d-4ae9-bcf3-ad83455bf8fd
md"""
## φᵢ
"""

# ╔═╡ f0d270b5-1b3c-46d0-bb0f-538ab9042ac0
for i in 1:4
	φ[i] = φ[i].subs(
		solve([
			Eq(α₀*φ[i].subs(Dict(x=>A)), 0),
			Eq(α₁*φ[i].subs(Dict(x=>B)), 0)
		], a[1:(i+2)])
	).subs(Dict(
		a[1:(i+2)] .=> rand(i+2)
	))
end

# ╔═╡ c51dd093-b4db-4dc6-b055-db283a9f38fb
φ

# ╔═╡ 9addb720-8959-4352-96e6-cc78565707db
md"""
## C, D, Α
"""

# ╔═╡ 18dc9a47-2b0a-4ba0-8d20-897e0a59cd59
C = [
	integrate(
		(diff(φ[k], x, 2) + p*diff(φ[k], x))*φ[i],
		(x, A, B)
	) |> N
	for k in 1:4, i in 1:4
]

# ╔═╡ cbdf522a-3f27-422d-b3a0-c117635830f6
D = [
	integrate(
		-(diff(φ[0], x, 2) + p*diff(φ[0], x))*φ[i],
		(x, A, B)
	) |> N
	for i in 1:4
]

# ╔═╡ d9a09fe3-d287-46ba-9c62-129a51df5ed9
Α = C \ D

# ╔═╡ 5680e56b-d32b-4bfa-880b-1d04353cc8ab
md"""
## Yₙ
"""

# ╔═╡ 55d2fee6-352f-4e66-85a0-4db0f7684f95
Y(X::Number; n=3) = (φ[0] + sum(Α[i] * φ[i] for i in 1:n)).subs(Dict(x=>X))

# ╔═╡ 730961f8-0817-4e00-8907-8e1abe40c193
Y(0.5; n=3)

# ╔═╡ 430400d2-065d-48c5-be9a-ef5782c37083
Y(0.5; n=4)

# ╔═╡ 4a50b380-b841-4f4c-8d75-a4a6db9821b4
[
	(φ[0] + sum(Α[i] * φ[i] for i in 1:j))
	for j in 1:4
]

# ╔═╡ Cell order:
# ╟─4c0586b8-4577-11ed-1644-933e14cee4c8
# ╟─dbe335ed-e652-437f-baf8-a962ad945684
# ╠═dc297401-0be2-4a55-8a6b-c6b030f6b0c2
# ╠═bd070091-3d90-4b13-909f-c8796c47ab4b
# ╠═939ada0b-ca09-4ff7-8d3e-02ffb9c581e4
# ╠═70345809-aa46-4019-8b70-5a7cb01463de
# ╠═0933c652-63b8-45c3-ba11-02165e33ba05
# ╠═ce78d978-ee8b-49cb-b3db-9f9744d03fc1
# ╟─349fbf01-e76c-40ba-96e3-fafe47a4416c
# ╠═36abfde9-1ba0-4082-a0fb-2059d6fc01d1
# ╠═b3d6efce-bd8a-4b53-b2fb-751566a4603f
# ╟─e9a082e9-a65d-4ae9-bcf3-ad83455bf8fd
# ╠═f0d270b5-1b3c-46d0-bb0f-538ab9042ac0
# ╠═c51dd093-b4db-4dc6-b055-db283a9f38fb
# ╟─9addb720-8959-4352-96e6-cc78565707db
# ╠═18dc9a47-2b0a-4ba0-8d20-897e0a59cd59
# ╠═cbdf522a-3f27-422d-b3a0-c117635830f6
# ╠═d9a09fe3-d287-46ba-9c62-129a51df5ed9
# ╟─5680e56b-d32b-4bfa-880b-1d04353cc8ab
# ╠═55d2fee6-352f-4e66-85a0-4db0f7684f95
# ╠═730961f8-0817-4e00-8907-8e1abe40c193
# ╠═430400d2-065d-48c5-be9a-ef5782c37083
# ╠═4a50b380-b841-4f4c-8d75-a4a6db9821b4
