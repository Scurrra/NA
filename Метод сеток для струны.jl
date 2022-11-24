### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ e82ff0d4-6677-11ed-3106-1fa3a3bcd132
begin
	import Pkg; 	Pkg.activate()

	using DataFrames

	using LinearAlgebra
end

# ╔═╡ 44ed6842-c347-4556-b1fc-31e451454fe5
md"""
# Метод сеток для струны

$\frac{∂^2 u}{∂ t^2} = \frac{∂^2 u}{∂ x^2} - \frac{1.28}{(x + 0.6t + 1)^3}$

$0≤x≤1, 0≤t≤1$

$u(0,t)=\frac{1}{0.6t+1}, u(1,t)=\frac{1}{0.6t+2}$

$u(x,0)=\frac{1}{x+1}, \frac{∂ u(x, 0)}{∂ t}=\frac{-0.6}{(x+1)^2}$
"""

# ╔═╡ 75b89c36-904f-455e-bedf-ce7706129e02
a(x, t) = 1

# ╔═╡ fd47ef09-d9b3-4859-a313-6cd37669feb1
φ(x, t) = -1.28 / (x + 0.6t + 1)^3

# ╔═╡ 9be24890-9c3a-4735-bdb5-cdfb2570354c
γ₀(t) = 1 / (0.6t + 1)

# ╔═╡ be1002cb-e8fa-49a7-ac6b-7831ee159f5a
γ₁(t) = 1 / (0.6t + 2)

# ╔═╡ b38a029e-4fbb-4278-a38a-ba921218d4dc
α(x) = 1 / (x + 1)

# ╔═╡ f8b114a5-79f1-4dc0-9c63-89429f3f1823
α₂(x) = 2 / (x + 1)^3

# ╔═╡ a75690a1-940e-44fe-9911-3f0538c8ca95
β(x) = -0.6 / (x + 1)^2

# ╔═╡ 04909472-89e6-441b-807a-d10327b435b7
md"""
## Решение
"""

# ╔═╡ 7b80c814-a0d2-42d0-95bd-bab8fd850524
h = .1

# ╔═╡ 13e3f436-bc18-4648-ab3e-5e284b620ff1
τ = 0.05

# ╔═╡ cc06079a-b137-44f9-b9de-10591dcc53b8
X = collect(0:h:1)

# ╔═╡ f1b4cfe0-5746-4710-96b6-e5b3c60508a3
T = collect(0:τ:1)

# ╔═╡ 7de1c12a-d741-4300-af75-eb335cdbe7ca
M, N = length(X), length(T)

# ╔═╡ a740d295-a981-4918-92a0-e859cad1b516
U = zeros(N, M)

# ╔═╡ d83ad64a-0db3-439a-8afa-f37fbe3c87c8
begin
	U[:, 1] = γ₀.(T)
	U[:, M] = γ₁.(T)
	
	U
end

# ╔═╡ ad756cc2-f379-4fab-889c-e80632c023f4
begin
	U[1, 2:(M-1)] = α.(X[2:(M-1)])
	U
end

# ╔═╡ 5ea9bab3-a7ce-4077-8688-c053954260bf
begin
	U[2, 2:(M-1)] = α.(X[2:(M-1)]) + 
		τ*β.(X[2:(M-1)])
		(τ^2)/2 * (a.(X[2:(M-1)], 0) .* α₂.(X[2:(M-1)]) + φ.(X[2:(M-1)], 0))
	U
end

# ╔═╡ 6de10daf-ab98-4560-9374-35acb2de6153
s = τ^2 / h^2

# ╔═╡ 4c09c3ec-dab6-41f1-aa49-9f90e02eb518
for m in 2:(M-1), n in 2:(N-1)
	U[n+1, m] = s * U[n, m+1] + 2*(1-s)*U[n, m] + s*U[n, m-1] + τ^2 * φ(X[m], T[n])
end

# ╔═╡ 88826f91-2725-4a65-86bf-e95e4e5029f8
U

# ╔═╡ d4ef29a8-76f0-4e33-988b-b181e6a6b8a8
DataFrame(
	[T U],
	[:t; Symbol.("x=" .* string.(X))]
)

# ╔═╡ Cell order:
# ╟─e82ff0d4-6677-11ed-3106-1fa3a3bcd132
# ╠═44ed6842-c347-4556-b1fc-31e451454fe5
# ╠═75b89c36-904f-455e-bedf-ce7706129e02
# ╠═fd47ef09-d9b3-4859-a313-6cd37669feb1
# ╠═9be24890-9c3a-4735-bdb5-cdfb2570354c
# ╠═be1002cb-e8fa-49a7-ac6b-7831ee159f5a
# ╠═b38a029e-4fbb-4278-a38a-ba921218d4dc
# ╠═f8b114a5-79f1-4dc0-9c63-89429f3f1823
# ╠═a75690a1-940e-44fe-9911-3f0538c8ca95
# ╟─04909472-89e6-441b-807a-d10327b435b7
# ╟─7b80c814-a0d2-42d0-95bd-bab8fd850524
# ╟─13e3f436-bc18-4648-ab3e-5e284b620ff1
# ╟─cc06079a-b137-44f9-b9de-10591dcc53b8
# ╟─f1b4cfe0-5746-4710-96b6-e5b3c60508a3
# ╠═7de1c12a-d741-4300-af75-eb335cdbe7ca
# ╠═a740d295-a981-4918-92a0-e859cad1b516
# ╠═d83ad64a-0db3-439a-8afa-f37fbe3c87c8
# ╠═ad756cc2-f379-4fab-889c-e80632c023f4
# ╠═5ea9bab3-a7ce-4077-8688-c053954260bf
# ╠═6de10daf-ab98-4560-9374-35acb2de6153
# ╠═4c09c3ec-dab6-41f1-aa49-9f90e02eb518
# ╠═88826f91-2725-4a65-86bf-e95e4e5029f8
# ╟─d4ef29a8-76f0-4e33-988b-b181e6a6b8a8
