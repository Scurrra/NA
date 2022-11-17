### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 158d720c-60f8-11ed-303b-8d95a09a1be6
begin
	import Pkg; 	Pkg.activate()

	using LinearAlgebra
end

# ╔═╡ f0ea5fd5-030f-47df-8b5e-0ffa2e9ff402
md"""
# Метод сеток

$\frac{∂ u}{∂ t} = \frac{∂^2 u}{∂ x^2} + 0.5*(x^2 - 2t)$

$0≤x≤1, 0≤t≤0.02$

$u(0,t)=0, u(1,t)=0.5t$

$u(x,0)=0$
"""

# ╔═╡ 5d737f47-e338-453d-8e4c-51f01fa7ea7f
a(x,t) = 1

# ╔═╡ 2b9cd96d-d650-4035-9e3f-22562a9604ce
φ(x,t) = .5 * (x^2 - 2*t)

# ╔═╡ 71b24ccd-edee-4fb6-b8b7-891af3ea1ca0
γ₀(t) = 0

# ╔═╡ 6d1873cf-375b-438f-8cbb-51562824a671
γ₁(t) = .5*t

# ╔═╡ 27d8cfdf-efdf-44b4-99f7-ee8ff56c7f87
ψ(t) = 0

# ╔═╡ ad2c17ef-99fc-44fc-834f-8f9c2c3d065d
h = .1

# ╔═╡ 187f3f83-af5c-42a7-a2e4-951a7c7dc246
x = 0:h:1 |> collect

# ╔═╡ f29e75d7-3f36-4d7a-ab27-438080503d9e
md"""
## По явной разрастной схеме
"""

# ╔═╡ c42b0030-2e1d-41f8-8b24-dcfe58b9de9a
τ_explicit = .005

# ╔═╡ b6d4e09a-df97-4976-84a3-a71b61a9451e
t_explicit = 0:τ_explicit:0.02 |> collect

# ╔═╡ 3db660c8-2899-48b3-a9ae-07a255a11b3c
u_explicit = zeros((length(t_explicit), length(x)))

# ╔═╡ b8374fb8-669b-43ba-8195-2a5d8c120008
for n in 2:length(t_explicit)-1
	u_explicit[n, 1], u_explicit[n, length(x)] = γ₀(t_explicit[n]), γ₁(t_explicit[n])
end

# ╔═╡ 1f486f33-b995-411b-994b-97ff7350b43c
for m in 2:length(x)-1
	u_explicit[1, m] = ψ(x[m])
end

# ╔═╡ 1aa60332-77ef-4dd7-ba24-68dfdeace261
u_explicit

# ╔═╡ 3fa3be0e-663c-42c2-92c4-4c8c26aa7e33
s_explicit = τ_explicit / h^2

# ╔═╡ 82c427af-14b6-4caa-abd3-03ff09644304
for m in 2:length(x)-1, n in 2:length(t_explicit)-1
	u_explicit[n+1, m] = s_explicit * u_explicit[n, m+1] + 
		(1 - 2 * s_explicit) * u_explicit[n, m] + 
		s_explicit * u_explicit[n, m-1] + τ_explicit * φ(x[m], t_explicit[n])
end

# ╔═╡ 5b6629e1-0843-4bbb-bb83-418139f704a8
u_explicit

# ╔═╡ 34e40968-9f47-45e1-b9ae-d5c5e169823c
md"""
## По неявной разрастной схеме
"""

# ╔═╡ 4bc77f03-9456-48de-8410-dda1e81f662c
τ_implicit = .02

# ╔═╡ 35b194f1-3235-43a0-bd62-1b2fcda11ad7
t_implicit = 0:τ_implicit:0.02 |> collect

# ╔═╡ 017c4f0f-a587-42f7-8ac7-7ffe701833b2
u_implicit = zeros((length(t_implicit), length(x)))

# ╔═╡ badcd11f-47a2-4f7d-8c10-82494da9435f
for n in 1:length(t_implicit)
	u_implicit[n, 1], u_implicit[n, length(x)] = γ₀(t_implicit[n]), γ₁(t_implicit[n])
end

# ╔═╡ 2fbaa7ab-fffa-4753-a4e6-de6095d21a22
for m in 2:length(x)-1
	u_implicit[1, m] = ψ(x[m])
end

# ╔═╡ c38390d3-1caa-4980-b855-a42dddecd596
u_implicit

# ╔═╡ d08c6b7d-5bb4-4661-9270-f75d413e9ebd
s_implicit = τ_implicit / h^2

# ╔═╡ d2bd5c63-679b-4b94-8208-19aa52b2b5b9
A = diagm(
	0=>fill(s_implicit, length(x)-2),
	1=>fill(-(1+2*s_implicit), length(x)-2),
	2=>fill(s_implicit, length(x)-2)
)[1:end-2, :]

# ╔═╡ 26f7d797-4aa1-43f8-8ae5-98e995162797
b = -u_implicit[1, 2:end-1] - 
	[τ_implicit * φ(x[m], t_implicit[1]) for m in 2:length(x)-1]

# ╔═╡ 25af3057-aa37-4e0c-ab2f-23a1418b25b5
u_implicit[2, :] = A \ b

# ╔═╡ b0e52b18-dc98-4e74-a290-377d84cf783b
u_implicit

# ╔═╡ Cell order:
# ╟─158d720c-60f8-11ed-303b-8d95a09a1be6
# ╠═f0ea5fd5-030f-47df-8b5e-0ffa2e9ff402
# ╠═5d737f47-e338-453d-8e4c-51f01fa7ea7f
# ╠═2b9cd96d-d650-4035-9e3f-22562a9604ce
# ╠═71b24ccd-edee-4fb6-b8b7-891af3ea1ca0
# ╠═6d1873cf-375b-438f-8cbb-51562824a671
# ╠═27d8cfdf-efdf-44b4-99f7-ee8ff56c7f87
# ╠═ad2c17ef-99fc-44fc-834f-8f9c2c3d065d
# ╠═187f3f83-af5c-42a7-a2e4-951a7c7dc246
# ╟─f29e75d7-3f36-4d7a-ab27-438080503d9e
# ╠═c42b0030-2e1d-41f8-8b24-dcfe58b9de9a
# ╠═b6d4e09a-df97-4976-84a3-a71b61a9451e
# ╠═3db660c8-2899-48b3-a9ae-07a255a11b3c
# ╠═b8374fb8-669b-43ba-8195-2a5d8c120008
# ╠═1f486f33-b995-411b-994b-97ff7350b43c
# ╠═1aa60332-77ef-4dd7-ba24-68dfdeace261
# ╠═3fa3be0e-663c-42c2-92c4-4c8c26aa7e33
# ╠═82c427af-14b6-4caa-abd3-03ff09644304
# ╠═5b6629e1-0843-4bbb-bb83-418139f704a8
# ╟─34e40968-9f47-45e1-b9ae-d5c5e169823c
# ╠═4bc77f03-9456-48de-8410-dda1e81f662c
# ╠═35b194f1-3235-43a0-bd62-1b2fcda11ad7
# ╠═017c4f0f-a587-42f7-8ac7-7ffe701833b2
# ╠═badcd11f-47a2-4f7d-8c10-82494da9435f
# ╠═2fbaa7ab-fffa-4753-a4e6-de6095d21a22
# ╠═c38390d3-1caa-4980-b855-a42dddecd596
# ╠═d08c6b7d-5bb4-4661-9270-f75d413e9ebd
# ╠═d2bd5c63-679b-4b94-8208-19aa52b2b5b9
# ╠═26f7d797-4aa1-43f8-8ae5-98e995162797
# ╠═25af3057-aa37-4e0c-ab2f-23a1418b25b5
# ╠═b0e52b18-dc98-4e74-a290-377d84cf783b
