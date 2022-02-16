### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 0cb92a82-8e93-11ec-091b-f1abb82c2225
begin
	import Pkg;
	Pkg.activate()

	using LinearAlgebra

	using PlutoUI
	TableOfContents()
end

# ╔═╡ 32f02311-d88b-47ef-a399-59752ce037a5
md"""
# Лабораторная работа #1
# Решение СЛАУ методом Гаусса
"""

# ╔═╡ 2fd3c08f-0768-4533-90be-d0f457f77bdf
begin
	∞(A::Matrix{<:Number}) = maximum( sum(abs, A, dims=2) )
	∞(v::Vector{<:Number}) = maximum( abs, v )
end

# ╔═╡ 90464d1e-75aa-423c-9f29-dce94cf84003
begin
	|₁(A::Matrix{<:Number}) = maximum( sum(abs, A, dims=1) )
	|₁(v::Vector{<:Number}) = sum( abs, v )
end

# ╔═╡ 8d82ca3e-b058-49b7-8398-064de7da35d3
begin
	|₂(A::Matrix{<:Number}) = maximum( eigen(A' * A).values )
	|₂(v::Vector{<:Number}) = sum( abs2, v ) |> sqrt
end

# ╔═╡ f3572ec8-4ec9-4b87-b39b-34f3f652c8a0
md"""
## Задание 1.
"""

# ╔═╡ f1857dc6-6a51-4bed-b2a0-749c325b7794
D = [
	6.22 1.42 -1.72 1.91
	1.44 5.33 1.11 -1.82
	-1.72 1.11 5.24 1.42
	1.91 -1.82 1.42 6.55
]

# ╔═╡ 6bbb66c7-5458-4482-93e3-06293bf7dbb7
A = D + I

# ╔═╡ 53a63945-6809-4fc4-8987-21bd7b7c5e90
f̄ = [7.53, 6.06, 8.05, 8.06]

# ╔═╡ d052cb09-e039-4c65-96ad-e14c32a0fbbc


# ╔═╡ 9ec40330-5abf-432a-98fe-a40ac88f4da1


# ╔═╡ 7a4ea9df-c08b-4ef1-8f1e-df476d8698a3


# ╔═╡ ba30acac-a40f-4034-b16d-8d4292477d1e


# ╔═╡ 98618ada-6cc9-4436-802b-b66acda14603


# ╔═╡ 0e3dde33-49e7-4b4d-ba0d-ab68f1b4e382
md"""
## Задание 4. Матричные нормы
### Кубическая норма
"""

# ╔═╡ c3145d32-6325-477b-bbcd-0d7809b837ca
∞(A)

# ╔═╡ ae51d816-32f5-4bd6-9805-a53b556ac05c
md"""
### Октоэдрическая норма
"""

# ╔═╡ 2d89109b-dd51-4a22-acab-ed8bea61880c
|₁(A)

# ╔═╡ 73340101-8c51-4f19-a200-58c99002138f
md"""
### Сферическая норма
"""

# ╔═╡ df3c1167-ed28-4b77-8ac7-0f8ebddee9ae
|₂(A)

# ╔═╡ Cell order:
# ╟─0cb92a82-8e93-11ec-091b-f1abb82c2225
# ╟─32f02311-d88b-47ef-a399-59752ce037a5
# ╠═2fd3c08f-0768-4533-90be-d0f457f77bdf
# ╠═90464d1e-75aa-423c-9f29-dce94cf84003
# ╠═8d82ca3e-b058-49b7-8398-064de7da35d3
# ╟─f3572ec8-4ec9-4b87-b39b-34f3f652c8a0
# ╠═f1857dc6-6a51-4bed-b2a0-749c325b7794
# ╠═6bbb66c7-5458-4482-93e3-06293bf7dbb7
# ╠═53a63945-6809-4fc4-8987-21bd7b7c5e90
# ╠═d052cb09-e039-4c65-96ad-e14c32a0fbbc
# ╠═9ec40330-5abf-432a-98fe-a40ac88f4da1
# ╠═7a4ea9df-c08b-4ef1-8f1e-df476d8698a3
# ╠═ba30acac-a40f-4034-b16d-8d4292477d1e
# ╠═98618ada-6cc9-4436-802b-b66acda14603
# ╟─0e3dde33-49e7-4b4d-ba0d-ab68f1b4e382
# ╠═c3145d32-6325-477b-bbcd-0d7809b837ca
# ╟─ae51d816-32f5-4bd6-9805-a53b556ac05c
# ╠═2d89109b-dd51-4a22-acab-ed8bea61880c
# ╟─73340101-8c51-4f19-a200-58c99002138f
# ╠═df3c1167-ed28-4b77-8ac7-0f8ebddee9ae
