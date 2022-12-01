### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 8e28362a-717e-11ed-1120-6918292191a0
begin
	import Pkg; 	Pkg.activate()

	using DataFrames
end

# ╔═╡ 3f2b10bd-7da1-4dcf-8072-70e805f4317e
md"""
# Метод "кривых" сеток для Пуассона

$\frac{∂^2 u(x, y)}{∂ x^2} + \frac{∂^2 u(x, y)}{∂ y^2} = 1.91$

$x^2 + y^2 = 1$

$u_Γ(x,y) = 0.7x^2 + 0.3y^2 + 0.6$
"""

# ╔═╡ 9be2deea-747f-4807-8459-432b52d07909
h = 0.2

# ╔═╡ 66827d38-9627-4de2-9f13-3509fe1adbbd
UΓ((x, y)) = 0.7*x^2 + 0.3*y^2 + 0.6

# ╔═╡ c2fc543f-ac8a-4bed-be9a-8814087e07e0
md"""
## Область значения функции
"""

# ╔═╡ 28d011a2-7059-45f8-b904-50b7f0855464
begin
	X = -1:h:1
	Y = -1:h:1

	X, Y
end

# ╔═╡ b86ee75b-80c8-48aa-86dd-3e5cd0fb831f
Γ = [
	x^2 + y^2 ≤ 1 ? (x, y) : missing
	for x in X, y in Y
]

# ╔═╡ 6b9a05af-a745-45b7-b180-c20a2d6ef755
Γ_ind = collect(CartesianIndices(Γ))[.!ismissing.(Γ)]

# ╔═╡ 0c589cfa-83c6-4433-94a2-90d35209c30a
begin
	Β = [
		isapprox(x^2 + y^2, 1; atol=0.2) ? (x, y) : missing
		for x in X, y in Y
	]

	Β_ind = collect(CartesianIndices(Β))[.!ismissing.(Γ) .&& .!ismissing.(Β)]
	Β[.!(.!ismissing.(Γ) .&& .!ismissing.(Β))] .= missing
	Β
end

# ╔═╡ 20b6d048-b454-4cec-8877-05a33e76796b
md"""
## Индексы точек внутренней области
"""

# ╔═╡ 1e231d97-2f5f-42cc-9e95-3138265cac69
Β_ind

# ╔═╡ 9971b77c-3b7b-4f99-9fc9-ace5edc4998b
I_ind = Γ_ind[.!(Γ_ind .∈ (Β_ind,))]

# ╔═╡ 63ea5ed7-2397-47ff-ae46-140e574c8661
isdisjoint(Β_ind, I_ind)

# ╔═╡ 011bc374-ea63-4a6e-b210-7cdaf1efab55
md"""
## Решение
"""

# ╔═╡ 80ce6b26-f425-4df9-9d5e-29a9615cae9e
begin
	U = fill(Inf, size(Γ))
	U[Β_ind] = UΓ.(Β[Β_ind])
	
	U
end

# ╔═╡ 0fb42765-fd3f-41eb-b082-888f05e9b724
U[I_ind] .= 0 

# ╔═╡ cf1e66c3-4f2d-485f-8df3-810b8eb7d653
U

# ╔═╡ 122e4823-4843-42fc-8012-531a300f880a
collect(CartesianIndices(U))[I_ind]

# ╔═╡ 58d83b6c-bd0f-4e55-b4d3-599afded82eb
begin
	M = hcat([
		(collect(LinearIndices(U)).-6)[I_ind],
		(collect(LinearIndices(U)).-1)[I_ind],
		collect(LinearIndices(U))[I_ind],
		(collect(LinearIndices(U)).+1)[I_ind],
		(collect(LinearIndices(U)).+6)[I_ind]
	]...)
	
	mapping = Dict(
		collect(CartesianIndices(U))[I_ind] .=> [
		M[i, :]
		for i in 1:size(M, 1)
	])
end

# ╔═╡ 63410f59-548c-4d25-a540-c51207f2d72e
Φ = fill(1.91 * h^2, 121)

# ╔═╡ f7af4ec0-5eda-4f58-92a9-c546f8af120f
begin
	A = zeros(121, 121)
	
	for ind in I_ind
		A[mapping[ind][3], mapping[ind]] = [1, 1, -4, 1, 1]
	end
	
	A
end

# ╔═╡ 04ee94a1-1a5b-4eb5-bae3-7bbaf2d53a20
# индексы точек, где будем решать систему
S_ind = unique(sort(vcat(last.(collect(mapping))...)))

# ╔═╡ 252b809b-aa0a-4add-84e3-1dcebd88b976
begin
	u = zeros(11, 11)
	u[1:121] = A[S_ind, :] \ Φ[S_ind]
end

# ╔═╡ 3eb7bf8e-6fea-4206-be13-4ee7dd7ed710
U[I_ind] = u[I_ind]

# ╔═╡ 90ff9884-e93b-41e2-a4af-d1419f9da9a4
DataFrame(
	[Y U],
	["Y"; "X=".*string.(X)]
)

# ╔═╡ Cell order:
# ╟─8e28362a-717e-11ed-1120-6918292191a0
# ╟─3f2b10bd-7da1-4dcf-8072-70e805f4317e
# ╠═9be2deea-747f-4807-8459-432b52d07909
# ╠═66827d38-9627-4de2-9f13-3509fe1adbbd
# ╟─c2fc543f-ac8a-4bed-be9a-8814087e07e0
# ╠═28d011a2-7059-45f8-b904-50b7f0855464
# ╠═b86ee75b-80c8-48aa-86dd-3e5cd0fb831f
# ╠═6b9a05af-a745-45b7-b180-c20a2d6ef755
# ╠═0c589cfa-83c6-4433-94a2-90d35209c30a
# ╟─20b6d048-b454-4cec-8877-05a33e76796b
# ╠═1e231d97-2f5f-42cc-9e95-3138265cac69
# ╠═9971b77c-3b7b-4f99-9fc9-ace5edc4998b
# ╠═63ea5ed7-2397-47ff-ae46-140e574c8661
# ╟─011bc374-ea63-4a6e-b210-7cdaf1efab55
# ╠═80ce6b26-f425-4df9-9d5e-29a9615cae9e
# ╠═0fb42765-fd3f-41eb-b082-888f05e9b724
# ╠═cf1e66c3-4f2d-485f-8df3-810b8eb7d653
# ╠═122e4823-4843-42fc-8012-531a300f880a
# ╠═58d83b6c-bd0f-4e55-b4d3-599afded82eb
# ╠═63410f59-548c-4d25-a540-c51207f2d72e
# ╠═f7af4ec0-5eda-4f58-92a9-c546f8af120f
# ╠═04ee94a1-1a5b-4eb5-bae3-7bbaf2d53a20
# ╠═252b809b-aa0a-4add-84e3-1dcebd88b976
# ╠═3eb7bf8e-6fea-4206-be13-4ee7dd7ed710
# ╠═90ff9884-e93b-41e2-a4af-d1419f9da9a4
