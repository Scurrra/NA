### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 6509dd74-6bf8-11ed-12f4-47b6df1891db
begin
	import Pkg; 	Pkg.activate()

	using DataFrames
end

# ╔═╡ 7792829f-c65b-41b9-a241-e41ab2fc7fa4
md"""
# Метод сеток для Пуассона

$\frac{∂^2 u(x, y)}{∂ x^2} + \frac{∂^2 u(x, y)}{∂ y^2} = \frac{2}{1+y} \left( 1 + \left( \frac{x + 0.2}{1+y} \right)^2 \right)$

$0≤x, y≤1$

$u(x,0)=(x+0.2)^2, u(x,1)=\frac{(x+0.2)^2}{2}$

$u(0,y)=\frac{0.2^2}{1+y}, u(1,y)=\frac{1.2^2}{1+y}$
"""

# ╔═╡ dc381d88-07eb-4106-8f62-88bf84f9a784
φ(x, y) = 2 / (1+y) * (1 + ((x+0.2)/(1+y))^2)

# ╔═╡ c2e4214d-9ede-482c-b2e3-382258a8be3d
h = .2

# ╔═╡ 202975ee-3e27-41d0-8095-61718f7e0ca1
X = 0:h:1

# ╔═╡ 1829b2a2-7107-4a43-9410-c85150cfbd04
Y = 0:h:1

# ╔═╡ 3499385a-84b5-4e08-b032-1fe317e56534
α₀(x) = (x+0.2)^2

# ╔═╡ fd9243da-2cfd-460a-8e9f-1ca11d9d4d81
α₁(x) = 0.5 * (x+0.2)^2

# ╔═╡ bf6e189a-3514-4a24-a59e-42c6e68b0eeb
β₀(y) = 0.04 / (1+y)

# ╔═╡ ca1e4c1b-2b92-485b-b393-bedafa42c1cb
β₁(y) = 1.44 / (1+y)

# ╔═╡ 04e40916-4e98-4b46-b60f-6ec433725788
U = Matrix(undef, 6, 6)

# ╔═╡ 2ae96f39-0976-4972-898d-03e5c1063e16
begin
	U[1, :] = α₀.(X)
	U[end, :] = α₁.(X)
	U
end

# ╔═╡ 815bee99-432a-4e45-92db-81f1c4255b28
begin
	U[:, 1] = β₀.(Y)
	U[:, end] = β₁.(Y)
	U
end

# ╔═╡ 1ea7e59a-5a96-4a1c-afc1-98d1da53b290
σ = 1

# ╔═╡ cb44e0e8-9b41-4423-a4fc-6700aaab39e0
md"""
## Подбираем отображение для построения системы
"""

# ╔═╡ 6c85fc4f-0e19-4571-845d-6d2e7c866d39
mapping = Dict(
	collect(CartesianIndices(U))[2:end-1, 2:end-1] 
		.=> 
	collect(CartesianIndices((4, 4)))
)

# ╔═╡ 3544bf54-48db-45e2-8f4a-2428aedc3715
M = cat([
	(collect(LinearIndices(U)).-6)[2:end-1, 2:end-1],
	(collect(LinearIndices(U)).-1)[2:end-1, 2:end-1],
	collect(LinearIndices(U))[2:end-1, 2:end-1],
	(collect(LinearIndices(U)).+1)[2:end-1, 2:end-1],
	(collect(LinearIndices(U)).+6)[2:end-1, 2:end-1]
]..., dims=3)

# ╔═╡ 1e1f1c5c-da8b-48ff-bf92-cdc6b8b4e7bc
M[mapping[CartesianIndex((2, 2))], :]

# ╔═╡ 954b05cd-6132-4da3-8964-5aff75a9a027
md"""
## Решаем систему
"""

# ╔═╡ a2652a7e-e9b3-4e9c-b16f-33a82bdba037
Φ = vcat([
	φ(x, y)
	for x in X, y in Y
]...) * h^2

# ╔═╡ d602ac4f-4da9-4425-bf07-c2ddb8aa2402
begin
	A = zeros(36, 36)
	
	for key in collect(keys(mapping))
		A[M[mapping[key], :][3], M[mapping[key], :]] = [1, σ, -2*(1+σ), σ, 1]
	end

	A
end

# ╔═╡ 6589dccb-4314-427a-99da-c3cd8ad69251
[1, σ, -2*(1+σ), σ, 1]

# ╔═╡ b0547d53-ebb8-4d9b-9717-1dcee1a68d76
ind = vcat([
	M[mapping[key], :]
	for key in collect(keys(mapping))
]...) |> unique |> sort

# ╔═╡ c1866a9e-9030-46fb-bed1-d0b3d226c315
begin
	u = zeros(6, 6)
	u[1:36] = A[ind, :] \ Φ[ind]
end

# ╔═╡ 62fff2a2-88d9-42c0-9dfc-34cc3b0c4002
begin
	U[2:end-1, 2:end-1] = u[2:end-1, 2:end-1]
	U
end

# ╔═╡ Cell order:
# ╟─6509dd74-6bf8-11ed-12f4-47b6df1891db
# ╟─7792829f-c65b-41b9-a241-e41ab2fc7fa4
# ╠═dc381d88-07eb-4106-8f62-88bf84f9a784
# ╠═c2e4214d-9ede-482c-b2e3-382258a8be3d
# ╠═202975ee-3e27-41d0-8095-61718f7e0ca1
# ╠═1829b2a2-7107-4a43-9410-c85150cfbd04
# ╠═3499385a-84b5-4e08-b032-1fe317e56534
# ╠═fd9243da-2cfd-460a-8e9f-1ca11d9d4d81
# ╠═bf6e189a-3514-4a24-a59e-42c6e68b0eeb
# ╠═ca1e4c1b-2b92-485b-b393-bedafa42c1cb
# ╠═04e40916-4e98-4b46-b60f-6ec433725788
# ╠═2ae96f39-0976-4972-898d-03e5c1063e16
# ╠═815bee99-432a-4e45-92db-81f1c4255b28
# ╠═1ea7e59a-5a96-4a1c-afc1-98d1da53b290
# ╟─cb44e0e8-9b41-4423-a4fc-6700aaab39e0
# ╠═6c85fc4f-0e19-4571-845d-6d2e7c866d39
# ╠═3544bf54-48db-45e2-8f4a-2428aedc3715
# ╠═1e1f1c5c-da8b-48ff-bf92-cdc6b8b4e7bc
# ╟─954b05cd-6132-4da3-8964-5aff75a9a027
# ╠═a2652a7e-e9b3-4e9c-b16f-33a82bdba037
# ╠═d602ac4f-4da9-4425-bf07-c2ddb8aa2402
# ╠═6589dccb-4314-427a-99da-c3cd8ad69251
# ╠═b0547d53-ebb8-4d9b-9717-1dcee1a68d76
# ╠═c1866a9e-9030-46fb-bed1-d0b3d226c315
# ╠═62fff2a2-88d9-42c0-9dfc-34cc3b0c4002
