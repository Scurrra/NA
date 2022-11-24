### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ a309d30c-5b79-11ed-1d7b-07abf949f105
begin
	import Pkg; 	Pkg.activate()
	
	using DataFrames
end

# ╔═╡ e7905776-c7bd-4117-8f7d-2249441c99eb
md"""
# Метод сеток решения ОДУ

$y''(x) - 0.21^3 y(x) = - 0.21^2 x, 0 ≤ x ≤ 1$

$y(0) = 1, y(1) = ℯ^5 + 1$
"""

# ╔═╡ 8277e971-f860-4001-bd6e-1555c2f2a5a5
p(x) = 0

# ╔═╡ b0bae4d1-938f-4866-9a83-0403605270d0
q(x) = 0.21^3

# ╔═╡ b506d777-7ab5-4f51-9115-10b6b66822c1
f(x) = -0.21^2 * x

# ╔═╡ f6ab7a23-4aef-4bb9-bda8-2ee085c9c8cb
α₀, α₁ = 1, 1

# ╔═╡ 4e089f72-1260-4c6d-b85d-cf10c108d597
β₀, β₁ = 0, 0

# ╔═╡ 373357f3-8e83-436b-b61a-0a1bd2321bda
γ₀, γ₁ = 1, exp(0.21)+1

# ╔═╡ d05b5f73-d591-484a-8cc4-7a1c2b0f4a07
h = 0.1

# ╔═╡ fb9e3d29-3318-4196-84c0-be53bb6065f3
x = 0:h:1

# ╔═╡ f86057bb-798c-4e45-bf2e-4ffa17ad45ae
md"""
# a, b, c, t
"""

# ╔═╡ 6e64eb2f-6fa6-4174-ad84-e392a1d04149
a = 1 .+ h/2 * p.(x)

# ╔═╡ 5339ffd0-37e2-4a37-8edf-5fc024850e07
b = -(2 .+ h^2 * q.(x))

# ╔═╡ f18d2930-dcbd-4cc7-a0d1-72ecdb9b690f
c = 1 .- h/2 * p.(x)

# ╔═╡ ec1813da-c6be-41d4-8ecf-3b27a35f1378
t = h^2 * f.(x)

# ╔═╡ 3a42d48e-39f7-4f66-b10a-6717de801f2a
DataFrame(
	a=a,
	b=b,
	c=c,
	t=t
)

# ╔═╡ c4ab72b5-ec0f-4010-9440-47822f402f96
α₀₁, α₀₂ = α₀-β₀/h, -β₁/h

# ╔═╡ 6f615847-895e-41e9-ab8c-cfd218b3a830
β₀₁, β₀₂ = β₀/h, α₁+β₁/h

# ╔═╡ d1e28d81-3a67-4a78-bd95-7659f09694ac
γ₀₁, γ₀₂ = γ₀, γ₁

# ╔═╡ 64c11d0a-3994-47c7-bb5e-4d427f659a53
md"""
# X, Z
"""

# ╔═╡ a65b895a-08e9-438f-bcbd-8551006fdce5
all(abs.(b) .≥ abs.(c) + abs.(a))

# ╔═╡ 845215de-646a-4f76-943c-568adaf6fb8a
X, Z = [-β₀₁/α₀₁], [γ₀₁/α₀₁]

# ╔═╡ 3e952630-62b1-4352-86ff-80d0cd51a604
for i in 1:10
	push!(
		X,
		-a[i] / (b[i] + c[i]*X[i])
	)
	push!(
		Z,
		(t[i] - c[i]*Z[i]) / (b[i] + c[i]*X[i])
	)
end

# ╔═╡ 35b65ac1-23b3-4984-8fb1-5ad7aa1b5db8
DataFrame(
	X=X,
	Z=Z
)

# ╔═╡ 89841273-f162-415a-ab28-df84f01ea115
md"""
# Y
"""

# ╔═╡ 12e0e8e4-2f52-4daf-80e3-0838b8b754e4
Y = zeros(11)

# ╔═╡ 0c132ec4-a1a4-4e7c-824f-e677a396ecfe
Y[end] = ([
	1 -X[end]
	α₀₂ β₀₂
] \ [γ₀₂, Z[end]])[1]

# ╔═╡ 63085aae-e384-4223-b86c-d7775b83ad20
for i in 11:-1:2
	Y[i-1] = Y[i] * X[i] + Z[i]
end

# ╔═╡ 33a65907-37f0-4a17-905b-605d30fb558a
DataFrame(
	x=x,
	Y=Y
)

# ╔═╡ Cell order:
# ╠═a309d30c-5b79-11ed-1d7b-07abf949f105
# ╟─e7905776-c7bd-4117-8f7d-2249441c99eb
# ╠═8277e971-f860-4001-bd6e-1555c2f2a5a5
# ╠═b0bae4d1-938f-4866-9a83-0403605270d0
# ╠═b506d777-7ab5-4f51-9115-10b6b66822c1
# ╠═f6ab7a23-4aef-4bb9-bda8-2ee085c9c8cb
# ╠═4e089f72-1260-4c6d-b85d-cf10c108d597
# ╠═373357f3-8e83-436b-b61a-0a1bd2321bda
# ╠═d05b5f73-d591-484a-8cc4-7a1c2b0f4a07
# ╠═fb9e3d29-3318-4196-84c0-be53bb6065f3
# ╟─f86057bb-798c-4e45-bf2e-4ffa17ad45ae
# ╠═6e64eb2f-6fa6-4174-ad84-e392a1d04149
# ╠═5339ffd0-37e2-4a37-8edf-5fc024850e07
# ╠═f18d2930-dcbd-4cc7-a0d1-72ecdb9b690f
# ╠═ec1813da-c6be-41d4-8ecf-3b27a35f1378
# ╠═3a42d48e-39f7-4f66-b10a-6717de801f2a
# ╠═c4ab72b5-ec0f-4010-9440-47822f402f96
# ╠═6f615847-895e-41e9-ab8c-cfd218b3a830
# ╠═d1e28d81-3a67-4a78-bd95-7659f09694ac
# ╟─64c11d0a-3994-47c7-bb5e-4d427f659a53
# ╠═a65b895a-08e9-438f-bcbd-8551006fdce5
# ╠═845215de-646a-4f76-943c-568adaf6fb8a
# ╠═3e952630-62b1-4352-86ff-80d0cd51a604
# ╠═35b65ac1-23b3-4984-8fb1-5ad7aa1b5db8
# ╟─89841273-f162-415a-ab28-df84f01ea115
# ╠═12e0e8e4-2f52-4daf-80e3-0838b8b754e4
# ╠═0c132ec4-a1a4-4e7c-824f-e677a396ecfe
# ╠═63085aae-e384-4223-b86c-d7775b83ad20
# ╠═33a65907-37f0-4a17-905b-605d30fb558a
