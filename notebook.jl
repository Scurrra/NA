### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ e3ebbb44-daa3-11ec-36c2-a5272e97694d
begin
	import Pkg;
	Pkg.activate()

	using LinearAlgebra
	using DataFrames, PrettyTables
	using Polynomials

	using PlutoUI
	TableOfContents()
end

# ╔═╡ 002d16fa-9a2a-477d-bfd7-4324655183cd
A = [
	1 1.5 2.5 3.5
	1.5 1 2 1.6
	2.5 2 1 1.7
	3.5 1.6 1.7 1
]

# ╔═╡ a83615c0-41b8-4fc3-afda-badc15241a86
md"""
## M₃
"""

# ╔═╡ 50139705-6ebb-465f-b92e-3c9629be0cf4
begin
	M₃ = diagm(ones(4))
	M₃[3, :] = -A[4, :]
	M₃[3, 3] = 1
	M₃[3, :] = M₃[3, :] ./ A[4, 3]
	M₃
end

# ╔═╡ 1ad77896-48da-4c0a-ac44-021012102c92
begin
	M₃⁻¹ = diagm(ones(4))
	M₃⁻¹[3, :] = A[4, :]
	M₃⁻¹
end

# ╔═╡ 4238440a-72f3-4a10-abd0-33df63bd24a9
A¹ = M₃⁻¹ * A * M₃

# ╔═╡ 65b0b1f4-995c-4f4f-a174-76a01f1860f4
md"""
## M₂
"""

# ╔═╡ c2209805-5b74-480b-a21e-a541a416ac05
begin
	M₂ = diagm(ones(4))
	M₂[2, :] = -A¹[3, :]
	M₂[2, 2] = 1
	M₂[2, :] = M₂[2, :] ./ A¹[3, 2]
	M₂
end

# ╔═╡ 5355e704-8267-4818-b23b-bc64ea09469f
begin
	M₂⁻¹ = diagm(ones(4))
	M₂⁻¹[2, :] = A¹[3, :]
	M₂⁻¹
end

# ╔═╡ 9bf294dd-0852-42d5-89b7-60297810158b
A² = M₂⁻¹ * A¹ * M₂

# ╔═╡ 1b024fb4-0711-4819-afec-0f4469f1dc29
md"""
## M₁
"""

# ╔═╡ 61f32ca4-f0b3-433b-9703-7031be7d6ec8
begin
	M₁ = diagm(ones(4))
	M₁[1, :] = -A²[2, :]
	M₁[1, 1] = 1
	M₁[1, :] = M₁[1, :] ./ A²[2, 1]
	M₁
end

# ╔═╡ 81b369b4-1d6b-49e1-83e4-2187dae242c1
begin
	M₁⁻¹ = diagm(ones(4))
	M₁⁻¹[1, :] = A²[2, :]
	M₁⁻¹
end

# ╔═╡ 25f66de9-086f-45a5-9332-6893fb40de39
A³ = M₁⁻¹ * A² * M₁

# ╔═╡ b66f37b1-961f-4e9e-b060-60297f9b1054
md"""
## S и S⁻¹
"""

# ╔═╡ ae70a894-0930-41e7-93ab-732da8d98397
S = M₃ * M₂ * M₁ 

# ╔═╡ 8f681374-2069-438d-8e60-53a41b8fe113
S⁻¹ = M₁⁻¹ * M₂⁻¹ * M₃⁻¹

# ╔═╡ 7f5ba1a7-46a2-41ff-9788-f533bdb11ba1
S⁻¹ * A * S

# ╔═╡ a56d6cdf-690d-433d-ad8a-d3a69b346b58
eigen(A).values

# ╔═╡ 6df3c70d-9a61-4629-bb6f-eaa499fc0c05
eigen(S⁻¹ * A * S).values

# ╔═╡ 9c4927c3-39de-4a21-ae87-856b19065f05
md"""
## Polynom
"""

# ╔═╡ 3f31618d-03c6-481b-be95-1fd8aa16532f
Polynomial([1; -A³[1, :]] |> reverse)

# ╔═╡ 8f85569c-269f-4b40-aef5-0c475c6b33e3
λ = Polynomial([1; -A³[1, :]] |> reverse) |> roots

# ╔═╡ c4703c87-c92d-4f0e-a37e-d1bc61e54407
ȳ(i) = [λ[i]^3, λ[i]^2, λ[i], 1]

# ╔═╡ daf71c2a-22a0-48bc-9dcc-29e1fbfbfdfe
ȳ(1)

# ╔═╡ 0e3ed9f8-b076-411a-a00e-ec091ad1ab78
x̄(i) = M₃ * M₂ * M₁ * ȳ(i)

# ╔═╡ e826069b-d16e-4d25-9f8e-0699e63694c6
M₁ * ȳ(1)

# ╔═╡ 12de2b7e-8d06-4e26-9fb1-eb325bda4b98
M₂ * M₁ * ȳ(1)

# ╔═╡ 6bead8fc-a2d0-459f-8bae-5d86a8408fe2
M₃ * M₂ * M₁ * ȳ(1)

# ╔═╡ d9c818e7-f7d1-407b-9ce4-e3fb2bf90bb5
x̄(1)

# ╔═╡ f995d12e-07fb-4d5e-949b-e3bbc40dccb4
md"""
## Check
"""

# ╔═╡ 433da834-cd3f-445e-bd9b-d4892e0db014
S⁻¹ * A * S * ȳ(1), λ[1] * ȳ(1)

# ╔═╡ 4899aada-9b4c-4e57-8de3-38dda92fa7f0
A * x̄(1), λ[1] * x̄(1)

# ╔═╡ 25d1cc43-8c2e-4fd1-92fc-1e3054ec8902


# ╔═╡ e4d571d4-1c4d-4011-91b3-13003da57141


# ╔═╡ Cell order:
# ╟─e3ebbb44-daa3-11ec-36c2-a5272e97694d
# ╟─002d16fa-9a2a-477d-bfd7-4324655183cd
# ╟─a83615c0-41b8-4fc3-afda-badc15241a86
# ╠═50139705-6ebb-465f-b92e-3c9629be0cf4
# ╠═1ad77896-48da-4c0a-ac44-021012102c92
# ╠═4238440a-72f3-4a10-abd0-33df63bd24a9
# ╟─65b0b1f4-995c-4f4f-a174-76a01f1860f4
# ╠═c2209805-5b74-480b-a21e-a541a416ac05
# ╠═5355e704-8267-4818-b23b-bc64ea09469f
# ╠═9bf294dd-0852-42d5-89b7-60297810158b
# ╟─1b024fb4-0711-4819-afec-0f4469f1dc29
# ╠═61f32ca4-f0b3-433b-9703-7031be7d6ec8
# ╠═81b369b4-1d6b-49e1-83e4-2187dae242c1
# ╠═25f66de9-086f-45a5-9332-6893fb40de39
# ╟─b66f37b1-961f-4e9e-b060-60297f9b1054
# ╠═ae70a894-0930-41e7-93ab-732da8d98397
# ╠═8f681374-2069-438d-8e60-53a41b8fe113
# ╠═7f5ba1a7-46a2-41ff-9788-f533bdb11ba1
# ╠═a56d6cdf-690d-433d-ad8a-d3a69b346b58
# ╠═6df3c70d-9a61-4629-bb6f-eaa499fc0c05
# ╟─9c4927c3-39de-4a21-ae87-856b19065f05
# ╠═3f31618d-03c6-481b-be95-1fd8aa16532f
# ╠═8f85569c-269f-4b40-aef5-0c475c6b33e3
# ╠═c4703c87-c92d-4f0e-a37e-d1bc61e54407
# ╠═daf71c2a-22a0-48bc-9dcc-29e1fbfbfdfe
# ╠═0e3ed9f8-b076-411a-a00e-ec091ad1ab78
# ╠═e826069b-d16e-4d25-9f8e-0699e63694c6
# ╠═12de2b7e-8d06-4e26-9fb1-eb325bda4b98
# ╠═6bead8fc-a2d0-459f-8bae-5d86a8408fe2
# ╠═d9c818e7-f7d1-407b-9ce4-e3fb2bf90bb5
# ╟─f995d12e-07fb-4d5e-949b-e3bbc40dccb4
# ╠═433da834-cd3f-445e-bd9b-d4892e0db014
# ╠═4899aada-9b4c-4e57-8de3-38dda92fa7f0
# ╠═25d1cc43-8c2e-4fd1-92fc-1e3054ec8902
# ╠═e4d571d4-1c4d-4011-91b3-13003da57141
