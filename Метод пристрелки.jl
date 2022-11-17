### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 6fbb17ae-55f7-11ed-2ee0-a99b1181f1dd
begin
	import Pkg; 	Pkg.activate()

	using Plots
	using DataFrames
end

# ╔═╡ d45451f0-7966-4507-9205-5ba2041e19cb
md"""
# Метод пристрелки

$y'(x) - \frac{2}{x}y'(x) - \frac{4}{x^2 + 2}y(x) = 8$

$0.5 ≤ x ≤ 1$

$y'(0.5) = 0.5, y(1) + y'(1) = 1$
"""

# ╔═╡ 0de6a492-11c0-4bd5-b4f7-841cd11c12ef
α₀, α₁ = 0, 1

# ╔═╡ ef3ce032-9a54-475c-8372-8205a50b618f
β₀, β₁ = 1, 1

# ╔═╡ 00259276-e8a6-41e9-86be-beffa0b202c6
γ₀, γ₁ = 0.5, 1

# ╔═╡ f3bbaa5c-63e9-484c-951f-7bc98f202130
z(x, y) = [
	y[2],
	-2/x*y[2] + 4/(x^2+2)*y[1] - 8
] 

# ╔═╡ cd74934c-541c-41d7-af06-37a8a62a2003
h = .01

# ╔═╡ 3fbebdcf-1fe4-44ef-aa3a-a5b3fbee0c9d
begin
	x₁ = [0.5]
	z₁ = [[γ₀, -1]]
end

# ╔═╡ 0094f09e-e386-4af9-9de5-8efac5f98f87
for i in 1:round(Int, (1-.5)/h)
	k₁ = h * z(x₁[end], z₁[end])
	k₂ = h * z(h/2 + x₁[end], k₁/2 + z₁[end])
	k₃ = h * z(h/2 + x₁[end], k₂/2 + z₁[end])
	k₄ = h * z(h/2 + x₁[end], k₃/2 + z₁[end])

	push!(x₁, x₁[end]+h)
	push!(z₁, z₁[end]+(k₁ + 2k₂ + 2k₃ + k₄)/6)
end

# ╔═╡ 49fe9460-a880-418d-898a-ecf3fa024a8b
first.(z₁)

# ╔═╡ 6ffc909f-dc30-493c-816a-4053f4b999f1
begin
	x₂ = [.5]
	z₂ = [[γ₀, 1]]
end

# ╔═╡ cbdadc39-cacd-4652-bb9f-49136b72a7b2
for i in 1:round(Int, (1-.5)/h)
	k₁ = h * z(x₂[end], z₂[end])
	k₂ = h * z(h/2 + x₂[end], k₁/2 + z₂[end])
	k₃ = h * z(h/2 + x₂[end], k₂/2 + z₂[end])
	k₄ = h * z(h/2 + x₂[end], k₃/2 + z₂[end])

	push!(x₂, x₂[end]+h)
	push!(z₂, z₂[end]+(k₁ + 2k₂ + 2k₃ + k₄)/6)
end

# ╔═╡ 2e9460e8-6f08-44e5-96a6-2cb93fbe32fb
first.(z₂)

# ╔═╡ 275674fb-eed9-47ed-9f6e-78736ba70440
x₂

# ╔═╡ 92117ba9-c667-4912-b2fe-6589d0a5acc3
t₂ = 1 - (first(z₂[end]) * 2)/(first(z₂[end]) - first(z₁[end]))

# ╔═╡ e92a5c72-32aa-41f6-9b94-cd025e311864
begin
	x₃ = [.5]
	z₃ = [[γ₀, t₂]]
end

# ╔═╡ b9be50c2-b3d7-48e7-9172-3e41128d663a
for i in 1:round(Int, (1-.5)/h)
	k₁ = h * z(x₃[end], z₃[end])
	k₂ = h * z(h/2 + x₃[end], k₁/2 + z₃[end])
	k₃ = h * z(h/2 + x₃[end], k₂/2 + z₃[end])
	k₄ = h * z(h/2 + x₃[end], k₃/2 + z₃[end])

	push!(x₃, x₃[end]+h)
	push!(z₃, z₃[end]+(k₁ + 2k₂ + 2k₃ + k₄)/6)
end

# ╔═╡ c0a50b7f-5509-4528-b374-9eeee59d6df1
first.(z₃)[1:10:end]

# ╔═╡ 679f5de4-7f99-4285-8679-e7006e8804f4
plot(
	x₃,
	[
		first.(z₁),
		first.(z₂),
		first.(z₃),
	];
	label=["z₀" "z₁" "z₂"]
)

# ╔═╡ da505992-f8a9-4c70-909f-b962743d32d5
DataFrame(
	x=x₁[1:10:end],
	z₁=first.(z₁)[1:10:end],
	z₂=first.(z₂)[1:10:end],
	z₃=first.(z₃)[1:10:end]
)

# ╔═╡ Cell order:
# ╟─6fbb17ae-55f7-11ed-2ee0-a99b1181f1dd
# ╟─d45451f0-7966-4507-9205-5ba2041e19cb
# ╠═0de6a492-11c0-4bd5-b4f7-841cd11c12ef
# ╠═ef3ce032-9a54-475c-8372-8205a50b618f
# ╠═00259276-e8a6-41e9-86be-beffa0b202c6
# ╠═f3bbaa5c-63e9-484c-951f-7bc98f202130
# ╠═cd74934c-541c-41d7-af06-37a8a62a2003
# ╠═3fbebdcf-1fe4-44ef-aa3a-a5b3fbee0c9d
# ╠═0094f09e-e386-4af9-9de5-8efac5f98f87
# ╠═49fe9460-a880-418d-898a-ecf3fa024a8b
# ╠═6ffc909f-dc30-493c-816a-4053f4b999f1
# ╠═cbdadc39-cacd-4652-bb9f-49136b72a7b2
# ╠═2e9460e8-6f08-44e5-96a6-2cb93fbe32fb
# ╠═275674fb-eed9-47ed-9f6e-78736ba70440
# ╠═92117ba9-c667-4912-b2fe-6589d0a5acc3
# ╠═e92a5c72-32aa-41f6-9b94-cd025e311864
# ╠═b9be50c2-b3d7-48e7-9172-3e41128d663a
# ╠═c0a50b7f-5509-4528-b374-9eeee59d6df1
# ╟─679f5de4-7f99-4285-8679-e7006e8804f4
# ╠═da505992-f8a9-4c70-909f-b962743d32d5
