### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 67d84dcc-5077-11ed-0ed2-c742f73ce6d5
begin
	import Pkg; 	Pkg.activate()

	using DataFrames
	using Plots
end

# ╔═╡ 9b13725c-0fdb-421f-bdc7-df4f02aadd9f
function euler(f, a, b, x₀, y₀, h)
	x = [x₀]
	y = [y₀]
	
	while x[end]+h < b
		push!(x, x[end] + h)
		push!(y, y[end] + h * f(x[end], y[end]))
	end

	return x, y
end

# ╔═╡ 9105fc1a-9f8a-4aa3-b242-106715871524
md"""
# Метод дифференциальной прогонки

31 m = 10

$y'' + (x+1)y' - \frac{20}{(x+1)^2}y = 5 (x+1)^5$

$0 ≤ x ≤ 1, 5y(0) + y'(0) = 0, y(1) = 32$
"""

# ╔═╡ 0ab23d66-befe-4d97-8444-03fa65b40c15
α₀, α₁ = 5, 1

# ╔═╡ 5ee16b5d-6e37-4e6a-a009-07bea790e168
β₀, β₁ = 1, 0

# ╔═╡ 23ecd47b-5068-450e-99d5-7579913fe341
γ₀, γ₁ = 0., 32.

# ╔═╡ a562981b-6532-46b3-b4d2-0399c832308a
md"""
## $Z₁, Z₂$
"""

# ╔═╡ 274aef89-c82f-4fe0-a295-b57c0820d5b8
md"""
$Z₁' = -Z₁² - (x+1) Z₁ + \frac{20}{(x+1)^2}$

$Z₁(0) = -5$
"""

# ╔═╡ b7679d80-65ab-4b94-8aa7-8eea19ee1f13
f₁(x, y) = -y^2 - (x+1) * y + 20/((x+1)^2)

# ╔═╡ ce095744-b0fa-4f93-9fff-c62267287e9b
x, Z₁ = euler(
	(x, y) -> -20/((x+1)^2)*y^2 + (x+1) * y + 1, 
	0., 1., 
	0., -5., 
	0.1
)

# ╔═╡ d83348b1-95f7-491e-ba0b-54d3e5bb44b1
20/((x.+1).^2)

# ╔═╡ 0e38eeec-e16e-4408-8459-a4cfc63db8c3
map(X->20/X, (x.+1).^2)

# ╔═╡ 85d8fcb0-e8ec-4c8f-aa15-1012d7b874c4
md"""
$Z₂' = -Z₂(Z₁ + (x+1)) + 5 (x+1)^5$

$Z₂(0) = 0$
"""

# ╔═╡ 872f4b04-91a4-4567-aae1-01a40d54f4d0


# ╔═╡ 8a1e40be-5516-4dac-a5f5-59f1da59bec0


# ╔═╡ 2e49e3c5-ebd2-4d21-bdf9-25c6c51dc29f


# ╔═╡ Cell order:
# ╠═67d84dcc-5077-11ed-0ed2-c742f73ce6d5
# ╠═9b13725c-0fdb-421f-bdc7-df4f02aadd9f
# ╠═9105fc1a-9f8a-4aa3-b242-106715871524
# ╠═0ab23d66-befe-4d97-8444-03fa65b40c15
# ╠═5ee16b5d-6e37-4e6a-a009-07bea790e168
# ╠═23ecd47b-5068-450e-99d5-7579913fe341
# ╟─a562981b-6532-46b3-b4d2-0399c832308a
# ╠═274aef89-c82f-4fe0-a295-b57c0820d5b8
# ╠═b7679d80-65ab-4b94-8aa7-8eea19ee1f13
# ╠═ce095744-b0fa-4f93-9fff-c62267287e9b
# ╠═d83348b1-95f7-491e-ba0b-54d3e5bb44b1
# ╠═0e38eeec-e16e-4408-8459-a4cfc63db8c3
# ╠═85d8fcb0-e8ec-4c8f-aa15-1012d7b874c4
# ╠═872f4b04-91a4-4567-aae1-01a40d54f4d0
# ╠═8a1e40be-5516-4dac-a5f5-59f1da59bec0
# ╠═2e49e3c5-ebd2-4d21-bdf9-25c6c51dc29f
