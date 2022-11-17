### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 163677d0-3a77-11ed-2e43-5d50845473cf
begin
	import Pkg; 	Pkg.activate()

	
end

# ╔═╡ 71d836dc-3160-4073-99f0-2fb198a61a4b
md"""
# Метод Эйлера с заданной точностью

Вариант 6

$y' = \frac{6 - x^2 y^2}{-x^2}$

$y(1) = 2, a = 1, b = 1.5, h = 0.05$

Отрезок: $[1, 1.5]$
"""

# ╔═╡ 49ec8cc8-0a22-4f48-9feb-f4fc2bb2b263
f(x, y) = (6 - x^2 * y^2) / (- x^2)

# ╔═╡ e1b21390-0acd-4cff-b69e-5744346ee320
x₀ = 1

# ╔═╡ e416cd52-4e65-462d-9d63-d603f62cdf1c
y₀ = 2

# ╔═╡ 86a5a41b-bb9f-407f-a3b3-93280f8f0660
a, b = 1, 1.5

# ╔═╡ ba9a6a9d-03c3-4069-88a7-10066920d417
ϵ = 10^-5

# ╔═╡ d075a9f3-9133-4ffc-802e-938e2b81f0cd
h = 0.05

# ╔═╡ ea0455f6-07e2-418a-b1f3-11c1d593a57d
x = [x₀]

# ╔═╡ 3f82c3e0-abd7-42d0-baa7-85a9649362bd
y = [(y₀, y₀)]

# ╔═╡ 94ceedfa-988e-48b9-b44b-9bb092e16bec
segment = []

# ╔═╡ f5ec37ee-e54e-4d51-96e8-33aff2cc0569
function sol(x=[1.], y=[(2., 2.)], b=1.5, h=0.05)
	segment = []
	H = []
	while x[end] <= b
		xₕ₂ = x[end] + h/2
		yₕ₂ = y[end][1] + h/2 * f(x[end], y[end][1])
		yₕ₂ = yₕ₂ + h/2 * f(xₕ₂, y[end][1])
	
		xₕ = x[end] + h
		yₕ = y[end][1] + h * f(x[end], y[end][1])

		if abs(yₕ - yₕ₂) < ϵ
			push!(segment, [x[end], xₕ])
			push!(H, h)
			push!(x, xₕ)
			push!(y, (yₕ, yₕ₂))
		else
			h /= 2
		end
	end

	return (segment, H, y[2:(end)])
end

# ╔═╡ 0fcec575-bf75-4251-ab55-62c8f9d188d9
solution = sol()

# ╔═╡ 1739502f-cc01-418b-9ca2-e0a354f9263d
solution .|> length

# ╔═╡ da196feb-3fd6-4886-a445-e23a37a14ac2
for i in 1:length(solution[1])
	println("Отрезок: $(solution[1][i]), h = $(solution[2][i]), yₕ = $(solution[3][i][1]), yₕ₂ = $(solution[3][i][2])")
end

# ╔═╡ 4757aba2-087d-4810-a697-2cb2e669a30e


# ╔═╡ 0f10cc4c-5a8c-4d05-adb9-fb7e5145888e


# ╔═╡ Cell order:
# ╠═163677d0-3a77-11ed-2e43-5d50845473cf
# ╟─71d836dc-3160-4073-99f0-2fb198a61a4b
# ╠═49ec8cc8-0a22-4f48-9feb-f4fc2bb2b263
# ╠═e1b21390-0acd-4cff-b69e-5744346ee320
# ╠═e416cd52-4e65-462d-9d63-d603f62cdf1c
# ╠═86a5a41b-bb9f-407f-a3b3-93280f8f0660
# ╠═ba9a6a9d-03c3-4069-88a7-10066920d417
# ╠═d075a9f3-9133-4ffc-802e-938e2b81f0cd
# ╠═ea0455f6-07e2-418a-b1f3-11c1d593a57d
# ╠═3f82c3e0-abd7-42d0-baa7-85a9649362bd
# ╠═94ceedfa-988e-48b9-b44b-9bb092e16bec
# ╠═f5ec37ee-e54e-4d51-96e8-33aff2cc0569
# ╠═0fcec575-bf75-4251-ab55-62c8f9d188d9
# ╠═1739502f-cc01-418b-9ca2-e0a354f9263d
# ╠═da196feb-3fd6-4886-a445-e23a37a14ac2
# ╠═4757aba2-087d-4810-a697-2cb2e669a30e
# ╠═0f10cc4c-5a8c-4d05-adb9-fb7e5145888e
