### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 1a132654-3ff7-11ed-2942-4f57d8b62b80
begin 
	import Pkg; 	Pkg.activate();	

	using PlutoUI;
	TableOfContents();
end

# ╔═╡ d21fc2dd-b81d-4984-af26-a69e517db48c
md"""
# Метод Рунге-Кутта

$\begin{cases}
	y' = -0.6 x y^2 + z^2 - x^2 - 1 \\
	z' = \frac{1}{0.3 z^2} - y - \frac{x}{z}
\end{cases}$

$\begin{cases}
	y(0) = \frac{1}{0.3} \\
	z(0) = 1
\end{cases}$

$a=0, b=0.5, h=0.1$
"""

# ╔═╡ 3f68e659-adb0-4327-9871-d604e8bc98ac
a, b = 0., 0.5

# ╔═╡ 00d31461-a20e-4d95-92e5-9c2cffe36929
h = 0.1

# ╔═╡ 086cc7dc-54e7-457c-92fb-34a013b1cd6c
f(x, y) = [
	-0.6*x*y[1]^2 + y[2]^2 - x^2 - 1,
	1/(.3*y[2]^2) - y[1] - x/y[2]
] 

# ╔═╡ c3f9d781-f7b1-4a7d-b2b6-f433d6ee35d8
x = [a]

# ╔═╡ 9aead913-6ff1-43f6-bbfb-1a079c7dde61
y = [[1/0.3, 1]]

# ╔═╡ 4ffefe5f-4532-4c9f-acfa-d0355c8cadbd
for i in 1:round(Int, (b-a)/h)
	k₁ = h * f(x[end], y[end])
	k₂ = h * f(h/2 + x[end], k₁/2 + y[end])
	k₃ = h * f(h/2 + x[end], k₂/2 + y[end])
	k₄ = h * f(h/2 + x[end], k₃/2 + y[end])

	push!(x, x[end]+h)
	push!(y, y[end]+(k₁ + 2k₂ + 2k₃ + k₄)/6)
end

# ╔═╡ 25ae0b17-0e6a-4ee4-83c7-b46362437d32
x

# ╔═╡ d6460a10-3db9-4e7b-b039-4399155c39ed
y

# ╔═╡ c1cedd78-c13e-460a-bd2f-25f94f39dcb2
x .=> y

# ╔═╡ 15cedf0e-d7cc-428f-8383-1df801b1a306


# ╔═╡ 07b616a9-1ae7-4173-8821-369d289f6c80


# ╔═╡ Cell order:
# ╠═1a132654-3ff7-11ed-2942-4f57d8b62b80
# ╟─d21fc2dd-b81d-4984-af26-a69e517db48c
# ╠═3f68e659-adb0-4327-9871-d604e8bc98ac
# ╠═00d31461-a20e-4d95-92e5-9c2cffe36929
# ╠═086cc7dc-54e7-457c-92fb-34a013b1cd6c
# ╠═c3f9d781-f7b1-4a7d-b2b6-f433d6ee35d8
# ╠═9aead913-6ff1-43f6-bbfb-1a079c7dde61
# ╠═4ffefe5f-4532-4c9f-acfa-d0355c8cadbd
# ╠═25ae0b17-0e6a-4ee4-83c7-b46362437d32
# ╠═d6460a10-3db9-4e7b-b039-4399155c39ed
# ╠═c1cedd78-c13e-460a-bd2f-25f94f39dcb2
# ╠═15cedf0e-d7cc-428f-8383-1df801b1a306
# ╠═07b616a9-1ae7-4173-8821-369d289f6c80
