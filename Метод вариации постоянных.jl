### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 637c75a4-4b2b-11ed-07f0-dbcdce4fb16f
begin
	import Pkg; 	Pkg.activate()

	using DataFrames
	using Plots
end

# ╔═╡ 7683f22b-3939-4f39-acba-82899b7c60a4
function euler(s, a, b, x₀, T₀, S₀, h)
	x = [x₀]
	T = [T₀]
	S = [S₀]

	while x[end]+h < b
		push!(T, T[end] + h * S[end])
		push!(S, S[end] + h * s(x[end], T[end-1], S[end]))
		push!(x, x[end] + h)
	end

	return x, T, S
end

# ╔═╡ ab4c1a90-3d8d-44cf-945a-b8540e1d4aaa
md"""
# Метод вариации постоянных

$y'' + \frac{1}{x}y' - 2y = x^2$

$0.5 ≤ x ≤ 1, y'(0.5) = -0.5, y'(1) = -1$
"""

# ╔═╡ 4841bc0e-e8a6-4b17-8ccb-24b6713cd905
α₀, α₁ = 0, 0

# ╔═╡ fcbb49fc-3f37-46e5-acee-30ed4478b716
β₀, β₁ = 1, 1 

# ╔═╡ 95a5c54a-af9c-45b8-8fbf-2a7ee2116abe
γ₀, γ₁ = -.5, -1.

# ╔═╡ 358ca7f3-c745-47ea-a209-54fcdea4d54e
md"""
## $Z, Z₁, Z₂$
"""

# ╔═╡ e301f4d4-7a73-4961-9910-d7102773e09b
md"""
$\begin{cases}
	Z''(x) + \frac{1}{x}Z'(x) - 2Z(x) = x^2, 0.5 ≤ x ≤ 1 \\
	Z(0.5) = 0, Z'(0.5) = 0
\end{cases}$

$\begin{cases}
	T'(x) = S(x) \\
	S'(x) = - \frac{1}{x}S(x) + 2T(x) + x^2, 0.5 ≤ x ≤ 1 \\
	T(0.5) = 0, S(0.5) = 0
\end{cases}$
"""

# ╔═╡ 6426609a-f117-4404-b5e2-6bae2164a547
x, T, S = euler(
	(x, t, s) -> -s/x + 2*t + x^2,
	0.5, 1.,
	0.5,
	.0, .0,
	.1
)

# ╔═╡ dee602f5-9406-4743-a57a-9fccc55d958b
md"""
$\begin{cases}
	Z₁''(x) + \frac{1}{x}Z₁'(x) - 2Z₁(x) = 0, 0.5 ≤ x ≤ 1 \\
	Z₁(0.5) = 0, Z₁'(0.5) = 1
\end{cases}$

$\begin{cases}
	T₁'(x) = S₁(x) \\
	S₁'(x) = - \frac{1}{x}S₁(x) + 2T₁(x), 0.5 ≤ x ≤ 1 \\
	T₁(0.5) = 0, S₁(0.5) = 1
\end{cases}$
"""

# ╔═╡ 14efdc2e-6737-4550-a610-4a57b1a6e3b1
T₁, S₁ = euler(
	(x, t, s) -> -s/x + 2*t,
	0.5, 1.,
	0.5,
	.0, 1.,
	.1
)[2:3]

# ╔═╡ ce637c87-1d2e-4bf4-96c7-f0d9f30043fd
md"""
$\begin{cases}
	Z₂''(x) + \frac{1}{x}Z₂'(x) - 2Z₂(x) = 0, 0.5 ≤ x ≤ 1 \\
	Z₂(0.5) = 1, Z₂'(0.5) = 0
\end{cases}$

$\begin{cases}
	T₂'(x) = S₂(x) \\
	S₂'(x) = - \frac{1}{x}S₂(x) + 2T₂(x), 0.5 ≤ x ≤ 1 \\
	T₂(0.5) = 1, S₂(0.5) = 0
\end{cases}$
"""

# ╔═╡ 410c92a7-0ee0-4de8-b90e-86b33ba2e630
T₂, S₂ = euler(
	(x, t, s) -> -s/x + 2*t,
	0.5, 1.,
	0.5,
	1., .0,
	.1
)[2:3]

# ╔═╡ f99495f0-a692-4ea1-9dbb-c20437f7ecf7
md"""
## $C₁, C₂$
"""

# ╔═╡ fe8d65ec-d598-4dfe-8abe-6f9bd4f44c79
A = [
	α₀*T₁[1] + β₀*S₁[1] α₀*T₂[1] + β₀*S₂[1] γ₀ - α₀*T[1] - β₀*S[1]
	α₁*T₁[end] + β₁*S₁[end] α₁*T₂[end] + β₁*S₂[end] γ₁ - α₁*T[end] - β₁*S[end]
]

# ╔═╡ 9705c731-9226-4b06-b502-51111326022f
C₁, C₂ = A[:, [1,2]] \ A[:, 3]

# ╔═╡ 89082044-4330-406f-8323-c33db2da3d35
md"""
## $yₖ$
"""

# ╔═╡ 6d46fda2-e349-4517-b9e2-3f8802285e74
yₖ = C₁.*T₁ + C₂.*T₂ + T

# ╔═╡ 3b609ca6-9310-4e06-bd37-bd9e62da5e07
DataFrame(
	xₖ=x,
	T=T, S=S,
	T₁=T₁, S₁=S₁,
	T₂=T₂, S₂=S₂,
	yₖ=yₖ
)

# ╔═╡ 5869e7b3-1b49-47e1-a132-3fc4f958a1b8
begin
	plot(
		x, yₖ;
		label=:none
	)
	scatter!(
		x, yₖ;
		label=:none
	)
end

# ╔═╡ Cell order:
# ╠═637c75a4-4b2b-11ed-07f0-dbcdce4fb16f
# ╠═7683f22b-3939-4f39-acba-82899b7c60a4
# ╟─ab4c1a90-3d8d-44cf-945a-b8540e1d4aaa
# ╠═4841bc0e-e8a6-4b17-8ccb-24b6713cd905
# ╠═fcbb49fc-3f37-46e5-acee-30ed4478b716
# ╠═95a5c54a-af9c-45b8-8fbf-2a7ee2116abe
# ╠═358ca7f3-c745-47ea-a209-54fcdea4d54e
# ╟─e301f4d4-7a73-4961-9910-d7102773e09b
# ╠═6426609a-f117-4404-b5e2-6bae2164a547
# ╠═dee602f5-9406-4743-a57a-9fccc55d958b
# ╠═14efdc2e-6737-4550-a610-4a57b1a6e3b1
# ╟─ce637c87-1d2e-4bf4-96c7-f0d9f30043fd
# ╠═410c92a7-0ee0-4de8-b90e-86b33ba2e630
# ╟─f99495f0-a692-4ea1-9dbb-c20437f7ecf7
# ╠═fe8d65ec-d598-4dfe-8abe-6f9bd4f44c79
# ╠═9705c731-9226-4b06-b502-51111326022f
# ╟─89082044-4330-406f-8323-c33db2da3d35
# ╠═6d46fda2-e349-4517-b9e2-3f8802285e74
# ╠═3b609ca6-9310-4e06-bd37-bd9e62da5e07
# ╟─5869e7b3-1b49-47e1-a132-3fc4f958a1b8
