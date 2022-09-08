### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ fc857e7a-2f75-11ed-2a90-7bffddfb116d
begin
	import Pkg;
	Pkg.activate();

	using Plots;
	using OrderedCollections;
end;

# ╔═╡ 195b797b-6626-4238-b636-aa0fbb54918c
md"""
# Метод Ньютона решения нелинейных уравнений

Вариант 85 (c 120)

$f(x) = x^3 - 0.1x^2 + 0.3x - 0.6 = 0$

$f'(x) = 3x^2 - 0.2x + 0.3$

$f''(x) = 6x - 0.2$
"""

# ╔═╡ e8b5f0f1-f0c1-401c-987b-caf58c6a359d
f(x) = x^3 - 0.1x^2 + 0.3x - 0.6

# ╔═╡ 21442f03-050b-4d0c-a747-368337a2326c
df(x) = 3x^2 - 0.2x + 0.3

# ╔═╡ 676d0dfd-4347-477a-b6b2-1131ba19a78e
ddf(x) = 6x - 0.2

# ╔═╡ dc96a872-4ba0-496f-b2ff-1aa08b1825cc
plot(
	0:.01:1,
	f.(0:.01:1);
	label=:none,
	ratio=:equal
)

# ╔═╡ 4bb3bcb8-e065-467d-a072-0b0546ba33bb
md"""
## Проверка условий сходимости
"""

# ╔═╡ d502a4dc-a0ff-43b7-b9da-4395d94bc688
x₀, δ = .7, .15

# ╔═╡ 45aeace4-a353-4139-b40f-6563fee44a3c
x₁, x₂ = x₀ - δ, x₀ + δ

# ╔═╡ 69c81197-093c-42ad-ba19-6656d978b325
md"""
### 1. $|f''(x)| \le K$
"""

# ╔═╡ 5dfd7e4b-dda1-4342-a179-2c6517327652
K = 5

# ╔═╡ 5c0cad29-ec82-42f1-a39a-510f416e5da5
begin
	plot(
		-1:.1:2,
		abs.(ddf.(-1:.1:2));
		label=:none
	)
	vline!([x₁, x₂];label=:none)
	hline!([K];label=:none)
end

# ╔═╡ bb375808-4182-445b-b71e-6db4af18ba05
md"""
### 2. $f'(x₀) ≠ 0, \frac{1}{|f'(x₀)|} ≤ B$
"""

# ╔═╡ bf7ff886-f14a-4ce0-8d3a-27e756ff9bf2
1 / abs(df(x₀))

# ╔═╡ e2ef6521-b6c9-4fd9-be79-6f7a1a63c375
B = .62

# ╔═╡ 35ef2457-24db-4c8f-a072-47f08e00b90c
md"""
### 3. $\left| \frac{f(x₀)}{f'(x₀)} \right| ≤ η$
"""

# ╔═╡ 6069b08f-2170-452f-9998-54a8d9035dae
abs(f(x₀) / df(x₀))

# ╔═╡ d7443960-704d-4719-a821-b0a5a925e1ba
η = .06

# ╔═╡ 1c8d870a-cdaa-41c5-bbe3-61a7799f3423
md"""
### 4. $h = K * B * η ≤ \frac{1}{2}$
"""

# ╔═╡ 796362b1-4004-4e35-8134-d6ec383ee35c
h = K * B * η

# ╔═╡ 63529b2e-fce3-4965-95fc-ef9d5bca1bba
md"""
### 5. $\frac{1 - \sqrt{1 - 2 h}}{h} η ≤ δ$
"""

# ╔═╡ 185fc047-aa36-4748-98f7-a24177f0f9a3
(1 - sqrt(1 - 2h)) / h * η

# ╔═╡ 5975bd35-7cfb-4b6d-aad4-ccf80360c224
(1 - sqrt(1 - 2h)) / h * η ≤ δ

# ╔═╡ 36e6ca87-adee-4a98-9773-fbd506a9dad9
md"""
## Итерационный процесс методом Ньютона
"""

# ╔═╡ 24763317-00dd-412b-91b3-f1fc19ce2db7
φ(x) = x - f(x) / df(x)

# ╔═╡ e21152a9-464a-4f48-87be-998c2ac8b6f9
XN = [x₀]

# ╔═╡ 8b892ef0-40ab-4f77-8d38-57f0dcf10d04
OrderedDict(XN .=> f.(XN))

# ╔═╡ 3a66f514-1325-4095-a278-f4daba14868d
f(XN[end])

# ╔═╡ 66afe089-a7c8-4d03-a2fd-651ec4b0822b
md"""
## Итерационный процесс методом секущих
"""

# ╔═╡ 48129eed-91b2-4da0-bfdb-11bfb577e6ca
φ(x::Vector) = x[end] - f(x[end]) * (x[end] - x[end-1]) / (f(x[end]) - f(x[end-1]))

# ╔═╡ 0585d71f-2819-471f-b23f-5bede575e761
push!(XN, φ(XN[end]))

# ╔═╡ 3ee935e8-1cc8-4e79-a042-550eaf6cba7b
while abs(XN[end] - XN[end-1]) > 10^(-5)
	push!(XN, φ(XN[end]))
end

# ╔═╡ d12a4ab9-3323-44a5-954b-9d79a250a55a
XS = [.5, 1]

# ╔═╡ 7636376a-4d54-439a-b643-ee4e066ce9a1
push!(XS, φ(XS))

# ╔═╡ f16a0103-997d-4151-8989-170f52aaa424
while abs(XS[end] - XS[end-1]) > 10^(-5)
	push!(XS, φ(XS[end]))
end

# ╔═╡ 0c7e2086-1e4c-497e-8cd7-b41b4861eebc
OrderedDict(XS .=> f.(XS))

# ╔═╡ 403eeeb2-d6e5-4258-afca-ddcae963be33
f(XS[end])

# ╔═╡ Cell order:
# ╠═fc857e7a-2f75-11ed-2a90-7bffddfb116d
# ╟─195b797b-6626-4238-b636-aa0fbb54918c
# ╠═e8b5f0f1-f0c1-401c-987b-caf58c6a359d
# ╠═21442f03-050b-4d0c-a747-368337a2326c
# ╠═676d0dfd-4347-477a-b6b2-1131ba19a78e
# ╠═dc96a872-4ba0-496f-b2ff-1aa08b1825cc
# ╟─4bb3bcb8-e065-467d-a072-0b0546ba33bb
# ╠═d502a4dc-a0ff-43b7-b9da-4395d94bc688
# ╠═45aeace4-a353-4139-b40f-6563fee44a3c
# ╟─69c81197-093c-42ad-ba19-6656d978b325
# ╠═5dfd7e4b-dda1-4342-a179-2c6517327652
# ╠═5c0cad29-ec82-42f1-a39a-510f416e5da5
# ╟─bb375808-4182-445b-b71e-6db4af18ba05
# ╠═bf7ff886-f14a-4ce0-8d3a-27e756ff9bf2
# ╠═e2ef6521-b6c9-4fd9-be79-6f7a1a63c375
# ╟─35ef2457-24db-4c8f-a072-47f08e00b90c
# ╠═6069b08f-2170-452f-9998-54a8d9035dae
# ╠═d7443960-704d-4719-a821-b0a5a925e1ba
# ╟─1c8d870a-cdaa-41c5-bbe3-61a7799f3423
# ╠═796362b1-4004-4e35-8134-d6ec383ee35c
# ╟─63529b2e-fce3-4965-95fc-ef9d5bca1bba
# ╠═185fc047-aa36-4748-98f7-a24177f0f9a3
# ╠═5975bd35-7cfb-4b6d-aad4-ccf80360c224
# ╟─36e6ca87-adee-4a98-9773-fbd506a9dad9
# ╠═24763317-00dd-412b-91b3-f1fc19ce2db7
# ╠═e21152a9-464a-4f48-87be-998c2ac8b6f9
# ╠═0585d71f-2819-471f-b23f-5bede575e761
# ╠═3ee935e8-1cc8-4e79-a042-550eaf6cba7b
# ╠═8b892ef0-40ab-4f77-8d38-57f0dcf10d04
# ╠═3a66f514-1325-4095-a278-f4daba14868d
# ╟─66afe089-a7c8-4d03-a2fd-651ec4b0822b
# ╠═48129eed-91b2-4da0-bfdb-11bfb577e6ca
# ╠═d12a4ab9-3323-44a5-954b-9d79a250a55a
# ╠═7636376a-4d54-439a-b643-ee4e066ce9a1
# ╠═f16a0103-997d-4151-8989-170f52aaa424
# ╠═0c7e2086-1e4c-497e-8cd7-b41b4861eebc
# ╠═403eeeb2-d6e5-4258-afca-ddcae963be33
