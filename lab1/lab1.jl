### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 0cb92a82-8e93-11ec-091b-f1abb82c2225
begin
	import Pkg;
	Pkg.activate()

	using LinearAlgebra
	using DataFrames, PrettyTables

	using PlutoUI
	TableOfContents()
end

# ╔═╡ 32f02311-d88b-47ef-a399-59752ce037a5
md"""
# Лабораторная работа #1
# Решение СЛАУ методом Гаусса
"""

# ╔═╡ 2fd3c08f-0768-4533-90be-d0f457f77bdf
begin
	∞(A::Matrix{<:Number}) = maximum( sum(abs, A, dims=2) )
	∞(v::Vector{<:Number}) = maximum( abs, v )
end

# ╔═╡ 90464d1e-75aa-423c-9f29-dce94cf84003
begin
	|₁(A::Matrix{<:Number}) = maximum( sum(abs, A, dims=1) )
	|₁(v::Vector{<:Number}) = sum( abs, v )
end

# ╔═╡ 8d82ca3e-b058-49b7-8398-064de7da35d3
begin
	|₂(A::Matrix{<:Number}) = maximum( eigen(A' * A).values ) |> sqrt
	|₂(v::Vector{<:Number}) = sum( abs2, v ) |> sqrt
end

# ╔═╡ f3572ec8-4ec9-4b87-b39b-34f3f652c8a0
md"""
## Задание 1. Метод Гаусса
"""

# ╔═╡ f1857dc6-6a51-4bed-b2a0-749c325b7794
D = [
	6.22 1.42 -1.72 1.91
	1.44 5.33 1.11 -1.82
	-1.72 1.11 5.24 1.42
	1.91 -1.82 1.42 6.55
]

# ╔═╡ 6bbb66c7-5458-4482-93e3-06293bf7dbb7
A = D + I

# ╔═╡ 53a63945-6809-4fc4-8987-21bd7b7c5e90
f̄ = [7.53, 6.06, 8.05, 8.06]

# ╔═╡ f2797a5d-8f42-4c03-ac10-4331c1075ff1
begin
	mutable struct SLE
		A::Matrix{<:Number}
		f̄::Vector{<:Number}

		solution::DataFrame

		x
		r̄

		det

		function SLE(A::Matrix{<:Number}, f̄::Vector{<:Number})
			(size(A, 1) != size(A, 2) || size(A, 1) != length(f̄))  && throw("Dimension mismatch")

			new(A, f̄)
		end
	end

	function solve(eq::SLE)
		eq.det = 1
		n = length(f̄)
		eq.solution = DataFrame(
			Matrix{Union{Missing, Number}}(missing, sum(2:n)+2n, 2+n+3),
			[["Step", "i"]; [ "aᵢ$(i)" for i=1:(n+1) ]; ["Σ", "Σᵢ"] ]
		)
		b = Int[]
		for step=1:n#(n-1)
			start = sum((n+3-step):(n+1))
			
			eq.solution[start+1, "Step"] = step
			eq.solution[(start+1):(n-step+start+1), "i"] = step:n

			if step == 1
				eq.solution[1:n, 3:(n+1+2)] = [eq.A eq.f̄]
			else
				push!(b, start)
				eq.solution[(start+1):(start+n+1-step) , (2+step):(n+1+2+1)] = (eq.solution[(start-n-1+step):(start-1), (2+step):(n+1+2+1)] |> Matrix) .- (eq.solution[(start-n-1+step):(start-1), step+1] |> Vector) * (eq.solution[start, (step+2):(n+1+2+1)] |> Vector)'
			end

			# det
			eq.det *= eq.solution[start+1, (2+step)]
			
			# Σ
			eq.solution[(start+1):(start+n+1-step), "Σ"] = sum( eq.solution[(start+1):(start+n+1-step), (2+step):(n+1+2)] |> Matrix, dims=2 )[:, 1]
			# b
			eq.solution[start+n+2-step, (2+step):(n+1+2+1)] = (eq.solution[start+1, (2+step):(n+1+2+1)] |> Vector) ./ eq.solution[start+1, (2+step)]
			# Σᵢ
			eq.solution[start+n+2-step, "Σᵢ"] = sum( eq.solution[start+n+2-step, (2+step):(n+1+2)] |> Vector)
			eq.solution[(start+1):(start+n+1-step), "Σᵢ"] = sum( eq.solution[(start+1):(start+n+1-step), (step+2):(n+1+2)] |> Matrix, dims=2)[:, 1]
		end

		start = sum(2:(n+1))-1
		eq.solution[start+1, :] .= missing
		eq.solution[start+1, "Step"] = n+1
		
		eq.solution[start+1, "aᵢ$(n)"] = 1
		# x
		eq.solution[start+1, "aᵢ$(n+1)"] = eq.solution[start, "aᵢ$(n+1)"]/eq.solution[start, "aᵢ$(n)"]
		# x̃
		eq.solution[start+1, "Σ"] = eq.solution[start, "Σ"]/eq.solution[start, "aᵢ$(n)"]
		for i=(n-1):-1:1
			eq.solution[start+n-i+1, 2+i] = 1
			# x
			eq.solution[start+n-i+1, "aᵢ$(n+1)"] = eq.solution[b[i], "aᵢ$(n+1)"] - eq.solution[(start+n-i):-1:(start+1), "aᵢ$(n+1)"] ⋅ (eq.solution[b[i], (i+3):(n+2)] |> Vector)
			# x̃
			eq.solution[start+n-i+1, "Σ"] = eq.solution[b[i], "Σ"] - eq.solution[(start+n-i):-1:(start+1), "Σ"] ⋅ (eq.solution[b[i], (i+3):(n+2)] |> Vector)
		end

		eq.x = eq.solution[(start+1):(start+n), "aᵢ$(n+1)"] |> Vector |> reverse .|> Float64

		eq.r̄ = (eq.A * eq.x - eq.f̄)

		with_terminal() do
			pretty_table(eq.solution, formatters=ft_nomissing)
		end
	end
end;

# ╔═╡ d052cb09-e039-4c65-96ad-e14c32a0fbbc
sle = SLE(A, f̄)

# ╔═╡ cd9d8f24-8bf7-49ee-a506-4a98c4c0a5d8
sle |> solve

# ╔═╡ cdab7141-216c-4f72-8c29-f2f5cd013670
pretty_table(sle.solution, formatters=ft_nomissing)

# ╔═╡ ba30acac-a40f-4034-b16d-8d4292477d1e
sle.x

# ╔═╡ 7a4ea9df-c08b-4ef1-8f1e-df476d8698a3
A \ f̄

# ╔═╡ 98618ada-6cc9-4436-802b-b66acda14603
md"""
## Задание 2. Невязка
### Кубическая норма
"""

# ╔═╡ e207693c-1d52-4ee4-92b2-8799a26decc2
∞(sle.r̄)

# ╔═╡ 3c693e33-245e-4eb5-a349-6eaa77dc5974
sle.r̄

# ╔═╡ d1712dc8-eeaa-4b55-b03b-de381f111ad6
md"""
### Октаэдрическая норма
"""

# ╔═╡ 0cc9be03-889f-411f-8a97-12f9ed0363f1
|₁(sle.r̄)

# ╔═╡ 3055ffbe-5012-4b65-8c6d-7fd33dd55661
md"""
### Сферическая норма
"""

# ╔═╡ 2159d142-9448-4ec9-ae27-abf34c797bbb
|₂(sle.r̄)

# ╔═╡ 330340cb-fc9a-456c-b292-5058ecdd968a
md"""
## Задание 3. Определитель матрицы $A$
"""

# ╔═╡ b1843d9f-4b38-43ec-a069-0a24a46a9773
sle.det

# ╔═╡ 4da8e181-f219-43e9-8747-58ed15eff730
det(sle.A)

# ╔═╡ 0e3dde33-49e7-4b4d-ba0d-ab68f1b4e382
md"""
## Задание 4. Матричные нормы
### Кубическая норма
"""

# ╔═╡ c3145d32-6325-477b-bbcd-0d7809b837ca
∞(A)

# ╔═╡ ae51d816-32f5-4bd6-9805-a53b556ac05c
md"""
### Октоэдрическая норма
"""

# ╔═╡ 2d89109b-dd51-4a22-acab-ed8bea61880c
|₁(A)

# ╔═╡ 73340101-8c51-4f19-a200-58c99002138f
md"""
### Сферическая норма
"""

# ╔═╡ df3c1167-ed28-4b77-8ac7-0f8ebddee9ae
|₂(A)

# ╔═╡ 84472daa-bd6b-4f58-a2d1-dbbfb6a222bb
md"""
# Лабораторная работа #2
# Нахождение обратной матрицы
"""

# ╔═╡ 436cf432-e6a5-4029-abf1-15e73833f689
md"""
## Задание 1. Нахождение обратной матрицы
"""

# ╔═╡ 073884d6-4126-4001-bb23-cc84b678cf7d
map(
	f -> begin
		sle = SLE(A, f |> Vector)
		solve(sle)
		sle
	end,
	eachcol(diagm(ones(size(A, 1))))
)

# ╔═╡ 8c0747d9-84b0-44a1-834e-ad5e02cd5a36
inv_sle(A::Matrix) = hcat(map(
	f -> begin
		sle = SLE(A, f |> Vector)
		solve(sle)
		sle.x
	end,
	eachcol(diagm(ones(size(A, 1))))
)...)

# ╔═╡ 5eea2411-64cc-40b9-a36b-3106f7ff26d0
A * inv_sle(A)

# ╔═╡ e8ced749-64a9-4335-921c-e8e78764a14d
sum(abs, A * inv(A) - I) 

# ╔═╡ c4604f18-cc84-4a85-96d1-456423c6edb5
md"""
## Задание 2. Число обусловленности
"""

# ╔═╡ a5fbe839-56cd-4338-be35-eded3046e3d9
cond(A::Matrix; norm::Function=|₂) =
	(norm(A), norm(inv_sle(A)), norm(A) * norm(inv_sle(A)))

# ╔═╡ 515370c4-8911-4aab-807a-89655aaaf341
cond(A)

# ╔═╡ b4d3a49d-368a-42db-afcc-4acce7c754eb
cond(A, norm=|₁)

# ╔═╡ 1977f686-9a9b-4e70-84b4-dd1dc394158b
cond(A, norm=∞)

# ╔═╡ Cell order:
# ╟─0cb92a82-8e93-11ec-091b-f1abb82c2225
# ╟─32f02311-d88b-47ef-a399-59752ce037a5
# ╠═2fd3c08f-0768-4533-90be-d0f457f77bdf
# ╠═90464d1e-75aa-423c-9f29-dce94cf84003
# ╠═8d82ca3e-b058-49b7-8398-064de7da35d3
# ╟─f3572ec8-4ec9-4b87-b39b-34f3f652c8a0
# ╠═f2797a5d-8f42-4c03-ac10-4331c1075ff1
# ╠═f1857dc6-6a51-4bed-b2a0-749c325b7794
# ╠═6bbb66c7-5458-4482-93e3-06293bf7dbb7
# ╠═53a63945-6809-4fc4-8987-21bd7b7c5e90
# ╠═d052cb09-e039-4c65-96ad-e14c32a0fbbc
# ╠═cd9d8f24-8bf7-49ee-a506-4a98c4c0a5d8
# ╠═cdab7141-216c-4f72-8c29-f2f5cd013670
# ╠═ba30acac-a40f-4034-b16d-8d4292477d1e
# ╠═7a4ea9df-c08b-4ef1-8f1e-df476d8698a3
# ╟─98618ada-6cc9-4436-802b-b66acda14603
# ╠═e207693c-1d52-4ee4-92b2-8799a26decc2
# ╠═3c693e33-245e-4eb5-a349-6eaa77dc5974
# ╟─d1712dc8-eeaa-4b55-b03b-de381f111ad6
# ╠═0cc9be03-889f-411f-8a97-12f9ed0363f1
# ╟─3055ffbe-5012-4b65-8c6d-7fd33dd55661
# ╠═2159d142-9448-4ec9-ae27-abf34c797bbb
# ╟─330340cb-fc9a-456c-b292-5058ecdd968a
# ╠═b1843d9f-4b38-43ec-a069-0a24a46a9773
# ╠═4da8e181-f219-43e9-8747-58ed15eff730
# ╟─0e3dde33-49e7-4b4d-ba0d-ab68f1b4e382
# ╠═c3145d32-6325-477b-bbcd-0d7809b837ca
# ╟─ae51d816-32f5-4bd6-9805-a53b556ac05c
# ╠═2d89109b-dd51-4a22-acab-ed8bea61880c
# ╟─73340101-8c51-4f19-a200-58c99002138f
# ╠═df3c1167-ed28-4b77-8ac7-0f8ebddee9ae
# ╟─84472daa-bd6b-4f58-a2d1-dbbfb6a222bb
# ╟─436cf432-e6a5-4029-abf1-15e73833f689
# ╠═073884d6-4126-4001-bb23-cc84b678cf7d
# ╠═8c0747d9-84b0-44a1-834e-ad5e02cd5a36
# ╠═5eea2411-64cc-40b9-a36b-3106f7ff26d0
# ╠═e8ced749-64a9-4335-921c-e8e78764a14d
# ╟─c4604f18-cc84-4a85-96d1-456423c6edb5
# ╠═a5fbe839-56cd-4338-be35-eded3046e3d9
# ╠═515370c4-8911-4aab-807a-89655aaaf341
# ╠═b4d3a49d-368a-42db-afcc-4acce7c754eb
# ╠═1977f686-9a9b-4e70-84b4-dd1dc394158b
