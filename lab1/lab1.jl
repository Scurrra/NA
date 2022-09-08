### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 6d05ed77-1924-4617-978b-c1fbfc0de550
using ImplicitEquations, Plots

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
		n = length(eq.f̄)
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

# ╔═╡ d052cb09-e039-4c65-96ad-e14c32a0fbbc
sle = SLE(A, f̄)

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

# ╔═╡ e8ced749-64a9-4335-921c-e8e78764a14d
sum(abs, A * inv(A) - I) 

# ╔═╡ c4604f18-cc84-4a85-96d1-456423c6edb5
md"""
## Задание 2. Число обусловленности
"""

# ╔═╡ f5271f10-5515-46e6-ae62-cec7eec1b49e
md"""
# Лабораторная работа #3
# Метод квадратных корней (м-д Холецкого)
"""

# ╔═╡ 2865dc6f-5c92-4752-88b8-5edb3e261e39
A₃ = [
	1.65 -1.76 0.77
	-1.76 1.04 -2.61
	0.77 -2.61 -3.18
]

# ╔═╡ 1db51949-df20-44eb-bbec-9047596a1284
f̄₃ = [
	2.15, 0.82, -0.73
]

# ╔═╡ 9f07676a-1777-42d5-aa04-68d1e3aebf86
A₃ \ f̄₃

# ╔═╡ 690da59c-dc1f-415f-89c4-9246216f1f41
md"""
## Задание 1. $A = Aᵀ > 0$
"""

# ╔═╡ b8ad41d5-88ee-4f75-9c54-5578501edf16
A₃ᵀ = A₃ |> transpose |> Matrix

# ╔═╡ d0e2df4d-df6a-4b2e-9843-58e9fcab2825
A₃ᵀ == A₃

# ╔═╡ 41776634-e2f1-427a-b516-ebd69d2293c4
[
	det(A₃[1:i, 1:i])
	for i=1:size(A₃, 1)
]

# ╔═╡ 8c76b049-469c-4663-9b93-af04e1821e92
ḡ₃ = A₃*f̄₃

# ╔═╡ 76d93737-1712-4cdc-ab45-95c70eb752b1
B = A₃ * A₃

# ╔═╡ 0d33adbc-2cc6-4333-9765-ca5a6877a9a9
[
	det(B[1:i, 1:i])
	for i=1:size(B, 1)
]

# ╔═╡ 65ec3210-63e4-45c6-bef6-e1bb45d3d16f
md"""
## Задание 2. $A = L * Lᵀ$
"""

# ╔═╡ 377c0a66-d2f7-445d-8425-ea5fa0d11bce
function LLᵀ(A::Matrix)
	L = zeros(size(A))
	n = size(A, 1)
	#return L

	L[1, 1] = sqrt( A[1, 1] )
	L[2:end, 1] .= A[2:end, 1] / L[1, 1]

	#return L
	
	for k=2:n
		L[k, k] = sqrt( A[k, k] - sum(a->a^2, L[k, 1:(k-1)]) )

		for i=(k+1):n
			L[i, k] = (A[i, k] - L[i, 1:(k-1)]' * L[k, 1:(k-1)]) / L[k, k]
		end
	end

	L, L' |> Matrix
end

# ╔═╡ 38701eb2-8750-4d3f-add0-3d7f0c23d1f5
L, Lᵀ = LLᵀ(B)

# ╔═╡ 5bfe6216-29fb-41be-ba7a-22b8cbd8d1cb
L * Lᵀ

# ╔═╡ 94e4d42c-8ca3-4488-9959-bffe12088477
B

# ╔═╡ b5e0fa99-40de-4ee9-983d-ce4dce9e6dcc
md"""
## Задание 3. Обратный метод Гаусса
"""

# ╔═╡ 50da520d-5de8-4565-b411-9ad2de8a3131
y₁ = ḡ₃[1] / L[1, 1]

# ╔═╡ 5d374e90-202b-4d69-8e4a-f10328af3b30
y₂ = (ḡ₃[2] - L[2, 1]*y₁) / L[2, 2]

# ╔═╡ 6ebc8f74-254f-4211-a34a-9e1f1b2dd93e
y₃ = (ḡ₃[3] - L[3, 1]*y₁ - L[3, 2]*y₂) / L[3, 3]

# ╔═╡ 6d1fc1d3-e69e-47d1-933a-c14e74e52496
x₃ = y₃ / Lᵀ[3, 3]

# ╔═╡ f507f669-e56c-4588-80bb-c47e4076c1c4
x₂ = (y₂ - Lᵀ[2, 3]*x₃) / Lᵀ[2, 2]

# ╔═╡ 8c6a9db4-5509-47b2-8cf1-7857ca1dd6ed
x₁ = (y₁ - Lᵀ[1, 3]*x₃ - Lᵀ[1, 2]*x₂) / Lᵀ[1, 1]

# ╔═╡ 2bad0e5c-b1a7-47c1-b7e4-5824c4cb331c
md"""
## Задание 4. Проверка решения
"""

# ╔═╡ be9cdc8a-f563-4d5f-9037-cb7c13d5dd13
f̄₃

# ╔═╡ 6c4ccf17-99cb-4dab-89a8-3cdf10c049ea
md"""
# Лабораторная работа #4
# Метод простых итераций

## Задание 1

1. Привести систему 
$\begin{cases}
	2x₁ + x₂ + αx₃ = -2.9 \\
	αx₁+5x₂+0.72x₃=-0.7 \\
	-1.2x₁ + 3x₂ + 1.7x₃ = -9.86
\end{cases}\;, α = 0.1*K$

к системе с диагональным преобладанием

$K = 1$

I + III - II -> III

2 III + I - II -> III
"""

# ╔═╡ b166cf1b-d5d7-4b0d-b4ac-54967d3cd8a9
α = 0.1

# ╔═╡ 95fa8cbe-9c17-4277-a17a-a1299378a368
A₄ = [
	2 1 α -2.9
	α 5 0.72 -0.7
	-1.2 3 1.7 -9.86
]

# ╔═╡ 159b490d-a0ee-4094-8bbf-4f5ed6671e6e
begin
	A₄[3,:] = 2 .* A₄[3,:] .+ A₄[1,:] .- A₄[2,:]
	A₄
end

# ╔═╡ 95aa6b37-dfcc-489f-8d9e-29dd424bbdf0
f₄̄ = A₄[:, 4]

# ╔═╡ f667d18f-a070-4921-9f6e-40a05abb13dc
A₄̂ = A₄[:, 1:3]

# ╔═╡ 9afaaf5d-e1a3-46f6-98dd-3e66f5bfa7d1
ϵ = 1/2*10^-4

# ╔═╡ 9615683c-bb4c-4d42-bbe6-5044ab2bedc7
C = inv(diagm(diag(A₄̂)))

# ╔═╡ 63083def-2687-4635-bf08-d82935c62b1f
∞(I - C*A₄̂)

# ╔═╡ dfa13584-a688-4ad4-b1d9-f8cdcc9c4727
function solve(A::Matrix, f̄::Vector; ϵ=ϵ)
	C = A |> diag |> diagm |> inv
	B = I - C * A
	ḡ = C * f̄
	
	x = [ḡ]
	push!(x, B*ḡ + ḡ)

	i = 2
	while true
		if ∞(B) ≤ 1//2
			∞(x[i] - x[i-1]) ≤ ϵ && return B, ḡ, i, x
		else
			∞(x[i] - x[i-1]) ≤ (1-∞(B)) / ∞(B) * ϵ && return B, ḡ, i, x
		end

		push!(x, B*x[i] + ḡ)
		i += 1
	end
	
end

# ╔═╡ cd9d8f24-8bf7-49ee-a506-4a98c4c0a5d8
sle |> solve

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

# ╔═╡ a5fbe839-56cd-4338-be35-eded3046e3d9
cond(A::Matrix; norm::Function=|₂) =
	(norm(A), norm(inv_sle(A)), norm(A) * norm(inv_sle(A)))

# ╔═╡ 515370c4-8911-4aab-807a-89655aaaf341
cond(A)

# ╔═╡ b4d3a49d-368a-42db-afcc-4acce7c754eb
cond(A, norm=|₁)

# ╔═╡ 1977f686-9a9b-4e70-84b4-dd1dc394158b
cond(A, norm=∞)

# ╔═╡ a4b44707-15b9-4533-befe-4681f011c665
ȳ = begin
	sle₃₁ = SLE(L, ḡ₃)
	sle₃₁ |> solve
	sle₃₁.x
end

# ╔═╡ 8d7419f9-2885-426e-9709-fdceac7537e5
sle₃₁.solution

# ╔═╡ dbf0aaf2-9132-4ab4-b2ce-3ae85b1986f5
ȳ

# ╔═╡ 35928c65-a4b5-4c9f-8c05-886fa514a061
x̄ = begin
	sle₃₂ = SLE(Lᵀ, ȳ)
	sle₃₂ |> solve
	sle₃₂.x
end

# ╔═╡ 850a3b72-1f63-463e-8833-9256062711ff
sle₃₂.solution

# ╔═╡ 36446fa9-2848-4147-822a-81ff4acb16e9
x̄

# ╔═╡ 95cd2ca2-4e51-489a-b184-a726c7080b01
A₃ * x̄ 

# ╔═╡ 1e232704-cbd9-4d4f-9731-7bf7d04c37ee
r̄ = A₃ * x̄ - f̄₃

# ╔═╡ e48e96e2-ce80-436f-b065-6abf7c18eeff
|₁(r̄)

# ╔═╡ 4f2181c2-ccc9-42e9-b50f-c3a4f5d44267
|₂(r̄)

# ╔═╡ 890060d2-4307-4c0f-b2c2-d5711c1e8cc1
∞(r̄)

# ╔═╡ 9020a0ba-4320-45b8-b1ae-6dc256929528
B̂, ḡ, it, x = solve(A₄̂, f₄̄)

# ╔═╡ abc4a380-2b37-47a4-afed-f3122f142827
B̂

# ╔═╡ ee3a783f-b20b-4fdb-8752-9f3ca17462fa
∞(B̂)

# ╔═╡ d22af669-5481-4689-b68d-41db97c3d8c4
∞(A₄̂ * x[end] - f₄̄)

# ╔═╡ f12203e5-8d90-456d-8bc9-afa87880e920
([
	2 1 α
	α 5 0.72
	-1.2 3 1.7
] * x[end] - [-2.9, -0.7, -9.86]) |> ∞

# ╔═╡ 505969fb-57e7-4a92-9e27-2668f6b8d7ca
A₄̂ \ f₄̄

# ╔═╡ 05e30414-3e2b-4642-b685-5d03aa5b4b20
md"""
## Задание 2
"""

# ╔═╡ a5e5afd1-b39c-49bb-a7eb-2c8c8a51e9c3
md"""
$B = \begin{bmatrix}
	p && q \\
	q && p
\end{bmatrix}$
"""

# ╔═╡ d5ee6675-2479-4f3e-9ba0-5240548ca69f
NORM1(x, y) = ∞([x y; y x])

# ╔═╡ 3fad14be-8441-43d6-b327-048cae944656
plot(
	NORM1 ≪ 1;
	xlims=(-2, 2),
	ylims=(-2, 2)
)

# ╔═╡ 11d3e4c6-6206-4b3b-a9d4-6c002cfe0953
md"""
$B = \begin{bmatrix}
	p && q \\
	-q && p
\end{bmatrix}$
"""

# ╔═╡ 1eafd2f4-9fa2-4dff-9bf3-031e11e907a6
NORM2(x, y) = ∞([x y; -y x])

# ╔═╡ 936fb5fd-b611-4998-a813-0b17ceb39251
plot(
	NORM2 ≪ 1;
	xlims=(-2, 2),
	ylims=(-2, 2)
)

# ╔═╡ Cell order:
# ╠═0cb92a82-8e93-11ec-091b-f1abb82c2225
# ╟─32f02311-d88b-47ef-a399-59752ce037a5
# ╠═2fd3c08f-0768-4533-90be-d0f457f77bdf
# ╠═90464d1e-75aa-423c-9f29-dce94cf84003
# ╠═8d82ca3e-b058-49b7-8398-064de7da35d3
# ╟─f3572ec8-4ec9-4b87-b39b-34f3f652c8a0
# ╟─f2797a5d-8f42-4c03-ac10-4331c1075ff1
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
# ╟─f5271f10-5515-46e6-ae62-cec7eec1b49e
# ╠═2865dc6f-5c92-4752-88b8-5edb3e261e39
# ╠═1db51949-df20-44eb-bbec-9047596a1284
# ╠═9f07676a-1777-42d5-aa04-68d1e3aebf86
# ╟─690da59c-dc1f-415f-89c4-9246216f1f41
# ╠═b8ad41d5-88ee-4f75-9c54-5578501edf16
# ╠═d0e2df4d-df6a-4b2e-9843-58e9fcab2825
# ╠═41776634-e2f1-427a-b516-ebd69d2293c4
# ╠═8c76b049-469c-4663-9b93-af04e1821e92
# ╠═76d93737-1712-4cdc-ab45-95c70eb752b1
# ╠═0d33adbc-2cc6-4333-9765-ca5a6877a9a9
# ╟─65ec3210-63e4-45c6-bef6-e1bb45d3d16f
# ╠═377c0a66-d2f7-445d-8425-ea5fa0d11bce
# ╠═38701eb2-8750-4d3f-add0-3d7f0c23d1f5
# ╠═5bfe6216-29fb-41be-ba7a-22b8cbd8d1cb
# ╠═94e4d42c-8ca3-4488-9959-bffe12088477
# ╠═a4b44707-15b9-4533-befe-4681f011c665
# ╠═35928c65-a4b5-4c9f-8c05-886fa514a061
# ╟─b5e0fa99-40de-4ee9-983d-ce4dce9e6dcc
# ╠═8d7419f9-2885-426e-9709-fdceac7537e5
# ╠═dbf0aaf2-9132-4ab4-b2ce-3ae85b1986f5
# ╠═50da520d-5de8-4565-b411-9ad2de8a3131
# ╠═5d374e90-202b-4d69-8e4a-f10328af3b30
# ╠═6ebc8f74-254f-4211-a34a-9e1f1b2dd93e
# ╠═850a3b72-1f63-463e-8833-9256062711ff
# ╠═36446fa9-2848-4147-822a-81ff4acb16e9
# ╠═6d1fc1d3-e69e-47d1-933a-c14e74e52496
# ╠═f507f669-e56c-4588-80bb-c47e4076c1c4
# ╠═8c6a9db4-5509-47b2-8cf1-7857ca1dd6ed
# ╟─2bad0e5c-b1a7-47c1-b7e4-5824c4cb331c
# ╠═95cd2ca2-4e51-489a-b184-a726c7080b01
# ╠═be9cdc8a-f563-4d5f-9037-cb7c13d5dd13
# ╠═1e232704-cbd9-4d4f-9731-7bf7d04c37ee
# ╠═e48e96e2-ce80-436f-b065-6abf7c18eeff
# ╠═4f2181c2-ccc9-42e9-b50f-c3a4f5d44267
# ╠═890060d2-4307-4c0f-b2c2-d5711c1e8cc1
# ╟─6c4ccf17-99cb-4dab-89a8-3cdf10c049ea
# ╠═b166cf1b-d5d7-4b0d-b4ac-54967d3cd8a9
# ╠═95fa8cbe-9c17-4277-a17a-a1299378a368
# ╠═159b490d-a0ee-4094-8bbf-4f5ed6671e6e
# ╠═95aa6b37-dfcc-489f-8d9e-29dd424bbdf0
# ╠═f667d18f-a070-4921-9f6e-40a05abb13dc
# ╟─9afaaf5d-e1a3-46f6-98dd-3e66f5bfa7d1
# ╠═9615683c-bb4c-4d42-bbe6-5044ab2bedc7
# ╠═63083def-2687-4635-bf08-d82935c62b1f
# ╠═dfa13584-a688-4ad4-b1d9-f8cdcc9c4727
# ╠═9020a0ba-4320-45b8-b1ae-6dc256929528
# ╠═abc4a380-2b37-47a4-afed-f3122f142827
# ╠═ee3a783f-b20b-4fdb-8752-9f3ca17462fa
# ╠═d22af669-5481-4689-b68d-41db97c3d8c4
# ╠═f12203e5-8d90-456d-8bc9-afa87880e920
# ╠═505969fb-57e7-4a92-9e27-2668f6b8d7ca
# ╟─05e30414-3e2b-4642-b685-5d03aa5b4b20
# ╠═6d05ed77-1924-4617-978b-c1fbfc0de550
# ╟─a5e5afd1-b39c-49bb-a7eb-2c8c8a51e9c3
# ╠═d5ee6675-2479-4f3e-9ba0-5240548ca69f
# ╠═3fad14be-8441-43d6-b327-048cae944656
# ╟─11d3e4c6-6206-4b3b-a9d4-6c002cfe0953
# ╠═1eafd2f4-9fa2-4dff-9bf3-031e11e907a6
# ╠═936fb5fd-b611-4998-a813-0b17ceb39251
