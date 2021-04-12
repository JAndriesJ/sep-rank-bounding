# ju-sep-rank

##  Code Overview:
  - ***Utils.jl***:
    - sd
  - ***Examples.jl*** :
    - Returns inbuilt examples for which the separable rank was computed.
  - ***Moments.jl*** :
    - Various...
  - ***sep_Constraints.jl*** :
    - The constriants of the optimization problem.
  - ***sep_Compute.jl*** :
    - The computation script


## Getting started:



## Theory:
The goal is to compute a lower bound for the separable rank of a density matrix representing a quantum mechanical state.
#### Sources:
- 1. [Qubit-qudit states with positive partial transpose](http://arxiv.org/abs/1210.0111v2)
- 2. [Separability of Hermitian Tensors and PSD Decompositions](https://arxiv.org/abs/2011.08132v1)
- 3. [Separability of mixed states: necessary and sufficient conditions](https://arxiv.org/pdf/quant-ph/9605038v2.pdf)

### Constraints
#### S₁
```
L ≥ 0 on M₂ₜ(S_ρ¹)
⟺
L(g⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0 for
g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d]
```

#### S₂
```
L ≥ 0 on M₂ₜ(S_ρ²)
⟺
L(g⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0 for
g = √(Tr(ρ)) - ∑xᵢ² , √(Tr(ρ)) - ∑yᵢ²
```

#### S₃
```
L ≥ 0 on M₂ₜ(S_ρ³)
⟺
L((Tr(ρ) - ∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
L((sqrt(∑yᵢ²) - 1)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) = 0
```

#### wG:
```
L((ρ - ([x,y]₌₁[x,y]₌₁ᵀ)) ⊗ ([x,y]₌ₗ[x,y]₌ₗᵀ)))for l ∈ 1,...,t-2.
= ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L(([x,y]₌₁[x,y]₌₁ᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))
```
#### sG:
```
M(Gρ ⊗ L) ⪰ 0 constraints
Where: Gρ := ρ - xxᵀ⊗yyᵀ
M(Gρ ⊗ L) = L(Gρ ⊗ [x, y]ₜ₋₂[x, y]ᵀₜ₋₂) ⪰ 0
output: ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
```

## Example States:
### Separable states:
#### Example sep1
|00><00| + |11><11|
```
 1  0  0  0
 0  0  0  0
 0  0  0  0
 0  0  0  1
```

#### Example sep2
|00><00| + |11><11| + (|0> + |1>)⊗(|0> + |1>)(<0| + <1|)⊗(<0| + <1|)
```
 2  1  1  1
 1  1  1  1
 1  1  1  1
 1  1  1  2
```

#### Example sep3
|00><00| + |11><11| + |12><12|
```
 1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  1  0  0  0  0
 0  0  0  0  0  1  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
```


#### Example sep4
|00><00| + |01><01| + |11><11| + |12><12|
```
 1  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  1  0  0  0  0
 0  0  0  0  0  1  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
```

#### Example sep5
|00><00| + |01><01| + |02><02| + |11><11| + |12><12|
```
 1  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0
 0  0  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  1  0  0  0  0
 0  0  0  0  0  1  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
```


#### Example sep6
Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁j₁ + i₂j₂ where i₁,j₁,i₂,j₂ in [3]
```
 2  3   4  3   4   5   4   5   6
 3  5   7  4   6   8   5   7   9
 4  7  10  5   8  11   6   9  12
 3  4   5  5   6   7   7   8   9
 4  6   8  6   8  10   8  10  12
 5  8  11  7  10  13   9  12  15
 4  5   6  7   8   9  10  11  12
 5  7   9  8  10  12  11  13  15
 6  9  12  9  12  15  12  15  18
```


### Entangled states:


#### Example ent1
(|01> + |10>)(<01| + <10|)
```
 1  0  0  1
 0  0  0  0
 0  0  0  0
 1  0  0  1
```



#### Example ent2
```
a  = 0.5; b = 0.8; # arbitary possitive numbers
p  = rand()
ψ₁ = a*ψ(2,2,2) + b*ψ(1,1,2)
ψ₂ = a*ψ(2,1,2) + b*ψ(1,2,2)
ρ[2,"ent2"] = p*sq(ψ₁) + (1 - p)*sq(ψ₂)

 0.479379  0.0       0.0        0.299612
 0.0       0.160621  0.100388   0.0
 0.0       0.100388  0.0627425  0.0
 0.299612  0.0       0.0        0.187258
```


#### Example ent3
```
a =  1/sqrt(2) ; p = rand();
ψₐ = a*(ψ(1,2,2) - ψ(2,1,2))
p*sq(ψₐ) + (1 - p)*psep(1,1,2)

 0.426205   0.0        0.0       0.0
 0.0        0.286897  -0.286897  0.0
 0.0       -0.286897   0.286897  0.0
 0.0        0.0        0.0       0.0
```

#### Example ent4
|00><00| + |02><02| + 2|11><11| + (|01> + |10>)(<01| + <10|)
```
 1  0  0  0  0  0  0  0  0
 0  1  0  1  0  0  0  0  0
 0  0  1  0  0  0  0  0  0
 0  1  0  1  0  0  0  0  0
 0  0  0  0  2  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
```


#### Example ent5  (untested)
Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁i₂ + j₁j₂ where i₁,j₁,i₂,j₂ in
```
  2   3   4   3   5   7   4   7  10
  3   4   5   4   6   8   5   8  11
  4   5   6   5   7   9   6   9  12
  3   4   5   4   6   8   5   8  11
  5   6   7   6   8  10   7  10  13
  7   8   9   8  10  12   9  12  15
  4   5   6   5   7   9   6   9  12
  7   8   9   8  10  12   9  12  15
 10  11  12  11  13  15  12  15  18
```
