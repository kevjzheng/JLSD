### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 8f2e1bc4-597d-4e02-b406-0a3908555350
using BenchmarkTools

# ╔═╡ b54766f6-263a-4fbe-991a-78a8207b5a0c
md"""
# Building SerDes Models in Julia, pt.1 - Background

![Cover UI](https://circuit-artists.com/wp-content/uploads/2024/06/cover_new.png)

This notebook shows how to build and simulate a *"simple"* SerDes model using Julia. It serves as a documentation of my own journey of learning Julia so far. Hopefully it can turn into something bigger in the future (i.e., a more generic and expandable SerDes simulation framework). There will be several parts

1. Background and getting started on Julia
2. SerDes simulation framework
2. Data structures, math, and more
3. Detailed transmitter example
4. Plotting (w/ Makie)

Due to limited time and resource, the entire source code for this project won't be available yet on GitHub, but enough will be explicitly shown in the notebook to demonstrate the framework (which I believe is the most important part) and what Julia can do
"""

# ╔═╡ 54a9ae9d-888f-4fc0-945c-fb4e5534c9e1
md"""
## Why Julia?

While there are many languages already available for scientific computing, like Python and MATLAB, no one is without its fault. Below are some discussions on their pros/cons, and why I believe Julia is worth investigating

"""

# ╔═╡ 8df26aaa-cf88-42ff-bbf8-50d445e5f858
md"""
### Python	
Despite its huge popularity, Python is almost a non-starter for me when it comes to buidling and simulating SerDes models. Python is a dynamic and interpreted language, which means the code is only converted to machine code at run time. Python's syntax is very human friendly and easy to read/learn, but this contributes to it being too slow compared to compiled languages like C/C++. Python is a great prototyping/scripting language, and when speed/performance doesn't matter as much, it's SWEET!

For computation heavy tasks, [Numpy](https://numpy.org/) hides away Python's slowness by essentially pre-compiling all array data structures and operations in C code. However, speed of your model can then be bottlenecked by some new algorithm you are experimenting with, an explicit and necessary for-loop somewhere, or post-sim data analysis. 

There are existing frameworks for SerDes modeling in Python, like [SerDesPy](https://github.com/richard259/serdespy) and [SymbaPy](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9293846) (not available on GitHub anymore). However, speed isn't the strongest point for these packages. SymbaPy quotes ~200 seconds for simulating 1 million bits (and ~2600 seconds w/ eye diagrams). 

Nevertheless, Python's ease of entry and open-source community remains key in its wide adoption among academia and industry experts alike. Therefore, it's definitely not going away any time soon. 🐍🐍🐍🐍🐍

"""

# ╔═╡ da497225-67b5-483b-b4f2-7cf5a4eac3e9
md"""
### MATLAB 
I love MATLAB for its plotting and tooling features (open source packages are still far away from what MATLAB plots are capable of). Its revamped Just-in-Time (JIT) compiler since R2015 gives a good trade-off between code performance and feeling "snappy" (like an interpreted language). MATLAB remains the go-to choice for industry and scientific communities due to its one-size-fit-all suite. 

But well... it's not free (and I don't blame them for having so many great functions and features out of the box).

However, because of its closed source nature, the langauge made conscious decisions to eliminate concepts like "pointers" in C/C++, and manage variable passing internally. This could lead to memory inefficiency and allocation issues if the programmer is not careful. In the case of simulating millions of bits in a SerDes system, I often found MATLAB to be a big memory hog and simulation could begin to slow down as time goes on (due to [garbage collection](https://en.wikipedia.org/wiki/Garbage_collection_(computer_science)))

Simulink, SerDes Toolbox, IBIS-AMI compatibility, etc. are other big pluses that come with MATLAB. Yet, most simulation models are still custom built to evaluate proprietary architecture and algorithm. Besides, codes are often compiled again to speed up simulation. And finally, for someone with a soft spot for open-source, I find more joy in learning a new language then figuring out how to using new proprietary tools. 😛

"""

# ╔═╡ 20e69978-795e-402c-98bc-ec69fd5c677e
md"""
### Julia
I came across [Julia](https://julialang.org/) because I googled "open-source alternative to MATLAB" one day like I have done many times before. 
![Julia overview](https://circuit-artists.com/wp-content/uploads/2024/06/julia_overview.png)

Essentially, Julia aims to solve the "two language problem" - we want a langauge that is easy/flexible like Python, but fast like C (especially for scientific computing). Julia's own introduction paragraph summarizes itself well:

!!! info "Julia Compared to Other Languages"
	Julia features optional typing, multiple dispatch, and good performance, achieved using type inference and just-in-time (JIT) compilation (and optional ahead-of-time compilation), implemented using LLVM. It is multi-paradigm, combining features of imperative, functional, and object-oriented programming. Julia provides ease and expressiveness for high-level numerical computing, in the same way as languages such as R, MATLAB, and Python, but also supports general programming. To achieve this, Julia builds upon the lineage of mathematical programming languages, but also borrows much from popular dynamic languages, including Lisp, Perl, Python, Lua, and Ruby. [link here](https://docs.julialang.org/en/v1/#man-julia-compared-other-languages)

I then decided to migrate some of my MATLAB models to Julia and do a benchmark comparison. The porting is relatively straight-forward once I picked up Julia's typing system, and below are code snippets of a PRBS generator function in Julia and MATLAB

"""

# ╔═╡ fc826d03-0652-437e-9e8e-ad515c585f42
md"""
This is a Julia function
```julia
function bist_prbs_gen(;poly, inv, Nsym, seed)
    seq = Vector{Bool}(undef,Nsym)
    for n = 1:Nsym
        seq[n] = inv
        for p in poly
            seq[n] ⊻= seed[p]
        end
        seed .= [seq[n]; seed[1:end-1]]
    end
    return seq, seed
end
```
"""

# ╔═╡ 088317b5-4718-49d8-8c9f-e171631d7f31
md"""
This is a MATLAB function
```julia
function [seq, seed] = bist_prbs_gen(poly,inv, Nsym, seed)
    seq = zeros(1,Nsym);
    for n = 1:Nsym
        seq(n) = inv;
        for p = poly
            seq(n) = xor(seq(n), seed(p));
        end
       seed = [seq(n), seed(1:end-1)];
    end
end
```

"""

# ╔═╡ 1bd48b0b-3a48-42dc-85db-7dbdcf3d5587
md"""
You can see the syntax between Julia and MATLAB is quite similar. The algorithm here is very straight-forward. Admittely, I didn't do much optimization on the register and bit level as a first pass, but the goal here is to really test the capability of the compilers.

A couple of things to highlight
1. Julia allows special unicode characters as part of the code! The symbol ⊻ means "xor". Though it doesn't make a difference in the code functionality, it does make one happy when you can use π and ℯ in your code.
2. There are some smaller nuiances when dealing with Julia arrays. The .= syntax is not a typo, but rather a broadcasting operation similar to MATLAB's element-wise operations (more on this in the future).
3. Variable typing in Julia is optional. Like Python, Julia can dynamically type variables at runtime. However, for the more advanced users, defining variable types early in Julia is also possible and often helps with performance.

Now, let's talk about speed. If I benchmark the functions by generating and timing 10 million bits using PRBS31 (i.e. bist\_prbs\_gen(polynomial = [28,31], inv = false, Nsym = 10e6, seed = ones(31)), below are the results. I needed to make sure I didn't make a mistake, but yes there is an almost **10x** difference between these seemingly identical functions!
![julia_vs_matlab](https://circuit-artists.com/wp-content/uploads/2024/06/julia_v_matlab_10x.png)


"""

# ╔═╡ 9bbe4518-cf4b-4125-8e39-e64fedd3fe19
md"""
For those interested, you can try implementing the same in Python. Surprisingly, the implementation using Python list instead of Numpy array is faster when I tried (~6s for list and ~12s for numpy arrays).

Assuming I didn't do anything blatantly false, this just adds to the confusion about how to optimize Python code when performance matters.

"""

# ╔═╡ 2b4fb604-a29d-43bf-bc1f-2fba4d4dea07
md"""
```python
import time
import numpy as np

def bist_prbs_gen(poly, inv, nsym, seed):
    seq = [False]*nsym
    for n in range(nsym):
        seq[n] = inv
        for p in poly:
            seq[n] ^= seed[p-1]
        seed[1:] = seed[0:-1]
        seed[0] = seq[n]

    return seq, seed


def bist_prbs_gen_arr(poly, inv, nsym, seed):
    seq = np.array([False]*nsym)
    for n in range(nsym):
        seq[n] = inv
        for p in poly:
            seq[n] ^= seed[p-1]
        seed[1:] = seed[0:-1]
        seed[0] = seq[n]

    return seq, seed



if __name__ == '__main__':
    start = time.time()
    seq1, seed1 = bist_prbs_gen([28, 31], False, int(10e6), [True]*31)
    end = time.time()
    print(end - start)   ## Showed 5.537s

    start = time.time()
    seq2, seed2 = bist_prbs_gen_arr([28, 31], False, int(10e6), np.array([True]*31))
    end = time.time()
    print(end - start)   ## Showed 13.543s
```
"""

# ╔═╡ e931847b-876b-4a91-b3cd-28ef6e45bf24
md"""
To drive the speed point home for Julia, the plot at the very top was generated after processing **one million bits** using a model that includes a jittered transmitter,  noisy channel, a simple RX with slicer adaptation and CDR loops. **The entire simulation only took ~10 seconds**. 

Of course, this example might be a bit artificial, but the **order of magnitude** speed improvement alone was enough to push me to spend more time learning Julia, so here we are 😄.
"""

# ╔═╡ 9305392e-d630-48ed-98ed-9ee1cc52b951
md"""
# Getting Started with Julia

Julia is relatively young compared to Python, so development tools might look much more "basic". However, it still offers an easy installation and a good extension in [Visual Studio](https://visualstudio.microsoft.com/) (my current IDE choice).

If you want to go through my notebooks or you are intrigued by what Julia can do for your field, follow the steps [**here**](https://computationalthinking.mit.edu/Spring21/installation/) to install Julia and Pluto (VS Code is optional if you want to expand upon what are shown in the notebooks)

"""

# ╔═╡ ad2a66fa-2e5a-4953-99bf-fb4369781d14
md"""
Once you can run Pluto from the Julia REPL like below, open this notebook using the GitHub URL (or copy the .jl file to your local director and open). 

![Julia REPL run Pluto](https://user-images.githubusercontent.com/6933510/91441094-eb01d080-e86f-11ea-856f-e667fdd9b85c.png)
![Pluto homepage](https://user-images.githubusercontent.com/6933510/91441391-6a8f9f80-e870-11ea-94d0-4ef91b4e2242.png)
"""

# ╔═╡ 1693bc36-3a9a-4947-9bbd-99b057672a62
md"""
If everything goes well, you should be able to run the small *bist\_prbs\_gen* benchmark next.

The BenchmarTools is the first helpful package we will learn to use. A macro is part of Julia's metaprogramming capability that simplies code (think of it as a function manipulating the string of your code), denoted by the @ sign. The @time macro is used to benchmark the execution time. What's neat is it also tells how many memory allocations are done, as well as compilation and garbage collection time.


Once you are used to it, it forces you to think a bit deeper down to the machine level, but you still remain writing in a high level language. **VERY SATISFYING!** 😍

"""

# ╔═╡ 086c7f93-17fb-4a74-a398-0a85f8afad52
function bist_prbs_gen(;poly, inv, Nsym, seed)
    seq = Vector{Bool}(undef,Nsym)
    for n = 1:Nsym
        for p in poly
            seq[n] ⊻= seed[p]
        end
        seed .= [seq[n]; seed[1:end-1]]
    end
    return seq.⊻inv, seed
end

# ╔═╡ 6f0e5fb3-0644-489e-809f-691ee56d7e2a
@time bist_prbs_gen(poly=[28,31], inv=false, Nsym=Int(10e6), seed=ones(Bool,31));

# ╔═╡ 86fc629f-f925-43e0-ae21-98593eea1833
md"""
For those who are more familiar with Python, Julia's JIT compilation might feel weird at first because the run time will be longer for the first run (due to compilation time) and you get the real time for subsequent runs. It's something you need to get used to if you are accustomed to Python's consistent behavior after hitting run. Nevertheless, this in my opinion makes Julia more suitable for running long simulations with more complicated models rather than simple scripting on lab benches.
"""

# ╔═╡ 69ca90d4-efc8-4555-839e-a4756b766fac
md"""
Once you are setup with Julia, I recommend introducing yourself to Julia through the following resources, and familiarize with these packages

- ["Why We Created Julia"](https://julialang.org/blog/2012/02/why-we-created-julia/)
- [Julia's official documentation](https://docs.julialang.org/en/v1/)
- [BenchmarkTools.jl](https://juliaci.github.io/BenchmarkTools.jl/stable/) - performance tracking
- [Revise.jl](https://timholy.github.io/Revise.jl/stable/) - makes Julia feel more like a runtime language when code changes
- [UnPack.jl](https://github.com/mauro3/UnPack.jl) - easier syntax to manipulate struct data
- [DataStructures.jl](https://juliacollections.github.io/DataStructures.jl/dev/) - special data structures like buffers, queues, trees, etc.
- [Parameters.jl](https://mauro3.github.io/Parameters.jl/stable/) - easier model parameter handling with default values and keywords
- [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/v0.14/) - interpolation in Julia
- [DSP.jl](https://docs.juliadsp.org/stable/contents/) - fast DSP functions like convolutions
- [Makie.jl](https://docs.makie.org/v0.21/) - my choice of plotting package in Julia
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"

[compat]
BenchmarkTools = "~1.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "61c30062a384602fe87efb278f2fa7342de39e7e"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╟─b54766f6-263a-4fbe-991a-78a8207b5a0c
# ╟─54a9ae9d-888f-4fc0-945c-fb4e5534c9e1
# ╟─8df26aaa-cf88-42ff-bbf8-50d445e5f858
# ╟─da497225-67b5-483b-b4f2-7cf5a4eac3e9
# ╟─20e69978-795e-402c-98bc-ec69fd5c677e
# ╟─fc826d03-0652-437e-9e8e-ad515c585f42
# ╟─088317b5-4718-49d8-8c9f-e171631d7f31
# ╟─1bd48b0b-3a48-42dc-85db-7dbdcf3d5587
# ╟─9bbe4518-cf4b-4125-8e39-e64fedd3fe19
# ╟─2b4fb604-a29d-43bf-bc1f-2fba4d4dea07
# ╟─e931847b-876b-4a91-b3cd-28ef6e45bf24
# ╟─9305392e-d630-48ed-98ed-9ee1cc52b951
# ╟─ad2a66fa-2e5a-4953-99bf-fb4369781d14
# ╟─1693bc36-3a9a-4947-9bbd-99b057672a62
# ╠═086c7f93-17fb-4a74-a398-0a85f8afad52
# ╠═8f2e1bc4-597d-4e02-b406-0a3908555350
# ╠═6f0e5fb3-0644-489e-809f-691ee56d7e2a
# ╟─86fc629f-f925-43e0-ae21-98593eea1833
# ╟─69ca90d4-efc8-4555-839e-a4756b766fac
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
