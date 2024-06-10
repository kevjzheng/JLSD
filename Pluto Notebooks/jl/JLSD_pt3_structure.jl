### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ e694af8f-8abf-408b-9df1-074a5c600256
using Parameters

# ╔═╡ 297d11ed-1eea-4c65-ab62-7773dbb8f103
module TrxStruct
using Parameters, DataStructures
	
@kwdef mutable struct Param
	const osr::Int64 			# needs initialization
	const data_rate::Float64 	# needs initialization
	const pam::Int8 = 2	
	const bits_per_sym::Int8 = Int(log2(pam))
	const fbaud::Float64 = data_rate/bits_per_sym
	const fnyq::Float64 = fbaud/2      
	const tui::Float64 = 1/fbaud			
	const dt::Float64 = tui/osr		

	const nsym_total::Int64 	# needs initialization
	const blk_size::Int64		# needs initialization
	const blk_size_osr::Int64 = blk_size*osr
	const nblk::Int64 = Int(round(nsym_total/blk_size))
	const subblk_size::Int64 = blk_size 
	const subblk_size_osr::Int64 = subblk_size*osr
	const nsubblk::Int64 = Int(blk_size/subblk_size)
	
	cur_blk::Int = 0
	# ... and more
end

@kwdef mutable struct Bist
	const param::Param

	#shared
	polynomial::Vector{UInt8}
	order::UInt8 = maximum(polynomial)
	inv = false

	#generator
	gen_seed = ones(Bool, order)
	gen_gray_map::Vector{UInt8} = []

	#checker
	chk_seed = zeros(Bool, order)
	chk_gray_map::Vector{UInt8} = gen_gray_map
	chk_lock_status = false
	chk_lock_cnt::Int64 = 0
	const chk_lock_cnt_threshold::Int64 = 256
	ber_err_cnt::Int64 = 0
	ber_bit_cnt::Int64 = 0

	#input/output vectors
	So_bits::Vector = zeros(Bool, param.bits_per_sym*param.blk_size)
	So::Vector{Float64} = zeros(param.blk_size)
	Si = CircularBuffer{UInt8}(param.blk_size) 
	#CircularBuffer is a special data type from the DataStructures package
end

@kwdef mutable struct Drv #dummy struct for TX driver
	param::Param
end

@kwdef mutable struct Ch #dummy struct for channel
	param::Param
end

@kwdef mutable struct Cdr #dummy struct for channel
	param::Param
end

#more structs definitions here
end

# ╔═╡ c7b78049-1cef-46d7-aec5-de8775cb29a9
module BlkBIST
using UnPack

export pam_gen_top!, ber_checker_top!

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


function pam_gen_top!(bist)
    @unpack pam, bits_per_sym, blk_size = bist.param
    @unpack polynomial, inv, gen_seed, gen_gray_map = bist
	@unpack So, So_bits = bist

	#generate PRBS bits
    So_bits, gen_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                	Nsym=bits_per_sym*blk_size, seed=gen_seed)

	#generate PAM symbols
    fill!(So, zero(Float64)) #reset So to all 0
    for n = 1:bits_per_sym
        @. So = So + 2^(bits_per_sym-n)*So_bits[n:bits_per_sym:end]
    end

    #gray encoding
    if ~isempty(gen_gray_map)
        for n in 1:blk_size
            So[n] = gen_gray_map[So[n] + 1]
        end
    end

    @. So = 2/(pam-1)*So - 1 #convert to analog voltage levels in +/-1 range

	return nothing
end

function ber_checker_top!(bist)
    @unpack cur_blk, pam, bits_per_sym = bist.param
    @unpack gen_gray_map, chk_start_blk, Si, Si_bits = bist 


    if cur_blk >= chk_start_blk #make start blk a parameter later
        if ~isempty(gen_gray_map)
            for n in 1:blk_size
                Si[n] = gen_gray_map[Si[n] + 1]
            end
        end
    
        Si_bits .= vec(stack(int2bits.(Si, bits_per_sym)))

        ber_check_prbs!(bist)
    end

	return nothing
end


function ber_check_prbs!(bist)
    @unpack polynomial, inv, chk_seed, Si_bits = bist
    nbits_rcvd = lastindex(Si_bits)

	# Uncomment if you want to add artifical BER
    # err_loc = rand(Uniform(0,1.0), nbits_rcvd).< 1e-4;
    # Si_bits .= Si_bits .⊻ err_loc

    if bist.chk_lock_status
        ref_bits, chk_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                                 Nsym=nbits_rcvd,seed=chk_seed)

        bist.ber_err_cnt += sum(Si_bits .⊻ ref_bits)
        bist.ber_bit_cnt += nbits_rcvd
    else
        for n = 1:nbits_rcvd
            brcv = Si_bits[n]
            btst = inv
			for p in polynomial
            	btst ⊻= chk_seed[p]
			end

            #need consecutive non-error for lock. reset when error happens
            bist.chk_lock_cnt = (btst == brcv) ? bist.chk_lock_cnt+1 : 0

            chk_seed .= [brcv; chk_seed[1:end-1]]

            if bist.chk_lock_cnt == bist.chk_lock_cnt_threshold
                bist.chk_lock_status = true
                println("prbs locked")
                ref_bits, chk_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                                        Nsym=nbits_rcvd-n, seed=chk_seed)
                bist.ber_err_cnt += sum(Si_bits[n+1:end] .⊻ ref_bits)
                bist.ber_bit_cnt += nbits_rcvd-n
                break
            end
        end
    end

	return nothing
end

function int2bits(num, nbit)
    return [Bool((num>>k)%2) for k in nbit-1:-1:0]
end

end

# ╔═╡ d2719225-8dd5-4385-b9d6-4e0f0d577a4e
using Makie, CairoMakie

# ╔═╡ 5cacc554-6c16-4f3e-a0bf-6acd8db7b1e7
using UnPack, DataStructures, Random, DSP

# ╔═╡ e6ef929a-1c1f-47c9-a1c2-888b84b641e4
using BenchmarkTools

# ╔═╡ 79fd0370-2446-11ef-2922-21820062726d
md"""
# Building SerDes Models in Julia, pt3 - Data and Code Structures
"""

# ╔═╡ c34485a7-3583-4c48-945f-f94435f90257
md"""
![code structure](https://circuit-artists.com/wp-content/uploads/2024/06/data_code_structure.png)
"""

# ╔═╡ 417676db-f650-4d85-84dc-465ee504a03d
md"""
Julia is not an Object Oriented Programming (OOP) language. Unlike Python where Class and Object definitions are supported, Julia's main focus is dealing with data (closer to MATLAB in this aspect even though MATLAB has OOP support too). OOP typically is good at managing and scaling complex software systems, but usually at the cost of performance (here is a famous [video rant](https://www.youtube.com/watch?v=QM1iUe6IofM&ab_channel=BrianWill) on OOP)
"""

# ╔═╡ 0d9aa8be-bf0d-41b2-91b5-d702195a3d18
md"""
In reality, OOP is just a programming paradigm that focuses on encapsulating *states* and handling (im)mutability at a large scale. For our purpose, we do have quite some states to take care of for each circuit module (e.g. seed in BIST, channel memory/ISI, CDR accumulators, and even jitter), so how can we achieve this without OOP in Julia?
"""

# ╔═╡ 7d46974f-5531-4408-a96f-4672385665e1
md"""
## Struct - packing parameters, data and states
"""

# ╔═╡ bb74fb24-4fb0-4fc6-b2c6-e9bf77f3d908
md"""
Enters **```struct```** - a [composite data type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) in Julia. Think of it as a collection of any data that you want to throw at it, or object with only *attributes/properties*. MATLAB also has ```struct```, and the closest native data type in Python might be dictionary (though they are still quite different). Here is a Julia struct definition:
"""

# ╔═╡ df1878ff-130d-4db4-8803-2772237a8505
struct Params_im
	osr::Int64 			#oversampling ratio for simulator
	data_rate::Float64 	#nominal data rate
	pam::Int8 			#number of data levels
	blk_size::Int64		#number of symbol in a block
	subblk_size::Int64  #number of symbols in a sub-block
	fbaud::Float64      #baud rate
	fnyq::Float64       #nyquist frequency
	tui::Float64 		#unit interval time
	dt::Float64 		#simulation time step
	# ... and more
end

# ╔═╡ f734a5c7-d0aa-4f65-b967-55bf859ad81e
p_im = Params_im(32, 10e9, 2, 2^10, 16, 10e9, 10e9/2, 1/10e9, 1/10e9/32);

# ╔═╡ 55db57f4-4918-4493-a71b-67cc4af609d6
md"""
!!! info "Julia Tips"
	Each field can be declared with a specific type (like Float64) or left empty (to be determined at runtime, but often less performant). Try giving the compiler as much as information for the code to be optimized, especially when for arrays/matrices
"""

# ╔═╡ 6466aa8f-c5bb-4643-b72c-ff687d992744
md"""
Just like that, we have successfully instantiated our first struct, packing all the global simulation parameters and convience variables together. Each field can be accessed with the "." operator.
"""

# ╔═╡ 5372f844-cedf-45e8-8799-106372630c10
p_im.data_rate #access field

# ╔═╡ a9ab933a-20a8-4948-a662-d5338b697034
md"""
However, struct in Julia is **immutable** by default - you can only read and not modify its internal fields once instantiated
"""

# ╔═╡ 282f989a-e9d5-45d2-a044-3d87df1c4827
p_im.osr = 16 #this will throw an error!

# ╔═╡ 04002f49-604b-42b1-9b1c-21b5bd014240
md"""
The reason is that immutability lets compilers use memory efficiently. But don't worry, Julia have a mutable struct! (note: they will be "slower" than immutable structs)
"""

# ╔═╡ 3ba76818-18a3-41cd-be27-fda0a7fede3d
mutable struct Params_m
	osr::Int64 			#oversampling ratio for simulator
	data_rate::Float64 	#nominal data rate
	pam::Int8 			#number of data levels
	blk_size::Int64		#number of symbol in a block
	subblk_size::Int64  #number of symbols in a sub-block
	fbaud::Float64      #baud rate
	fnyq::Float64       #nyquist frequency
	tui::Float64 		#unit interval time
	dt::Float64 		#simulation time step
	# ... and more
end

# ╔═╡ e65ff729-04ac-4303-bfe4-864f7066b7fe
p_m = Params_m(32, 10e9, 2, 2^10, 16, 10e9, 10e9/2, 1/10e9, 1/10e9/32);

# ╔═╡ 1136c126-c359-4c46-95d5-7d8b92c672ef
p_m.data_rate #so far so good

# ╔═╡ ad83cf45-9b84-45bc-874f-5efbcc3bfa71
p_m.osr = 16 #yay!

# ╔═╡ fcbbbab5-5ed8-498b-b9af-005e1e112af7
md"""
But wait, there are certain things that likely won't change once we set them. Can we have the best of both worlds? Fortunately, Julia now supports immutable fields inside a mutable struct w/ the ```const``` keyword.
"""

# ╔═╡ d812042d-e16d-4535-a358-ce6fdb710933
mutable struct Params_m2
	const osr::Int64 			#oversampling ratio for simulator
	const data_rate::Float64 	#nominal data rate
	const pam::Int8 			#number of data levels
	const blk_size::Int64		#number of symbol in a block
	const subblk_size::Int64  	#number of symbols in a sub-block
	const fbaud::Float64      	#baud rate
	const fnyq::Float64       	#nyquist frequency
	const tui::Float64 			#unit interval time
	const dt::Float64 			#simulation time step
	cur_blk::Int
	# ... and more
end

# ╔═╡ d9522b57-cfb2-4108-a3da-846153141b12
p_m2 = Params_m2(32, 10e9, 2, 2^10, 16, 10e9, 10e9/2, 1/10e9, 1/10e9/32, 0);

# ╔═╡ 5d934cc6-0049-48e5-b4c2-92a7c7ae09fb
p_m2.osr = 16 #oops error again

# ╔═╡ 9ec3d251-4ac3-4262-a1af-3f2b057bf155
p_m2.cur_blk = 1 #this field is mutable

# ╔═╡ 9335e8fb-9fde-414b-8c7c-b618b3b308e5
md"""
Cool! This is all pretty intuitive, but something feels off. When we instantiated a struct instant, the arguments need to be in order and w/o keywords. If the struct becomes bigger, the code isn't exactly easy to read and maintain. Also, some variables are dependent on some other variables (e.g. dt is calculated from tui and osr). A default calculation would be nice. The @kwdef macro from the Parameters.jl package solves both problems.
"""

# ╔═╡ 4755226e-d845-4f70-ae0e-fc2329d0c4c0
@kwdef mutable struct Param
	const osr::Int64 			# needs initialization
	const data_rate::Float64 	# needs initialization
	const pam::Int8 = 2	
	const bits_per_sym::Int8 = Int(log2(pam))
	const fbaud::Float64 = data_rate/bits_per_sym
	const fnyq::Float64 = fbaud/2      
	const tui::Float64 = 1/fbaud			
	const dt::Float64 = tui/osr		

	const nsym_total::Int64 	# needs initialization
	const blk_size::Int64		# needs initialization
	const blk_size_osr::Int64 = blk_size*osr
	const nblk::Int64 = Int(round(nsym_total/blk_size))
	const subblk_size::Int64 = blk_size 
	const subblk_size_osr::Int64 = subblk_size*osr
	const nsubblk::Int64 = Int(blk_size/subblk_size)
	
	cur_blk::Int = 0
	# ... and more
end

# ╔═╡ 51e571fd-367a-41bd-ab56-325b0bc626c5
md"""
```@kwdef``` allows fields to have default values, and the struct constructor becomes keyword based. Now the way to instantiate the struct becomes much self explanatory:
"""

# ╔═╡ 231e7fcc-ecc5-4dd4-b07d-1782413228ad
param = Param(
				osr = 32,
				data_rate = 10e9,
				pam = 4,
				blk_size = 2^10,
				nsym_total = Int(1e6));

# ╔═╡ 934f66ca-6866-485a-aaab-68ab48e01574
param.fbaud #access the instantiated fields here and see their values

# ╔═╡ ecc52900-8d05-4104-a4c4-53f27b899a0f
md"""
Note that the default calculations are done just once during instantiation. That's why using ```const``` is also important to make sure the calculated values remain intact. For example, if ```pam``` field is not constant and instantiated to be 2 and later changed to 4, any dependent field (like ```bits_per_sym```) won't be updated automatically. There is a way to allow mutability and automatic updates, which we will cover in the future, but it's still best practice to keep constants, well, constant. 
"""

# ╔═╡ 50f2f85a-9ec2-40c5-aa36-7ac8646148ac
md"""
Let's now define a Bist struct that contains parameters and states for our PRBS generator/checker in the previous notebook
"""

# ╔═╡ 4e20103f-ef37-4172-9679-1a8ef1f8d4f8
@kwdef mutable struct Bist
	const param::Param

	#shared
	const polynomial::Vector{UInt8}
	const order::UInt8 = maximum(polynomial)
	const inv = false

	#generator
	gen_seed = ones(Bool, order)
	gen_gray_map::Vector{UInt8} = []

	#checker
	chk_seed = zeros(Bool, order)
	chk_gray_map::Vector{UInt8} = gen_gray_map
	chk_lock_status = false
	chk_lock_cnt::Int64 = 0
	const chk_lock_cnt_threshold::Int64 = 256
	ber_err_cnt::Int64 = 0
	ber_bit_cnt::Int64 = 0

	#input/output vectors
	So_bits::Vector = zeros(Bool, param.bits_per_sym*param.blk_size)
	So::Vector{Float64} = zeros(param.blk_size)
	Si = CircularBuffer{UInt8}(param.blk_size) 
	#CircularBuffer is a special data type from the DataStructures package
	Si_bits::Vector = zeros(Bool, param.bits_per_sym*param.blk_size)

end


# ╔═╡ 5d762756-2740-46c7-a45a-62e667d96ef6
bist = Bist(
			param = param, #the param variable defined above
			polynomial = [28,31]); #that's it!

# ╔═╡ 469aa4e7-50f4-488c-883e-6ef8acce7975
md"""
Oh btw, a struct can have a field that points to a struct too. When we instantiate the Bist struct, the ```params``` passed into the Bist field is the reference (or pointer) to the variable. We can check the equivalence of the ```param``` instance and the ```param``` inside ```bist```
"""

# ╔═╡ 68884750-d247-4dec-9019-7e0c179ebea4
param === bist.param #equivalence check

# ╔═╡ 48af9bdd-960b-4e7c-8c62-4b04d6b8a82e
md"""
If something inside ```param``` is changed, it can be accessed through ```bist``` as well
"""

# ╔═╡ 6ce6f7f4-e6d1-4fc5-8391-707279f35696
param.cur_blk = 1;

# ╔═╡ a3128c73-fcad-46d2-bf05-75feb1e3ab16
println(bist.param.cur_blk) 
#Pluto can't tract changes in a struct field, so if you changed params.cur_blk above, remember to rerun this code cell to see the change (like Jupyter)

# ╔═╡ 54117b87-85c5-4757-a736-5dd0bfab8a14
md"""
This is useful when global states need to be passed around, but we don't want our specialized functions (say a prbs_gen function operating only on the ```bist``` object) to take the ```param``` object as an explicit argument too.
"""

# ╔═╡ b52b20fc-7939-430f-bc1b-38e50e4ec1f5
md"""
!!! info "Julia Tips"
	Julia passes mutable variables by their references/pointers into functions, known as "passing by sharing". In contrast, MATLAB passes values of the arguments ("passing by value"), so a NEW copy is made every time the ```param``` instance is passed into another constructor or function (thus leading to memory and speed penalties).

	Julia still creates copies of a immutable variable (like numbers). For more information, [click here](https://www.reddit.com/r/Julia/comments/ll62uu/confusion_about_the_argument_passing_behaviour/)
"""

# ╔═╡ 7c433885-9c9e-48ad-a331-1ee96f3e8119
md"""
## Named Tuples
"""

# ╔═╡ 6d95615e-a9d0-40d6-9c54-e6a158d411f8
md"""
Another useful and flexible data structure in Julia is Named Tuple. It's essentially a tuple with keyword fieldnames but it's immutable. You instantiate with (kw1=val1,kw2=val2...)
"""

# ╔═╡ fa01061f-c472-4ab8-93b1-645d84c8194d
np_test = (a = 1, b=2, c="Foo");

# ╔═╡ f3cf6e75-a0bd-4692-9430-37c39105148f
np_test.a

# ╔═╡ ebf443e0-286d-437d-8bc6-9c52617b39ec
np_test.c = "can't change"

# ╔═╡ b6b8344f-fbb0-4050-8917-824ca7c5545a
md"""
For NP, there is no need to predefine the field names like a struct. So we can just pack different things into a NP at will. Eventually, we will create a big ```trx``` (for transceiver) named tuple that stores all the structs that correspond to each circuit module, like ```bist```, ```drv```, ```ch```, ```cdr```, etc.
"""

# ╔═╡ 5e83caa9-44a2-4f2c-9948-f64eb087cad5
md"""
!!! info "Julia Tips"
	Julia has a shorthand for constructing named tuple using keywords. Any statement after ";" will be expanded into (kw1 = kw1, kw2 = kw2, ...)
	The line above is the same as ```trx = (param=param, bist=bist, drv=drv ...)```. 
"""

# ╔═╡ 7f921c7d-f11e-46f9-99aa-14c55ad082b6
md"""
## Mutating functions
"""

# ╔═╡ be7e4564-630a-4617-bc37-dc270892c8bf
md"""
The simulation framework fully exploits the mutability and "passing by sharing" in Julia, so the input/output vectors are also stored in the struct (i.e., So_bits, So, Si in Bist). Let's define two test functions to see the advantage of mutating functions
"""

# ╔═╡ c306ecd3-f521-420e-868a-2f7ca41b5fa6
function double_return(S)
	return 2*S
end

# ╔═╡ 54c11b28-4550-4603-83e9-cb6188215f3e
function double_mutate!(S)
	@. S = 2 * S
	# the @. macro broadcasts the . to all operations. Equivalent to S .= 2 .* S
	return nothing
end

# ╔═╡ 65c2f05e-b119-4869-89c2-7638c48f2494
S1 = randn(1000000);

# ╔═╡ b59b371a-1e7e-4b9c-86b3-49493de40b34
@time Sx2 = double_return(S1);

# ╔═╡ bdc0c012-541f-4f6b-9504-55e43bd75b8a
S2 = randn(1000000);

# ╔═╡ e5746626-6860-480e-a45d-5640394f12cf
@time double_mutate!(S2)

# ╔═╡ 2cc8af32-7897-46b4-88e1-2e4f40576823
md"""
As we can see, the non-mutating version takes an input, creates and returns a new array with the modified value. It takes 2 allocations for the computation and ~3x longer. That's why you will see many functions with "!" in Julia, like sort!, push!, append!, etc along with their non-mutating counterparts.
"""

# ╔═╡ a7551ac2-c94f-48ff-8bea-9fe1172402e5
md"""
Since our framework simulates on a per-block basis and the vector size for each block is fixed, we will directly operate on the struct's pre-allocated input/output vectors w/ mutating functions. Here are example ```pam_gen_top!``` and ```ber_checker_top!``` functions taking a Bist struct as an input
"""

# ╔═╡ 8e43bdcb-7515-471c-9756-f5bfa8e6ed6d
md"""
## Modules
"""

# ╔═╡ 08035fd1-79d0-424b-b700-01ed1b10a156
md"""
As we discussed in the beginning, struct is similar to a class with only attributes and no member functions. This means that the functions operating on the struct instances need to live somewhere else.
"""

# ╔═╡ 6dd4ab1f-4b80-48a6-a53a-e397b4eb580e
md"""
Julia uses modules to organize code. It's best to go through the official documentation on Modules [here](https://docs.julialang.org/en/v1/manual/modules/). 
We will create two modules to demonstrate how the framework code is structured. (expand the code below to view details)
"""

# ╔═╡ 1a47e550-62b4-4b5d-871b-07b2e3a5d1a0
md"""
!!! info "Julia Tips"
	The ```export``` keyword in a module defines the internal functions that can be directly called at parent level without using the namespace. For example, if BlkBIST is used in the testbench file, ```pam_gen_top!``` can be directly called, but we need to use ```BlkBIST.bist_prbs_gen``` for the non-exported function.
"""

# ╔═╡ 0a622317-4a56-4b83-8aa5-8e808bb46136
md"""
Here a design decision is made to put all custom structs under a single module, and functions in another. The main reason is perhaps due to the No. 1 problem I have so far with Julia: it doesn't completely feel like a interpreted runtime language yet. Revise.jl package does a good job tracking real time updates in *functions*, but it can't do it for structs in a module (detailed [here](https://timholy.github.io/Revise.jl/stable/limitations/)). The solution is to put all structs in a dedicated module and include at the Main namespace [\[ref\]](https://www.youtube.com/watch?v=_3vJSBk0Bls&ab_channel=jmms). (It's ok if you don't follow fully - the key point is that it's some fundamental limitation to Julia's development environment right now).  

Another reason is due to Julia's multiple dispatch capability. It is encouraged to have functions outside of the scope of the custom data type definitions ([in constrast to conventional OOP](https://stackoverflow.com/questions/33755737/julia-oop-or-not?rq=4)).

Of course, once development is almost complete (or at least when you data types are not changing anymore), we can then fully modularize the custom structs for the final performance squeeze.

"""

# ╔═╡ 48d30a20-4cbc-4d14-be04-2f82938f000e
md"""
## File/Code Structure
"""

# ╔═╡ 669b166a-7d07-4741-b1b6-f266ef432b93
md"""
The ```src``` folder in the code base is then structured as the following
```bash
src
│
├─ blks
│ 	├─BlkBIST.jl
│ 	├─BlkCH.jl
│ 	├─BlkRX.jl
│ 	├─BlkTX.jl
│ 	└─WvfmGen.jl
│
├─ structs
│ 	└─TrxStruct.jl
│
├─ tb
│ 	└─TB.jl
│
├─ util
│  	└─Util_JLSD.jl
│
└─ Main.jl 
```

"""

# ╔═╡ e897900d-3df2-4ec0-8512-7303ba057091
md"""
The ```blks``` folders contain the modules of functions for each circuit block. ```structs``` contains the modules of data structs for each circuit block. ```tb``` has the run_sim functions that define how the circuit blocks are connected. The ```Main.jl``` file is main script to be run.
"""

# ╔═╡ e544f744-8111-438b-9094-08a14e6598f1
md"""
At this point, you should be able to understand the code structure and hopefully the run_sim() functions in TB.jl as well. Try running Main.jl (after Julia is setup) and see the following plot show up (with transmitter jitter turned on). In the next notebook, we will walk through the specifics of the ```BlkTX.jl``` to build a reasonably complex transmitter model.
"""

# ╔═╡ 2d9576ac-aefa-4ba3-bb21-24f0806abab4
md"""
![plot_w_jitter](https://circuit-artists.com/wp-content/uploads/2024/06/plot_w_jitter.png)
"""

# ╔═╡ 705f1580-cf83-4fe9-8091-f94b73a69c53
md"""
## Helper functions
"""

# ╔═╡ 436ebf41-8e28-466c-842b-cab12cc58ff8
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

# ╔═╡ 980169aa-1521-4721-bca6-e65f8a0cbf52
function pam_gen_top!(bist)
    @unpack pam, bits_per_sym, blk_size = bist.param
    @unpack polynomial, inv, gen_seed, gen_gray_map = bist
	@unpack So, So_bits = bist

	#generate PRBS bits
    So_bits, gen_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                	Nsym=bits_per_sym*blk_size, seed=gen_seed)

	#generate PAM symbols
    fill!(So, zero(Float64)) #reset So to all 0
    for n = 1:bits_per_sym
        @. So = So + 2^(bits_per_sym-n)*So_bits[n:bits_per_sym:end]
    end

    #gray encoding
    if ~isempty(gen_gray_map)
        for n in 1:blk_size
            So[n] = gen_gray_map[So[n] + 1]
        end
    end

    @. So = 2/(pam-1)*So - 1 #convert to analog voltage levels in +/-1 range

	return nothing
end

# ╔═╡ f1153611-34bb-4486-b606-3dd0245b0c8c
function ber_check_prbs!(bist)
    @unpack polynomial, inv, chk_seed, Si_bits = bist
    nbits_rcvd = lastindex(Si_bits)

	# Uncomment if you want to add artifical BER
    # err_loc = rand(Uniform(0,1.0), nbits_rcvd).< 1e-4;
    # Si_bits .= Si_bits .⊻ err_loc

    if bist.chk_lock_status
        ref_bits, chk_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                                 Nsym=nbits_rcvd,seed=chk_seed)

        bist.ber_err_cnt += sum(Si_bits .⊻ ref_bits)
        bist.ber_bit_cnt += nbits_rcvd
    else
        for n = 1:nbits_rcvd
            brcv = Si_bits[n]
            btst = inv
			for p in polynomial
            	btst ⊻= chk_seed[p]
			end

            #need consecutive non-error for lock. reset when error happens
            bist.chk_lock_cnt = (btst == brcv) ? bist.chk_lock_cnt+1 : 0

            chk_seed .= [brcv; chk_seed[1:end-1]]

            if bist.chk_lock_cnt == bist.chk_lock_cnt_threshold
                bist.chk_lock_status = true
                println("prbs locked")
                ref_bits, chk_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                                        Nsym=nbits_rcvd-n, seed=chk_seed)
                bist.ber_err_cnt += sum(Si_bits[n+1:end] .⊻ ref_bits)
                bist.ber_bit_cnt += nbits_rcvd-n
                break
            end
        end
    end

	return nothing
end

# ╔═╡ 25c78989-e38f-4525-b997-1641b2d47641
function int2bits(num, nbit)
    return [Bool((num>>k)%2) for k in nbit-1:-1:0]
end

# ╔═╡ c0f32334-53e5-43d4-9910-f134e994ec1c
function ber_checker_top!(bist)
    @unpack cur_blk, pam, bits_per_sym = bist.param
    @unpack gen_gray_map, chk_start_blk, Si, Si_bits = bist 


    if cur_blk >= chk_start_blk #make start blk a parameter later
        if ~isempty(gen_gray_map)
            for n in 1:blk_size
                Si[n] = gen_gray_map[Si[n] + 1]
            end
        end
    
        Si_bits .= vec(stack(int2bits.(Si, bits_per_sym)))

        ber_check_prbs!(bist)
    end

	return nothing
end

# ╔═╡ 80255c95-17dd-4ec5-b7a5-457b743deff8
@kwdef mutable struct Drv #dummy struct for TX driver
	param::Param
end

# ╔═╡ 1f797581-299e-4e52-b095-4893e499d401
drv = Drv(
		param = param,
		#TX driver params here
	);

# ╔═╡ 46d16141-5d63-4ede-a615-48aba2827a49
@kwdef mutable struct Ch #dummy struct for channel
	param::Param
end

# ╔═╡ e55cbdf5-c8ef-4372-935f-1203777fcb8a
ch = Ch(
		param = param,
		#Channel params here
	);

# ╔═╡ 5e584a09-a95f-4647-99ee-6f9b938ec705
@kwdef mutable struct Cdr #dummy struct for CDR
	param::Param
end

# ╔═╡ ee066af3-38d2-4ba2-8a71-dd1f6501c97d
cdr = Cdr(
		param = param,
		#CDR params here
	);

# ╔═╡ e22188b0-c747-4723-8e6e-874f2a564bab
trx = (;param, bist, drv, ch, cdr); #add more circuit modules if needed
#We can already "feel" that our transceiver is being built up ✌

# ╔═╡ b45061db-acb8-4f5e-be4e-bb8de116cc1b
trx.param.osr #example of accessing parameters

# ╔═╡ 79c7784c-e247-46d3-8810-8304f9766f17
trx.param === drv.param   #still referencing the same param

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
BenchmarkTools = "~1.5.0"
CairoMakie = "~0.12.0"
DSP = "~0.7.9"
DataStructures = "~0.18.20"
Makie = "~0.21.0"
Parameters = "~0.12.3"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "e62344340a33b1c9ae21b83a5b5b573fc75e70cd"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "588e0d680ad1d7201d4c6a804dcb1cd9cba79fbb"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "aec444a07f2b3df8d41a47fabd02841b32be2dc5"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.12.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "1755070db557ec2c37df2664c75600298b0c1cfc"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.0.3"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "2493cdfd0740015955a8e46de4ef28f49460d8bc"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.3"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "801aef8228f7f04972e596b09d4dba481807c913"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.4"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "fc713f007cff99ff9e50accba6373624ddd33588"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be50fe8df3acbffa0274a744f1a99d29c45a57f4"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "23ddd329f4a2a65c7a55b91553b60849bd038575"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.11"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "e96f6e1dba3c008d95b97103a330be6287411c67"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.21.0"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "f23e301d977e037ff8df4e1f5d8035cd78a1e250"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.8.0"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "1865d0b8a2d91477c8b16b49152a32764c7b1f5f"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "a14a99e430e42a105c898fcc7f212334bc7be887"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.4"

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

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "763a8ceb07833dd51bb9e3bbca372de32c0605ad"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.0"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d483cd324ce5cf5d61b77930f0bbd6cb61927d21"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.2+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "5d54d076465da49d6746c647022f3b3674e64156"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.8"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─79fd0370-2446-11ef-2922-21820062726d
# ╟─c34485a7-3583-4c48-945f-f94435f90257
# ╟─417676db-f650-4d85-84dc-465ee504a03d
# ╟─0d9aa8be-bf0d-41b2-91b5-d702195a3d18
# ╟─7d46974f-5531-4408-a96f-4672385665e1
# ╟─bb74fb24-4fb0-4fc6-b2c6-e9bf77f3d908
# ╠═df1878ff-130d-4db4-8803-2772237a8505
# ╠═f734a5c7-d0aa-4f65-b967-55bf859ad81e
# ╟─55db57f4-4918-4493-a71b-67cc4af609d6
# ╟─6466aa8f-c5bb-4643-b72c-ff687d992744
# ╠═5372f844-cedf-45e8-8799-106372630c10
# ╟─a9ab933a-20a8-4948-a662-d5338b697034
# ╠═282f989a-e9d5-45d2-a044-3d87df1c4827
# ╟─04002f49-604b-42b1-9b1c-21b5bd014240
# ╠═3ba76818-18a3-41cd-be27-fda0a7fede3d
# ╠═e65ff729-04ac-4303-bfe4-864f7066b7fe
# ╠═1136c126-c359-4c46-95d5-7d8b92c672ef
# ╠═ad83cf45-9b84-45bc-874f-5efbcc3bfa71
# ╟─fcbbbab5-5ed8-498b-b9af-005e1e112af7
# ╠═d812042d-e16d-4535-a358-ce6fdb710933
# ╠═d9522b57-cfb2-4108-a3da-846153141b12
# ╠═5d934cc6-0049-48e5-b4c2-92a7c7ae09fb
# ╠═9ec3d251-4ac3-4262-a1af-3f2b057bf155
# ╟─9335e8fb-9fde-414b-8c7c-b618b3b308e5
# ╠═e694af8f-8abf-408b-9df1-074a5c600256
# ╠═4755226e-d845-4f70-ae0e-fc2329d0c4c0
# ╟─51e571fd-367a-41bd-ab56-325b0bc626c5
# ╠═231e7fcc-ecc5-4dd4-b07d-1782413228ad
# ╠═934f66ca-6866-485a-aaab-68ab48e01574
# ╟─ecc52900-8d05-4104-a4c4-53f27b899a0f
# ╟─50f2f85a-9ec2-40c5-aa36-7ac8646148ac
# ╠═4e20103f-ef37-4172-9679-1a8ef1f8d4f8
# ╠═5d762756-2740-46c7-a45a-62e667d96ef6
# ╟─469aa4e7-50f4-488c-883e-6ef8acce7975
# ╠═68884750-d247-4dec-9019-7e0c179ebea4
# ╟─48af9bdd-960b-4e7c-8c62-4b04d6b8a82e
# ╠═6ce6f7f4-e6d1-4fc5-8391-707279f35696
# ╠═a3128c73-fcad-46d2-bf05-75feb1e3ab16
# ╟─54117b87-85c5-4757-a736-5dd0bfab8a14
# ╟─b52b20fc-7939-430f-bc1b-38e50e4ec1f5
# ╟─7c433885-9c9e-48ad-a331-1ee96f3e8119
# ╟─6d95615e-a9d0-40d6-9c54-e6a158d411f8
# ╠═fa01061f-c472-4ab8-93b1-645d84c8194d
# ╠═f3cf6e75-a0bd-4692-9430-37c39105148f
# ╠═ebf443e0-286d-437d-8bc6-9c52617b39ec
# ╟─b6b8344f-fbb0-4050-8917-824ca7c5545a
# ╠═1f797581-299e-4e52-b095-4893e499d401
# ╠═e55cbdf5-c8ef-4372-935f-1203777fcb8a
# ╠═ee066af3-38d2-4ba2-8a71-dd1f6501c97d
# ╠═e22188b0-c747-4723-8e6e-874f2a564bab
# ╟─5e83caa9-44a2-4f2c-9948-f64eb087cad5
# ╠═b45061db-acb8-4f5e-be4e-bb8de116cc1b
# ╠═79c7784c-e247-46d3-8810-8304f9766f17
# ╟─7f921c7d-f11e-46f9-99aa-14c55ad082b6
# ╟─be7e4564-630a-4617-bc37-dc270892c8bf
# ╠═c306ecd3-f521-420e-868a-2f7ca41b5fa6
# ╠═54c11b28-4550-4603-83e9-cb6188215f3e
# ╠═65c2f05e-b119-4869-89c2-7638c48f2494
# ╠═b59b371a-1e7e-4b9c-86b3-49493de40b34
# ╠═bdc0c012-541f-4f6b-9504-55e43bd75b8a
# ╠═e5746626-6860-480e-a45d-5640394f12cf
# ╟─2cc8af32-7897-46b4-88e1-2e4f40576823
# ╟─a7551ac2-c94f-48ff-8bea-9fe1172402e5
# ╠═980169aa-1521-4721-bca6-e65f8a0cbf52
# ╠═c0f32334-53e5-43d4-9910-f134e994ec1c
# ╟─8e43bdcb-7515-471c-9756-f5bfa8e6ed6d
# ╟─08035fd1-79d0-424b-b700-01ed1b10a156
# ╟─6dd4ab1f-4b80-48a6-a53a-e397b4eb580e
# ╟─297d11ed-1eea-4c65-ab62-7773dbb8f103
# ╟─c7b78049-1cef-46d7-aec5-de8775cb29a9
# ╟─1a47e550-62b4-4b5d-871b-07b2e3a5d1a0
# ╟─0a622317-4a56-4b83-8aa5-8e808bb46136
# ╟─48d30a20-4cbc-4d14-be04-2f82938f000e
# ╟─669b166a-7d07-4741-b1b6-f266ef432b93
# ╟─e897900d-3df2-4ec0-8512-7303ba057091
# ╟─e544f744-8111-438b-9094-08a14e6598f1
# ╟─2d9576ac-aefa-4ba3-bb21-24f0806abab4
# ╟─705f1580-cf83-4fe9-8091-f94b73a69c53
# ╠═436ebf41-8e28-466c-842b-cab12cc58ff8
# ╠═f1153611-34bb-4486-b606-3dd0245b0c8c
# ╠═25c78989-e38f-4525-b997-1641b2d47641
# ╠═80255c95-17dd-4ec5-b7a5-457b743deff8
# ╠═46d16141-5d63-4ede-a615-48aba2827a49
# ╠═5e584a09-a95f-4647-99ee-6f9b938ec705
# ╠═d2719225-8dd5-4385-b9d6-4e0f0d577a4e
# ╠═5cacc554-6c16-4f3e-a0bf-6acd8db7b1e7
# ╠═e6ef929a-1c1f-47c9-a1c2-888b84b641e4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
