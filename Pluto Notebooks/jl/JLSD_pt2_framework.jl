### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d6346d6d-7513-4d3f-a4e6-123cca197bd1
using Makie, CairoMakie

# ╔═╡ 2313e492-dfde-4e4e-9bfd-59e6a2a54b25
using DSP, Random, Distributions, UnPack, BenchmarkTools

# ╔═╡ 80f26a50-22aa-11ef-2af1-7799c685dae0
md"""
# Building SerDes Models in Julia, pt2 - Simulation Framework

![framework cover](https://circuit-artists.com/wp-content/uploads/2024/06/framework.png)

For the uninitiated, SerDes systems have evolved to become extremely complicated as the transmission speed increases (>100Gbps) and transmission medium changes (from copper to fiber). The challenge lies in the intersection between continuous- and discrete-time signal processing. One example is when clocks are used to sample signals (either implicitly with a data slicer or explicitly with data converters nowadays). Complex digital adaptation and clock-data recovery (CDR) algorithms become more interwined because all loops are running at the same time. How to effectively model and simulate all these components together has always been an interesting topic.

While many different analysis methods exist, including frequency and statistcal analysis, time domain results remain the final sign-off. Therefore, this framework focuses on transient simulation first and hopefully other analysis tools will be added in the near future.


"""

# ╔═╡ 59159363-df04-4e00-97b7-08dc74e6029f
md"""
## Key Simulation Parameters - time step and block size
"""

# ╔═╡ e42eb1ee-f93a-45ee-873e-57a282525b93
md"""
In order to simulate time domain waveforms, we resort to arrays (or vectors) to represent such waveforms. Now in order to build an *efficient* model, we first need to decide the granuarity of a simulation step, as well as the block size to simulate in batch. 
If we want to simulate a sequence of random bits, they need to be converted to continuous time waveform first with a finite time step. Then the step size choice directly impacts the accuracy and the speed of your sim. Allow me to illustrate
"""

# ╔═╡ 61cc9495-0c51-4058-b167-c758e40ab9c0
md"""
### Simulation time step
"""

# ╔═╡ adc81355-9249-4e2a-905f-7c9f65923000
tui = 1/10e9; # Unit Interval time for 10Gbps data

# ╔═╡ 7647e49a-5731-4772-8458-2585fe7f0980
osr1 = 4; # Oversampling ratio per symbol, i.e. # of points to simulate per symbol

# ╔═╡ 55d0d654-a061-44d6-9ae6-778281b00b16
osr2 = 32;

# ╔═╡ 8cb8ec3b-cb8c-4bb8-8378-489d304c8af5
dt1 = tui/osr1; # Simulation time step

# ╔═╡ 8f2cbd28-9d7e-4c4d-af93-881d59bb5864
dt2 = tui/osr2;

# ╔═╡ fe9b92ec-afe4-4d4a-af11-7d45f271f2cc
bits = [0,1,0,1,0,0,1,1,1]; #arbitrary bits; use bitrand(N) to generate random bits

# ╔═╡ 22905a47-35c5-4fea-92d8-5e8306141dbe
md"""
When ```osr``` is higher, the resultant waveform looks more like perfect square waveforms. Lower ```osr``` results in waveform distortion.
"""

# ╔═╡ 22021e9a-d3ed-444c-a722-4a03a7114a84
md"""
Here is what happens when we pass this waveform through a "channel" by means of convolution. With a larger time step, the convolution finishes faster because the array size is obviously smaller. However, there are noticeable kinks in the waveform, which might lead to simulation inaccuracy later. The accuracy increase also plateaus once the curve looks smooth enough.

Note that the increase in computation time doesn't increase linearly w.r.t array size becasue the convolution method in DSP.jl package uses FFT. Therefore, it's often more beneficial to just pick a bigger ```osr``` because Nlog(N) is cool! 
"""

# ╔═╡ 2aeb31ea-d262-4a58-b808-3dd66e413bb4
bw_ir = 8e9;

# ╔═╡ ca62496a-7f7c-428f-958b-1c317bd0f160
tlen_ir = 20*tui; #20UI long impulse response

# ╔═╡ c495f9db-3afa-4611-ab17-d9d5952a05b6
@time Vout1 = conv(ir1, Vbits1)*dt1;

# ╔═╡ c6f7940f-52f6-4bde-82d1-74d95ab066ef
@time Vout2 = conv(ir2, Vbits2)*dt2;

# ╔═╡ fa545192-e0c3-4188-a980-6cee173f99f5
md"""
### Simulation (sub-)block size
"""

# ╔═╡ 060687a7-a32b-4d14-aaf2-6b1b9b4bd1f3
md"""
Once the simulation time step is chosen, the next question is the batch size of the signal we want to process. The extreme examples would be when simulating 1 million bits, we can either generate and process *1 bit at a time* or *1 million at once*. This choice will equally impact simulation speed and accuracy. Like all things in engineering, it's always a spectrum and our framework should give the designer the freedom to adjust these parameters.
"""

# ╔═╡ 8da8b7c0-b36d-466d-87f4-5a7007cac5f1
md"""
Here we define two types of blocks: an outer loop block and an inner loop sub-block. Each block will contain different modules corresponding to different circuits in your system. Let's look at the example block diagram below.
"""

# ╔═╡ 0df54cc1-5c45-4b31-9274-5d7a2148400f
md"""
![sim blks](https://circuit-artists.com/wp-content/uploads/2024/06/sim_blks.png)
"""

# ╔═╡ b4c5c9ef-cd7d-46a5-b619-0ebf9f1e77b8
md"""
The outter loop block is suitable for circuits that are "quasi-static" and can process many bits at once. For example, a transmitter can generate waveforms for a large number of bits up front instead of simulating multiple short sequences, since we don't expect the tranmistter parameters to change during simulation (e.g. FIR coefficients and bandwidth). This will save time by doing repetitve convolutions (as we have seen before, longer convolutions become faster due to FFT algorithm).
"""

# ╔═╡ 0deea94e-fb8c-43bb-b0f8-f79f2fce6ae8
md"""
The inner loop sub-block is more suitable for circuits that require more accuracy and have more complicated loop dynamics. For example, the adaptation and CDR loops normally operate on a parallel bus of data and error slicer outputs. This vector then effectively determine the sub-block size because we want to capture the real dynamics of the feedback loops. Another way to think about this is that we implicitly created a "digital clock" through this sub-block loop, i.e. each loop iteration is a digital clock cycle.
"""

# ╔═╡ 657379de-8e3f-4eb5-9f8f-1a16ef4ed93e
md"""
![sim flow](https://circuit-artists.com/wp-content/uploads/2024/06/sim_flow.png)
"""

# ╔═╡ a86122f9-c706-4ec8-a1d8-62151ce6f983
md"""
The figure above illustrate the relationships between the big block and sub-blocks. We first batch process a long bit sequence, then iterate over the waveform through the sub-blocks.  

Note that although this picture implies that the block size has to be an integer multiple of sub-block size, but in reality it could be extended through an elastic buffer layer in between. The current code only assumes integer relationship between these two numbers, but the elastic buffer approach would be the most general for modeling say frequency offsets between TX and RX (will be addressed in future development)
"""

# ╔═╡ 27a5d616-267b-4689-b2f3-7e04d38514ac
md"""
## The run_sim() skeleton
"""

# ╔═╡ 910ad862-0baf-4ff0-8109-f19ae0688282
md"""
With this in mind, we can start with a run_sim() function pseudo-code already as
"""

# ╔═╡ 8775eb83-56ac-4874-8950-7cf872db4dd9
function run_sim(params)
	@unpack Nblks, Nsubblks, blk_size, subblk_size = params
	@unpack bist_params, tx_params, rx_params = params
	
	for n = 1:Nblks
		# run batch processing blocks here
		# e.g. 
		# S = bist_gen(bist_params)
		# Vtx = tx_drv(tx_params, Vtx)

		for m = 1:Nsubblks
			# run shorter loops here
			# e.g.
			# Φ = clkgen(pi_code)
			# Ssub = sampler(Vtx, Φ)
			# dvec = slicers(Ssub)
			# ...

			# record waveforms here
		end

		# bist_check(bist_params, data_rcvd)

	end
	
end

# ╔═╡ 1b8a3cbb-5b0b-47b2-945a-6e28a07c4ca6
md"""
**That's it!** This is pretty much the big picture of the framework, and now it's a matter of adding in the functional model of each circuit module **AS LONG AS** they satisfy this ```osr```, ```blk_size```, ```subblk_size``` framework. This usually means dealing with boundary conditions internally, remembering states from the previous (sub-)block, etc. 
"""

# ╔═╡ 80caa4e9-33b6-4a5e-b59b-ebc59fb1ec78
md"""
## World's best SerDes
"""

# ╔═╡ 4946bfa5-7041-453a-b3de-1242a4a8fa62
md"""
Let's build the world's most ideal SerDes using this framework. What we need are some BIST functions that generates and checks PRBS, an ideal transmitter and slicers. Let's start with the bist\_prbs\_gen from our [first notebook](https://github.com/kevjzheng/JLSD/blob/main/Pluto%20Notebooks/jl/JLSD_pt1_background.jl)
"""

# ╔═╡ bb44eaea-1452-491e-a4d8-8ffdcf62f533
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

# ╔═╡ 361ed8e3-5be7-4021-8ad1-7e7f67cf2255
md"""
And here is an example ber\_check\_prbs function
"""

# ╔═╡ c6610a7c-8a6a-4573-8eb7-adf8647a37ea
function ber_check_prbs(rcvd_bits; poly, inv, seed, lock_status, lock_cnt, lock_threshold, ber_err_cnt, ber_tot_cnt)
	
	nbits_rcvd = lastindex(rcvd_bits)
	
    if lock_status #if prbs already locked, use prbs_gen for reference bits
        ref_bits, seed = bist_prbs_gen(poly=poly, inv=inv,
										Nsym=nbits_rcvd, seed=seed)
        ber_err_cnt += sum(rcvd_bits .⊻ ref_bits)
		ber_tot_cnt += nbits_rcvd

    else # if not locked yet, use received bits as seed
        for n = 1:nbits_rcvd
            brcv = rcvd_bits[n]
            btst = inv
			for p in poly
            	btst ⊻= seed[p]
			end
            seed .= [brcv; seed[1:end-1]]

			#need consecutive non-error for lock. reset when error happens
            lock_cnt = (btst == brcv) ? lock_cnt+1 : 0

            if lock_cnt == lock_threshold
                lock_status = true
                println("prbs locked")
				#run prbs till end
                ref_bits, seed = bist_prbs_gen(poly=poly, inv=inv,
										Nsym=nbits_rcvd-n, seed=seed)
				ber_err_cnt += sum(rcvd_bits[n+1:end] .⊻ ref_bits)
				ber_tot_cnt += nbits_rcvd - n

                break
            end
        end
    end

	return seed, lock_status, lock_cnt, ber_err_cnt, ber_tot_cnt
end

# ╔═╡ a2fdc63e-9dbc-4d31-bb59-23575c0f3697
md"""
We can verify these two functions with direct feedthrough
"""

# ╔═╡ 3040f918-4a8c-4b3c-a497-75d730cc9365
poly = [28,31]; #PRBS31 polynomial

# ╔═╡ 744e09cb-a6aa-440b-82b8-00bb91bb1be0
inv = false;

# ╔═╡ 67304e59-8c2c-4b2f-a915-c387067c3f0a
Nsym = Int(10e6); #Try changing number of bits

# ╔═╡ 5030b487-04ed-4ecb-b643-6cc9ca076bf6
gen_seed = ones(Bool, 31);

# ╔═╡ 9760ad97-ef24-4eca-9cf6-2e180cd48840
chk_seed = zeros(Bool, 31);

# ╔═╡ 177e53b4-8291-478e-b7dc-a50bbaf518fe
bits_test, gen_seed1 = bist_prbs_gen(poly=poly, inv=inv, Nsym=Nsym, seed=gen_seed);

# ╔═╡ adcf57be-95a4-4f8b-802c-29f26c3ea78a
chk_seed1, lock_status1, lock_cnt1, ber_err_cnt1, ber_tot_cnt1 =
	ber_check_prbs(bits_test; poly=poly, inv=inv, seed=chk_seed, 
					lock_status=false, lock_cnt=0, lock_threshold=256,
					ber_err_cnt=0, ber_tot_cnt=0);

# ╔═╡ 5ba2a617-0b23-411d-bc57-5a4a6f2453b0
md"""
Try printing lock\_status, lock\_cnt, etc. below
"""

# ╔═╡ e450e9bb-bb24-4957-971e-720bcfb42c54
println(ber_err_cnt1)

# ╔═╡ c67a0e29-b615-4ffe-88a2-6a5bdc71f0a6
md"""
We can now artificially add bit errors to verify our checker
"""

# ╔═╡ 87b1c700-84e0-4649-b67f-82af4bea7054
ber_target = 1e-4; #play around with this number and Nsym

# ╔═╡ 9060e6e5-8265-44b5-8024-dd981f089253
err_loc = rand(Uniform(0,1.0),Nsym).< ber_target;
#Uniform is from Distributions package

# ╔═╡ 4826cfb3-70cc-409a-a004-533f9938df32
bits_werr = bits_test .⊻ err_loc; #flip the bits at error locations

# ╔═╡ 9b9244d0-bd60-4901-8ce4-6326af249a5a
chk_seed2, lock_status2, lock_cnt2, ber_err_cnt2, ber_tot_cnt2 =
	ber_check_prbs(bits_werr; poly=poly, inv=inv, seed=chk_seed, 
					lock_status=false, lock_cnt=0, lock_threshold=256,
					ber_err_cnt=0, ber_tot_cnt=0);

# ╔═╡ e6de1c80-083d-4ca2-a957-0fa91c5f449d
println("BER ≈ " * string(ber_err_cnt2/ber_tot_cnt2)) 

# ╔═╡ e6a8a14c-649c-4f76-b59b-ff84a5f9e095
md"""
Let's now build an ideal transmitter and receiver. The ideal transmitter will simply take the bit sequence and oversample according to ```osr``` to create the time domain waveform. Our ideal receiver will then take samples from the waveform with direct indexing. Because the waveform is perfect, it doesn't really matter which index we take the samples at.
"""

# ╔═╡ 6d4b04cc-97a9-4ede-a243-6c15d96522c7
function tx_drv_top(bits; osr, swing)
	#just like the gen_wvfm helper function, but no need for the time vector
	#change output waveform to analog voltage in +/- swing/2 
	#the dots are needed to do elementwise operation on the bits vector
	return kron(sign.(bits.-0.5), swing/2*ones(osr)) 

end

# ╔═╡ 44d7f2d8-ea37-42f0-b884-a6b16b013b92
function rx_slice_top(Vin; sample_idx, osr, slice_lvl)
	# sample_idx is a number between 1 and osr, will wrap around for now
	Vsamples =  Vin[mod(Int(sample_idx)-1, osr)+1:osr:end]
	return Vsamples .> slice_lvl
end

# ╔═╡ c5dd05d4-d11c-4592-9a0b-b7aafe00205f
md"""
Now let's write a simple run\_sim\_ideal() function to plug all these functions in. In this example, we don't have a sub-block loop, which would make more sense when CDR and adaptation is introduced.
"""

# ╔═╡ 471cadad-a7ae-44ea-aa16-278ca3fdd31d
@time Vtx = run_sim_ideal(32, 128, Int(1e6));
#play around with the osr and blk_size parameters and see how the run time is impacted

# ╔═╡ 8cc4b284-c473-4b83-9871-4ffefbba492e
begin
#plotting. We will learn more Makie plots later
fig_tx = Figure()
ax_tx = Axis(fig_tx[1,1])
lines!(ax_tx,Vtx, label = "TX output")
axislegend()
fig_tx
end

# ╔═╡ 920279e2-760f-491a-bacc-f563c4f59c68
md"""
Congratulations! We have just built your first SerDes model. I admit that this model isn't that useful, and the code itself is not looking that pretty. Trust me, it will get better 🙂.

The main issue here is that we have quite a few *states* that needs tracking (e.g. seed, ber_cnt, etc.), and spelling them out explicitly will get out of hand pretty quickly. In the next notebook, we will introduce **struct** in Julia, which is a lighter weight data structure than "classes" for those who are familiar with object oriented programming.

As a preview, we will use structs to store relevant parameters and states for the corresponding circuit module, and pass the structs around in functions. The TX driver function can then be written as the following
"""

# ╔═╡ 9fb818e6-42f6-414e-bb5e-98cbc22839eb
function tx_drv_top(tx, bits)
	@unpack osr, blk_size, swing = tx #grab the parameters we need from the tx struct

	return kron(sign.(bits.-0.5), swing/2*ones(osr))
end

# ╔═╡ 4e535c7f-4a51-46c8-8ca7-6fbc294a6fac
function run_sim_ideal(osr, blk_size, nsym_tot)

	# initialize parameters
	gen_seed = ones(Bool, 31)
	chk_seed = zeros(Bool, 31)
	poly = [28,31]
	inv = false

	lock_status = false
	lock_cnt = 0
	ber_err_cnt = 0
	ber_tot_cnt = 0

	ber_target = 0e-5;

	Vtx = zeros(blk_size*osr)
	
	#round nblk to the nearest total # of symbols
	nblk = Int(round(nsym_tot/blk_size)) 
	
	# main loop
	for n = 1:nblk
		#generate bits
		bits, gen_seed = bist_prbs_gen(poly=poly, inv=inv, 
										Nsym=blk_size, seed=gen_seed);
		
		#generate TX waveform
		Vtx = tx_drv_top(bits, osr=osr, swing=1)

		#sample TX waveform
		bits_rcvd = rx_slice_top(Vtx, sample_idx=osr/2, osr=osr, slice_lvl=0)

		#add random bit error
		err_loc = rand(Uniform(0,1.0),blk_size).< ber_target;
		bits_werr = bits_rcvd .⊻ err_loc;

		#check error
		chk_seed, lock_status, lock_cnt, ber_err_cnt, ber_tot_cnt =
			ber_check_prbs(bits_werr; poly=poly, inv=inv, seed=chk_seed, 
					lock_status=lock_status, lock_cnt=lock_cnt, lock_threshold=256,
					ber_err_cnt=ber_err_cnt, ber_tot_cnt=ber_tot_cnt);
		
		#record waveforms here if needed
	
	end

	println("Sim done")
	println("BER: " * string(ber_err_cnt/ber_tot_cnt)) 

	return Vtx
end

# ╔═╡ 30c55f3a-c070-405b-b70b-23f017d31f87
md"""
!!! info "Julia tips"
	Oops! I named the function tx\_drv\_top again, but with different arguments. No worries! Julia implements something called [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/) (learn more [here](https://www.youtube.com/watch?v=kc9HwsxE1OY&ab_channel=TheJuliaProgrammingLanguage)). A function can take arguments of different types, and share the same name. This would be particularly helpful when we want to optimize function performance and/or allow different usage cases (e.g. explicit arguments vs. struct).
	We won't go too much into the details of how to use multiple dispatch and typing to optimize Julia performance, but it will be relevant for improving the framework in the long run.
"""

# ╔═╡ a4703815-c47b-4ad7-a943-aafbc0e1cc7e
function tx_drv_top!(tx, bits)
	#grab the parameters we need from the tx struct
	@unpack osr, blk_size, swing = tx 

	#store output directly into the struct
	tx.Vout .= kron(sign.(bits.-0.5), swing/2*ones(osr)) 
end

# ╔═╡ 64450e30-819b-4dcf-8928-4fc2ec396bad
md"""
!!! info "Julia tips"
	Julia uses a convention of putting "!" in the function name to denote that it's a mutating function. This means that the argument passed into the function will be directly changed and no explicit values will be returned. 
	Depending on your programming philosophy, there are pros and cons to this approach. Directly mutating the arguments allows more efficient memory usage (i.e., no need to allocate new variables every time a function is called), but one needs to be careful about the sequence of function calls and keep track of the states.
	This simulation framework choses to mainly to use mutating functions, which we will cover more in the next notebooks.
"""

# ╔═╡ 5728e2cd-b8eb-4979-9f77-19f50f5c89e5
md"""

In the next notebook, we will focus more on the data structures and how the code base will be structured. We will start to build a more complete TX module as we discuss more Julia details. 

"""

# ╔═╡ 414814e4-6461-41ec-a313-d8a7ebf2aa0d
md"""
## Helper function section
"""

# ╔═╡ b2f16df4-9575-4234-b380-bbd9cd0e6230
function gen_wvfm(bits; tui, osr)
	#helper function used to generate oversampled waveform
	
	Vbits = kron(bits, ones(osr)) #Kronecker product to create oversampled waveform
	dt = tui/osr
	tt = 0:dt:(length(Vbits)-1)*dt

	return tt, Vbits
end

# ╔═╡ 533a9068-cf15-4075-a84e-0fa57a2f9ecd
tt1, Vbits1 = gen_wvfm(bits, tui=tui, osr=osr1);

# ╔═╡ 71870941-cbc9-4d79-a6b9-fd3ce55d21b9
tt2, Vbits2 = gen_wvfm(bits, tui=tui, osr=osr2);

# ╔═╡ 10bab6fb-e6ab-4581-b17b-489227372266
begin
#plotting. We will learn more Makie plots later
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, tt1, Vbits1, label = "osr=4")
lines!(ax, tt2, Vbits2, label = "osr=32")
axislegend()
fig
end

# ╔═╡ 609b8929-e9a2-477d-bbd6-8fdbfb63be6c
begin
#plotting. We will learn more Makie plots later
fig2 = Figure()
ax2 = Axis(fig2[1,1])
lines!(ax2, tt1, Vout1[1:length(Vbits1)], label = "osr=4")
lines!(ax2, tt2, Vout2[1:length(Vbits2)], label = "osr=32")
axislegend()
fig2
end

# ╔═╡ 7cae8f4e-2c32-4edc-9de8-d3411ff632ff
md"""
!!! info "Julia tips"
	The arguments in the gen_wvfm function are separated by a semicolon. This syntax in Julia means that anything after the semicolon is a keyword argument. An explicit "kw = something" statement is needed during function calls. Having keyword argument allows a function to have default parameters as well. 
"""

# ╔═╡ e8774cf5-a499-4d33-91a5-b57bc08793fb
function gen_ir_rc(dt,bw,t_len)
	#helper function that directly calculates a first order RC response, normalized to the time step
    tt = [0:dt:t_len-dt;]

	#checkout the intuitive symbols!
    ω = (2*π*bw)
    ir = ω*exp.(-tt*ω)
    ir .= ir/sum(ir*dt)

    return ir
end

# ╔═╡ 5834ebea-c215-4da8-b5de-14a4349177ee
ir1 = gen_ir_rc(dt1, bw_ir, tlen_ir);

# ╔═╡ e5f3b3e4-5db3-4343-97ad-0407f02fef49
ir2 = gen_ir_rc(dt2, bw_ir, tlen_ir);

# ╔═╡ 5f98a625-e315-4846-a121-ad0a5dca6fd8
md"""
!!! info "Julia tips"
	Unicode is supported in Julia, thus one can write cleaner code using commonly known symbols like in the gen\_ir\_rc function.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
BenchmarkTools = "~1.5.0"
CairoMakie = "~0.12.0"
DSP = "~0.7.9"
Distributions = "~0.25.108"
Makie = "~0.21.0"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "0afc2734d3211a1fc132a1f5fddb034ba85d3425"

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
deps = ["Random"]
git-tree-sha1 = "82b912bb5215792fd33df26f407d064d3602af98"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.2"

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
git-tree-sha1 = "22c595ca4146c07b16bcf9c8bea86f731f7109d2"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.108"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

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
# ╟─80f26a50-22aa-11ef-2af1-7799c685dae0
# ╟─59159363-df04-4e00-97b7-08dc74e6029f
# ╟─e42eb1ee-f93a-45ee-873e-57a282525b93
# ╟─61cc9495-0c51-4058-b167-c758e40ab9c0
# ╠═adc81355-9249-4e2a-905f-7c9f65923000
# ╠═7647e49a-5731-4772-8458-2585fe7f0980
# ╠═55d0d654-a061-44d6-9ae6-778281b00b16
# ╠═8cb8ec3b-cb8c-4bb8-8378-489d304c8af5
# ╠═8f2cbd28-9d7e-4c4d-af93-881d59bb5864
# ╠═fe9b92ec-afe4-4d4a-af11-7d45f271f2cc
# ╠═533a9068-cf15-4075-a84e-0fa57a2f9ecd
# ╠═71870941-cbc9-4d79-a6b9-fd3ce55d21b9
# ╟─22905a47-35c5-4fea-92d8-5e8306141dbe
# ╟─10bab6fb-e6ab-4581-b17b-489227372266
# ╟─22021e9a-d3ed-444c-a722-4a03a7114a84
# ╠═2aeb31ea-d262-4a58-b808-3dd66e413bb4
# ╠═ca62496a-7f7c-428f-958b-1c317bd0f160
# ╠═5834ebea-c215-4da8-b5de-14a4349177ee
# ╠═c495f9db-3afa-4611-ab17-d9d5952a05b6
# ╠═e5f3b3e4-5db3-4343-97ad-0407f02fef49
# ╠═c6f7940f-52f6-4bde-82d1-74d95ab066ef
# ╟─609b8929-e9a2-477d-bbd6-8fdbfb63be6c
# ╟─fa545192-e0c3-4188-a980-6cee173f99f5
# ╟─060687a7-a32b-4d14-aaf2-6b1b9b4bd1f3
# ╟─8da8b7c0-b36d-466d-87f4-5a7007cac5f1
# ╟─0df54cc1-5c45-4b31-9274-5d7a2148400f
# ╟─b4c5c9ef-cd7d-46a5-b619-0ebf9f1e77b8
# ╟─0deea94e-fb8c-43bb-b0f8-f79f2fce6ae8
# ╟─657379de-8e3f-4eb5-9f8f-1a16ef4ed93e
# ╟─a86122f9-c706-4ec8-a1d8-62151ce6f983
# ╟─27a5d616-267b-4689-b2f3-7e04d38514ac
# ╟─910ad862-0baf-4ff0-8109-f19ae0688282
# ╠═8775eb83-56ac-4874-8950-7cf872db4dd9
# ╟─1b8a3cbb-5b0b-47b2-945a-6e28a07c4ca6
# ╟─80caa4e9-33b6-4a5e-b59b-ebc59fb1ec78
# ╟─4946bfa5-7041-453a-b3de-1242a4a8fa62
# ╠═bb44eaea-1452-491e-a4d8-8ffdcf62f533
# ╟─361ed8e3-5be7-4021-8ad1-7e7f67cf2255
# ╠═c6610a7c-8a6a-4573-8eb7-adf8647a37ea
# ╟─a2fdc63e-9dbc-4d31-bb59-23575c0f3697
# ╠═3040f918-4a8c-4b3c-a497-75d730cc9365
# ╠═744e09cb-a6aa-440b-82b8-00bb91bb1be0
# ╠═67304e59-8c2c-4b2f-a915-c387067c3f0a
# ╠═5030b487-04ed-4ecb-b643-6cc9ca076bf6
# ╠═9760ad97-ef24-4eca-9cf6-2e180cd48840
# ╠═177e53b4-8291-478e-b7dc-a50bbaf518fe
# ╠═adcf57be-95a4-4f8b-802c-29f26c3ea78a
# ╟─5ba2a617-0b23-411d-bc57-5a4a6f2453b0
# ╠═e450e9bb-bb24-4957-971e-720bcfb42c54
# ╟─c67a0e29-b615-4ffe-88a2-6a5bdc71f0a6
# ╠═87b1c700-84e0-4649-b67f-82af4bea7054
# ╠═9060e6e5-8265-44b5-8024-dd981f089253
# ╠═4826cfb3-70cc-409a-a004-533f9938df32
# ╠═9b9244d0-bd60-4901-8ce4-6326af249a5a
# ╠═e6de1c80-083d-4ca2-a957-0fa91c5f449d
# ╟─e6a8a14c-649c-4f76-b59b-ff84a5f9e095
# ╠═6d4b04cc-97a9-4ede-a243-6c15d96522c7
# ╠═44d7f2d8-ea37-42f0-b884-a6b16b013b92
# ╟─c5dd05d4-d11c-4592-9a0b-b7aafe00205f
# ╠═4e535c7f-4a51-46c8-8ca7-6fbc294a6fac
# ╠═471cadad-a7ae-44ea-aa16-278ca3fdd31d
# ╟─8cc4b284-c473-4b83-9871-4ffefbba492e
# ╟─920279e2-760f-491a-bacc-f563c4f59c68
# ╠═9fb818e6-42f6-414e-bb5e-98cbc22839eb
# ╟─30c55f3a-c070-405b-b70b-23f017d31f87
# ╠═a4703815-c47b-4ad7-a943-aafbc0e1cc7e
# ╟─64450e30-819b-4dcf-8928-4fc2ec396bad
# ╟─5728e2cd-b8eb-4979-9f77-19f50f5c89e5
# ╟─414814e4-6461-41ec-a313-d8a7ebf2aa0d
# ╠═b2f16df4-9575-4234-b380-bbd9cd0e6230
# ╟─7cae8f4e-2c32-4edc-9de8-d3411ff632ff
# ╠═e8774cf5-a499-4d33-91a5-b57bc08793fb
# ╟─5f98a625-e315-4846-a121-ad0a5dca6fd8
# ╠═d6346d6d-7513-4d3f-a4e6-123cca197bd1
# ╠═2313e492-dfde-4e4e-9bfd-59e6a2a54b25
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
