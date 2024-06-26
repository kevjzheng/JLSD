### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 907c712d-8fc0-49ac-9c72-ae53f8c78794
using StatsBase, DSP, Interpolations, FFTW, MAT

# ╔═╡ 1fcec875-7dda-4afe-bd9e-a9fea816fc32
using Parameters, DataStructures

# ╔═╡ 1417102e-62a3-4999-a278-9314e24a5edb
using UnPack, Random, Distributions, BenchmarkTools

# ╔═╡ a2334355-f50c-456f-8a51-578885f07528
using GLMakie, Makie, PlutoUI

# ╔═╡ 0f491d87-d14a-4bfe-8a2d-821797926590
include("../../src/structs/TrxStruct.jl");

# ╔═╡ 031f8558-b375-469f-8c62-c7d92a75441f
include("../../src/util/Util_JLSD.jl");

# ╔═╡ ae7a47ae-2f31-426f-88ff-1cce13ff9622
include("../../src/blks/BlkBIST.jl");

# ╔═╡ 9c951989-97db-4ebd-9cd6-233beec23664
include("../../src/blks/WvfmGen.jl");

# ╔═╡ 248b7a30-25df-11ef-11c5-d33b88daa42b
md"""
# Building SerDes Models in Julia, pt4 - Detailed Transmitter Example
"""

# ╔═╡ a934351d-c8c8-475e-a7de-5dee9e45d549
md"""
![tx blk diagram](https://circuit-artists.com/wp-content/uploads/2024/06/tx_blk_diagram.png)
"""

# ╔═╡ 938cc465-1b5a-4e02-9956-a57cb9e04ccb
md"""
Let's build some meaningful models! We will use a transmitter as an example to show some important concepts in this simulation framework and Julia. 

To set the stage better, this notebook will have codes both in a cell (which you can run), and shown as just text. We have included/imported relevant source files into the notebook. The text block is just showing what the source code is to avoid variable definition collision, like below. The actual ```Drv``` struct is instantiated from the ```TrxStruct```  module.
"""

# ╔═╡ 42718f5c-44a5-4d0f-b6f0-b6db8a499d14
md"""
```Julia
@kwdef mutable struct Drv
    const param::Param

    ir::Vector{Float64}
    swing = 0.7

    fir::Vector{Float64} = [1,0]
    fir_norm = fir/sum(abs.(fir))

	#not used in the current model
    rlm_en = false
    rlm = 1.0
    quantize = false
    dac_res = 7 

    jitter_en = false
    dcd = 0.0
    rj_s = 0.0
    sj_amp_ui = 0.0
    sj_freq = 0.0
    last_sj_phi = 0.0

    Sfir_conv::Vector = zeros(param.blk_size+length(fir)-1)
    Sfir = @views Sfir_conv[1:param.blk_size]
    Sfir_mem = @views Sfir_conv[param.blk_size+1:end]
	Vfir = zeros(param.blk_size_osr)

    prev_nui = 4
    Δtt_ext = zeros(prev_nui+param.blk_size+1)
    Δtt = zeros(param.blk_size)
    Δtt_prev_nui = @views Δtt_ext[end-prev_nui:end]
    Vext::Vector = zeros(prev_nui*param.osr+param.blk_size_osr)
    V_prev_nui = @views Vext[end-prev_nui*param.osr+1:end]
    tt_Vext::Vector = zeros(prev_nui*param.osr+param.blk_size_osr)
    tt_uniform::Vector = (0:param.blk_size_osr-1) .+ prev_nui/2*param.osr

    Vo_conv::Vector = zeros(param.blk_size_osr+lastindex(ir)-1) 
    Vo = @views Vo_conv[1:param.blk_size_osr] 
    Vo_mem = @views Vo_conv[param.blk_size_osr+1:end]

end

```
"""

# ╔═╡ 01820462-23f1-498e-95ac-3843101fe4eb
md"""
## What's in the TX model?
"""

# ╔═╡ ef180811-bb81-4b2d-955f-3785308a9bc9
md"""
Let's go through what's included in the transmitter/driver before discussing how each portion is modeled.
"""

# ╔═╡ bdf79336-77db-491b-ba6d-049b33f088b6
md"""
```Julia
    ir::Vector{Float64}
    swing = 0.7
```
"""

# ╔═╡ fcbd062f-d907-4691-9365-7fa3fe68790d
md"""
First, the driver will have its own continuous-time impulse response to model its bandwidth. Here a time-domain vector ```ir``` is used since we can also use SPICE simulated impulse response if desired. The ```ir``` vector has to use the simulation time step in ```param``` (i.e., use the same ```osr``` per symbol). This can be done during initialization. ```swing``` specifies the peak-to-peak magnitude of the driver's output signal.
"""

# ╔═╡ b09b12d2-e7fb-4944-a680-a5f676176390
md"""
```Julia
    Vo_conv::Vector = zeros(param.blk_size_osr+lastindex(ir)-1) 
    Vo = @views Vo_conv[1:param.blk_size_osr] 
    Vo_mem = @views Vo_conv[param.blk_size_osr+1:end]
```
"""

# ╔═╡ f371095b-daa4-4b57-9c6d-274d32f801bd
md"""
The ```Vo_*``` vectors are used for storing the TX's output waveforms and will be passed to the next block. Here we introduce a ```"view"``` of an array. A ```view``` is essentially a pointer to a sub-section of another vector, but not a standalone vector itself. For example:
"""

# ╔═╡ 916ae0ce-04e6-4d78-9521-ff5501725bd2
begin
one2ten = collect(1:10); #this is to convert UnitRange to Array

#this creates a copy of one2ten[1:5] and assign it to one2five
one2five = one2ten[1:5]; 

#this says one2five_view points to the section 1:5 of one2ten, but it's not an independent vector
one2five_view = view(one2ten, 1:5);
end;

# ╔═╡ 1ae615a4-7a81-4abb-b3e7-2435d8ea2616
md"""
If now we change the first element of ```one2ten```, ```one2five``` will not change, but ```one2five_view``` will because it's pointing to a sub-array ```one2ten```
"""

# ╔═╡ efa51aac-75c3-4621-bf04-ff42a87f5bd4
begin
one2ten[1] = 10; #change the first element here to see how the other two changes

println(one2five[1])
println(one2five_view[1])
end

# ╔═╡ 5b16cc5d-6de8-4126-ab11-bfb0cabf44b7
md"""
!!! info "Julia Tips"
	```@views``` is a macro that converts sliced arrays into views (pointers are much cheaper than creating copies of arrays). For more information on how to use the ```view``` syntax correctly, check [here](https://www.juliabloggers.com/the-view-and-views-macros-are-you-sure-you-know-how-they-work/)
"""

# ╔═╡ f0209a9a-ddfd-419f-85ae-fecf4882d7ab
md"""
With views, we can declare only one memory space for ```Vo_conv```, and create convience variables to make the first ```blk_size_osr``` length vector the actual output vector, and the rest becomes the memory vector that overalps with the beginning of the next block. The length for ```Vo_conv``` is known beforehand given the length of the ```ir``` vector and ```blk_size_osr```. 

This way, there is no need to explicitly copy and define ```Vo = Vo_conv[1:blk_size_osr]``` during our simulation to save space and time (something that's not possible in MATLAB). We will revisit this again when writing the convolution function later.
"""

# ╔═╡ 0f373aa6-55c9-4e48-abf0-019373f36a9b
md"""
![vector views for convolution](https://circuit-artists.com/wp-content/uploads/2024/06/vo_views.png)
"""

# ╔═╡ 0e155843-85fe-4573-afce-0fc11f56d3e8
md"""
```Julia
    fir::Vector{Float64} = [1,0]
    fir_norm = fir/sum(abs.(fir))
```
"""

# ╔═╡ e36746f1-19f4-4ecf-a6a0-29be5c307e26
md"""
```Julia
	Sfir_conv::Vector = zeros(param.blk_size+length(fir)-1)
    Sfir = @views Sfir_conv[1:param.blk_size]
    Sfir_mem = @views Sfir_conv[param.blk_size+1:end]
	Vfir = zeros(param.blk_size_osr)
```
"""

# ╔═╡ 4e8daf25-af1c-4ac8-8436-7d971fd9c004
md"""
```fir``` is a small vector containing the discrete-time coefficients for TX de-emphasis. ```fir_norm``` is the normalized coefficients to model the peak power constraint. Here we will just calculate it internally during initialization. Going along with ```fir``` are the internal vectors used for convolution, ```Sfir_conv```. Similarly, ```Sfir``` and ```Sfir_mem``` are the sub-array views.
```Vfir``` is just the oversampled version (by ```osr``` times) of the filtered sequence ```Sfir```.
"""

# ╔═╡ cbcd9551-da4c-4c91-bec6-a698f5ec6fa5
md"""
```Julia
    jitter_en = false
    dcd = 0.0      		#0.01 = 1%
    rj_s = 0.0 			#random jitter in seconds
    sj_amp_ui = 0.0 	#sinusoidal jitter amplitude in UI
    sj_freq = 0.0      	#sinusoidal jitter frequency in Hertz
    last_sj_phi = 0.0 	#state variable for previous blk's last SJ phase
```
"""

# ╔═╡ f0495ec1-04f5-4d72-82f9-35c768927e8e
md"""
```Julia
    prev_nui = 4
    Δtt_ext = zeros(prev_nui+param.blk_size+1)
    Δtt = zeros(param.blk_size)
    Δtt_prev_nui = @views Δtt_ext[end-prev_nui:end]
    Vext::Vector = zeros(prev_nui*param.osr+param.blk_size_osr)
    V_prev_nui = @views Vext[end-prev_nui*param.osr+1:end]
    tt_Vext::Vector = zeros(prev_nui*param.osr+param.blk_size_osr)
    tt_uniform::Vector = (0:param.blk_size_osr-1) .+ prev_nui/2*param.osr
```
"""

# ╔═╡ bb49149b-db9c-43d2-bf17-6862c26ebc11
md"""
These parameters/variables are used for modeling TX jitter. Currently the model supports duty cycle distortion (```dcd```), normally distributed random jitter (```rj```) and sinuosoidal jitter (```sj```). The parameters are defined in their "most comfortable" units.

```prev_nui``` denotes the number of previous symbols to be stitched to the current block's signal to prevent overflow/underflow when jitter is introduced. ``` Δtt*``` vectors store the jitter information at each edge location. The ```tt_Vext``` vector is the jittered time grid vector. ```tt_uniform``` is the convience vector to remap the jittered waveform back to our simulation grid. ```Vext``` and ```V_prev_nui``` are the extended block vectors and the previous N-UI symbols. More details on how these are used to model jitter later. 
"""

# ╔═╡ 34319afe-4569-4cf1-9df4-21f976388d1c
md"""
```Julia
	#not used in the current model
    rlm_en = false
    rlm = 1.0
    quantize = false
    dac_res = 7 
```
"""

# ╔═╡ 19c8f09b-e1f0-4117-8c6e-3623988b05fa
md"""
What's also possible is to model nonlinearity and quantization errors in the TX driver (if a DAC based model is desired). In fact, it would be a good exercise to do to extend this model and learn Julia on your own 😃. 
"""

# ╔═╡ c109ab44-1cc8-41f2-8f9b-2c04dbf824b7
md"""
Ok, time to extend the relevant structs
"""

# ╔═╡ 6a733af2-63f2-418c-8100-57195d9d790a
md"""
## The ```drv_top!``` function
"""

# ╔═╡ eb15958d-acd1-484f-972c-b51f49352dfd
md"""
Let's begin with some pseudo-code according to the block diagram at the very top. The ```drv_top!``` function would take the ```drv``` struct as an input (which will contain all the necessary parameters, internal states and pre-allocated output vectors), and a bit sequence vector from the BIST block.
"""

# ╔═╡ 111848d3-e97a-4869-807c-c9f927d2ca58
md"""
```Julia
function drv_top!(drv, input)
	@unpack all parameters and vectors

    apply_fir_filter!(Sfir, input, fir, kwargs...)

	oversample!(Vfir,Sfir)

	if jitter_en
		add_jitter!(drv, Vfir)
	end

	convolve!(Vo_conv, Vfir, ir, kwargs...)

end
```
"""

# ╔═╡ 4221ddc6-5958-4e75-9052-55929f7f6d74
md"""
We will begin with the convolution function since it can belong to the utility module and called by many other circuit blocks (anything with an impulse response really). 
"""

# ╔═╡ e103e442-69d3-44c6-8670-053534322412
function u_conv!(Vo_conv, input, ir; Vi_mem = Float64[], gain = 1)
    Vo_conv[eachindex(Vi_mem)] .= Vi_mem
    Vo_conv[lastindex(Vi_mem)+1:end] .= zero(Float64)

    Vo_conv .+= conv(gain .* input, ir)
    return nothing
end
#we will use a "u_" prefix for all utility functions to avoide name collision in the future

# ╔═╡ 022cf285-2a25-407a-a322-5cbab294aeba
md"""
A quick recap: the ```!``` is a Julian *convention* that denotes the function will mutate one or more of the arguments (usually the first argument). Here we pass in the pre-allocated output vector,```Vo_conv```, ```input``` and ```ir```. Optional arguments include a memory vector (from the previous block) ```Vi_mem```, and a ```gain``` factor. The gain factor comes in handy for general normalization (i.e. multiply by ```dt``` in convolution).
"""

# ╔═╡ 6ac67d2d-79ae-4f9e-9328-c5a5981e7680
md"""
Note that I didn't need to assign the specific sub-arrays of ```Vo_conv``` to ```Vo``` and ```Vo_mem```. This is automatically maintained with ```views```. We also directly use the ```conv``` function from DSP.jl since it's quite optimized with FFT. For continuous-time convolutions, the input and ir vectors could be pretty long, so FFT-based convolution is more suitable.
"""

# ╔═╡ 3adc12e0-1ec7-4ea2-bc2a-dc35a987d96f
#we will also create a non-mutating u_conv function for other uses
function u_conv(input, ir; Vi_mem = zeros(1), gain = 1)
    vconv = gain .* conv(ir, input)
    vconv[eachindex(Vi_mem)] += Vi_mem

    return vconv
end

# ╔═╡ f389a6fd-c388-4e9c-8d66-37c00305c1a6
md"""
Let's test it out:
"""

# ╔═╡ 2ab349bc-07af-4056-a729-37db20705389
md"""
For those who want to see eye diagrams, a helper eye diagram generation (heatmap based) function is included and we can plot the TX output like below (make sure to increase ```blk_size``` above to get more samples. Re-run cells if no change is seen). If the TX waveform is pretty "clean", the generated eye diagram might be a bit hard to see.
"""

# ╔═╡ 5062c917-f265-4e63-98a3-544bb9482bc1
begin #common eye diagram params
x_npts_ui = 256; #number of points per ui
x_nui = 2;
x_npts = x_nui*x_npts_ui; #plot two UI eye diagram
x_grid = -x_nui/2 : 1/x_npts_ui: x_nui/2-1/x_npts_ui;
y_range = 2; 
y_npts = 256; #number of points on y axis
y_grid = -y_range/2: y_range/y_npts: y_range/2 - y_range/y_npts;
end;

# ╔═╡ 85f612fa-1889-4631-a334-27530881d7b8
md"""
That wasn't too bad, was it? Actually, we can directly use the ```u_conv!``` function for applying the FIR filter as well. However, the FIR filter typically is much shorter (<10 taps) than the symbol vector, using FFT convolution might be an overkill. For optimization, a simple shift-and-add filter function can be written as below
"""

# ╔═╡ 91469577-cb17-4505-bdcd-bafce396d564
function u_filt!(So_conv, input, fir; Si_mem=Float64[])
    So_conv[eachindex(Si_mem)] .= Si_mem
    So_conv[lastindex(Si_mem)+1:end] .= zero(Float64)
    s_in = lastindex(input)
    
    for n=eachindex(fir)
        So_conv[n:s_in+n-1] .+= fir[n] .* input
    end

    return nothing
end

# ╔═╡ 4651d7ef-3d52-4275-a0ce-13051d1a1cdf
#non-mutating version
function u_filt(input, fir; Si_mem=Float64[])
	sconv = zeros(length(input) + length(fir) - 1)
	
    s_in = lastindex(input)
    
    for n=eachindex(fir)
        sconv[n:s_in+n-1] .+= fir[n] .* input
    end

	sconv[eachindex(Si_mem)] .+= Si_mem

    return sconv
end

# ╔═╡ 03ec64bb-7bb7-4e7f-9b4e-00ec713d739c
md"""
Now let's apply some FIR to open up the eye. Here we will use the non-mutating functions to maintain the states of the drv struct above.
"""

# ╔═╡ f96baca0-637c-4b02-855d-d758121fa96e
md"""
Driver bandwidth: $(@bind drv_bw PlutoUI.Slider(1e9:0.5e9:10e9, 10e9, true))

Post tap value: $(@bind posttap PlutoUI.Slider(-0.5:0.02:0, 0.0, true) )
"""

# ╔═╡ 01549352-1220-4e5b-ae1f-143a5bcada0c
fir_coeffs = [1.0, posttap];

# ╔═╡ bfee37c6-fa4c-4e08-9b76-cfa298fff745
fir_coeffs_norm = fir_coeffs ./ sum(abs.(fir_coeffs));

# ╔═╡ c67fbbff-19e3-4ffd-bafd-d934e4882268
md"""
Cool! Actually, that was the easy part. Because we have defined a good ```drv``` struct and had clever uses of ```views```, the functions seem quite straightforward and gave us interesting results to play with immediately. 
"""

# ╔═╡ dc08b7fa-f12a-4aa1-ac9d-049583cb91c1
md"""
## Let's be jittery
"""

# ╔═╡ 98fb5603-fe2e-455d-8d28-f96a8113121d
md"""
Too much coffee when I wrote this notebook? You bet! That's why our TX needs to be jittery now. For a refresher on the main types of TX jitter, [click here](https://people.engr.tamu.edu/spalermo/ecen689/lecture10_ee720_jitter.pdf) for Prof. Palermo's lecture slides on jitter.
"""

# ╔═╡ 55a76e1a-b520-4863-8864-ce2e772742ee
md"""
Jokes aside, modeling time domain phenomenon is always a challenge. How do we model jitter when our simulation time step is "fixed"? The key insight is that in our simulation framework, the zero crossing in the oversampled waveform is implicitly between the finite steps, i.e. @ t = N*osr+0.5. Despite our simulation time step being discrete, we can still use interpolation to "encode" where our zero crossing should be. Note that by using intermediate voltages, we can embed jitter information at fractional time steps.

![embed_jitter](https://circuit-artists.com/wp-content/uploads/2024/06/emdedded_jitter2.png)
"""

# ╔═╡ 24f1ef71-2d2f-42ff-a018-e2922a185de4
md"""
The trick then is to **warp or remap** a "jittery time grid" onto our "uniform time grid", for lack of better terms. We will first generate the Δt at each nominal edge transition.
"""

# ╔═╡ 93413781-bc43-45d2-8d35-0017b948b189
md"""
Let's start with duty cycle distortion. When there is a distortion of ```dcd```, it means all even edges shift by ```dcd/2*osr```, and all odd edges shift by ```-dcd/2*osr``` (the sign here doens't really matter). It's also best to set the ```blk_size``` as an even number because we are modeling duty cycle.
"""

# ╔═╡ 74ffbe27-ae1b-4dac-84b8-83f074038a2e
md"""
Now we add random jitter
"""

# ╔═╡ 22373810-8aa9-4950-9c41-e1cf08372421
md"""
Time for sinusoidal jitter. We first define the phase of jitter (```phi_sj```), then pass it into a sine function.
"""

# ╔═╡ 4c28d8ef-9169-4c9a-8a03-bc0f67e52991
md"""
!!! info "Julia Tips"
	Julia supports syntax like ```2π```! It will understand it as ```2*π```. Makes mathematical programming much cleaner. Latex like syntax (e.g. Φₛⱼ, created by \Phi + Tab + \ \_s + Tab + \ \_j + Tab) is also supported.
"""

# ╔═╡ 6be5384f-b06e-494e-a659-f772438e25d6
md"""
Now we have generated all the jitter information for each transition edge. Next step is to create the finer grid to go with our voltage waveform. First, we need to extend our Δtt vector with some more samples from previous UIs to avoid overflow/underflow.
"""

# ╔═╡ 510dec4a-07eb-4fa7-addb-876e3c171587
md"""
Note this is a similar style as our convolution with memory. The Δtt\_prev\_nui vector is automatically taken care of for the next block due to views. We will also define our ```Vext``` here.
"""

# ╔═╡ c66d9219-4f70-4917-855f-b3ff1c1c73f4
md"""
We now build a helper function to create a linear grid in between the transition times.
"""

# ╔═╡ 63b400e2-c7d2-4873-ad93-6a23d14b760b
function drv_jitter_tvec!(tt_Vext, Δtt_ext, osr)
    for n = 1:lastindex(Δtt_ext)-1
        tt_Vext[(n-1)*osr+1:n*osr] .= LinRange((n-1)*osr+Δtt_ext[n], n*osr+Δtt_ext[n+1], osr+1)[1:end-1]
    end

    return nothing
end

# ╔═╡ b94044fa-193a-4bb9-998b-02e74bbee98c
md"""
Pictorially, this function is doing the following
![jitter_time_axis](https://circuit-artists.com/wp-content/uploads/2024/06/jitter_time_axis.png)
"""

# ╔═╡ f063e527-2363-4692-b79a-9646a9b68723
md"""
Great! Now we created a jittered time axis for our voltage waveform. If we plot ```tt_Vext``` vs. ```Vext```, this will represent the waveform with jitter! That's because we are viewing this plot from a uniform time grid perspective. So the last step then is to remap this signal back onto our uniform simulation time grid, through (drum roll....) **interpolation**. 
"""

# ╔═╡ 14ce9a21-09de-4ec6-9511-5a4281459032
md"""
Julia has a interpolation package, and slightly different way of doing interpolation compared to MATLAB. We will stick with the simple linear interpolation since accuracy can always be adjusted with increasing ```osr```. 
"""

# ╔═╡ cebb0066-edf7-41c2-99ea-bd5c966fe04b
md"""
Julia's interpolation return a *function object* that can operate on any values you throw at it. Compared to MATLAB, the ```interp``` function takes the old axis/value and new axis together. Julia's interp function object comes handy when repeated interpolation is needed.
"""

# ╔═╡ a9785c06-83d2-46a4-aec1-63d0644c19df
md"""
!!! info "Julia Tips"
	Julia can be C/C++ like if you want to optimize for performance. Using the built-in ```linear_interpolation``` can make your code functional at first, but might have too much overhead (i.e. checking for conditions and inputs to make the right internal call). It's possible for you to write a specialized interpolation function for this case. As an example, check out the ```drv_interp_jitter!``` function in the appendix.
"""

# ╔═╡ d8da233c-81ef-4916-88e3-2b9b0fac5a57
md"""
To summarize and not interfere with previous codes in the notebook, you can play around with the code cell below to see how the jittered eye diagram change.
"""

# ╔═╡ 690ea232-0ef4-45ab-a158-6ff355d8101f
md"""
Parameters for you to play around with

DCD: $(@bind dcd PlutoUI.Slider(0:0.01:0.1, 0.0, true))

RJ (seconds): $(@bind rj_s PlutoUI.Slider(0:500e-15:5e-12, 0.0, true))

SJ amplitude (UI): $(@bind sj_amp_ui PlutoUI.Slider(0:0.02:0.3, 0.0, true))

SJ frequency (Hertz): $(@bind sj_freq PlutoUI.Slider(0:200e3:10e6, 2e6, true))

"""

# ╔═╡ 0c05c7de-ff91-4927-8f8e-44340f014d55
md"""
It's important to note that the impulse response of the driver (and subsequent channel, RX front end, etc.) plays a crucial role in low-pass filtering the jittered waveform to give it a "smoother look".
"""

# ╔═╡ 6d5ddcaf-bab8-44d8-a2f9-9a90bbe25fdf
md"""
## Analyzing jitter
"""

# ╔═╡ a296059d-7228-48da-b538-cac3e0c06ab7
md"""
The utility module includes some helper functions to analyze jitter. We can use the find_0x (find zero crossing) function to generate jitter statistics of any waveform.
"""

# ╔═╡ d130e739-b6d2-445e-ad27-510fb34da8bc
md"""
Because of normalization, there could be big jumps in the zero crossing points (try increasing the SJ amplitude and see). The unwrap_0x function in utility module can help unwrap these jumps.
"""

# ╔═╡ 24be3bdb-c3db-42ab-b2ef-c313b7108ff8
md"""
Lastly, we can now plot the distribution of the jitter! Ideally, there could be another family of functions that take an distribution and extract dcd, RJ, SJ, etc. (just like a test instrument). All these functionalities could be added in the future (by anyone!).
"""

# ╔═╡ 78de21b5-984b-4e1f-9759-dfef88885c48
md"""
## Putting it all together
"""

# ╔═╡ 1dfc2b52-1749-4782-9c5b-7c3289069803
md"""
We now reach the point to put everything together. Below is the ```dac_drv_top!``` function showing everything we have covered. Refer to the source code in repository for more details on the ```drv_add_jitter!``` (which is just an encapsulation of the jitter section)
"""

# ╔═╡ bb39273b-820f-472a-a6ae-c6ed9094044f
function dac_drv_top!(drv, Si)
    @unpack osr, dt, blk_size, blk_size_osr = drv.param
    @unpack ir, fir_norm, swing, Sfir_mem, Vo_mem = drv

	#FIR filter
    u_filt!(drv.Sfir_conv, Si, fir_norm, Si_mem = Sfir_mem)
	
	#Oversample
    kron!(drv.Vfir, drv.Sfir, ones(osr))

	#Add jitter
    if drv.jitter_en
        drv_add_jitter!(drv, drv.Vfir)
    end

	#Convolve w/ ir
    u_conv!(drv.Vo_conv, drv.Vfir, ir, Vi_mem=Vo_mem, gain=dt*swing/2)

end

# ╔═╡ bb408a87-8ded-403e-87a8-2a9cf42e79fb
md"""
## Conclusions
"""

# ╔═╡ 555d708b-ab85-4345-9e90-d6d5f480cbb1
md"""
In this notebook, we covered the important concepts in a transmitter of model, including de-emphasis FIR, impulse responses, and jitter. A lot more is yet to be done, like nonlinearity and quantization.

The key modeling methods also directly apply to other circuit modules as well. Take advantage of views in Julia for more efficient use of memory. Convolution is a common function that will be shared among pretty much all blocks. Interpolation is important for transforming time axis as well as sampling (on the RX side). 

Though not directly related to the simulation framework itself, plots are important tools when it comes to visualization and debug. We will cover more in depth about the plotting packages in Julia in the next one so that you can start generating great plots too! 
"""

# ╔═╡ 4deb9728-b83e-4f7e-8fd9-2cdf215c6408
md"""
## Helper Functions
"""

# ╔═╡ 3bf00ba7-1593-4aad-bd92-49095c1f65e5
function u_find_0x(input) #vectorized implementation
    sign_input = sign.(input)
    diff_sign = @views sign_input[2:end] .- sign_input[1:end-1]
    x_idx_crs = findall(abs.(diff_sign) .> 1 )
    x_idx_fine = Vector{Float64}(undef, lastindex(x_idx_crs))

    @. x_idx_fine = x_idx_crs+input[x_idx_crs]/(input[x_idx_crs]-input[x_idx_crs+1])
    
    return x_idx_fine
end

# ╔═╡ 8f763793-3e0c-4d72-988c-e5fe61585135
md"""
```Julia
function u_find_0x(input) #explicit for-loop implementation
    sign_input = sign.(input)
    diff_sign = @views sign_input[2:end] .- sign_input[1:end-1]
    x_idx_crs = findall(abs.(diff_sign) .> 1 )
    x_idx_fine = Vector{Float64}(undef, lastindex(x_idx_crs))

    for n = eachindex(x_idx_crs)
        x_crs = x_idx_crs[n]
        x_idx_fine[n] = x_crs+input[x_crs]/(input[x_crs]-input[x_crs+1])
    end
    return x_idx_fine
end
```
"""

# ╔═╡ 41cb22b3-f000-4326-91ce-e52e95048cf0
function u_unwrap_0x(xpts; tol_Δui = 0.5) #assumes 0-1UI range, vectorized
    nwrap = 0
    xpts_unwrap = zeros(lastindex(xpts))
    xpts_unwrap[1] = xpts[1]

    ΔΦ = @views xpts[1:end-1] .- xpts[2:end]

    nwrap = cumsum( (abs.(ΔΦ) .> tol_Δui) .* sign.(ΔΦ))

    xpts_unwrap[2:end] .= nwrap .+ @views xpts[2:end]

    return xpts_unwrap
end

# ╔═╡ 724709c6-1deb-4dd9-8d83-6db16a9bd2b6
md"""
```Julia
function u_unwrap_0x(xpts; tol_Δui = 0.5) #explicit for-loop implementation
    nwrap = 0
    xpts_unwrap = similar(xpts)
    xpt_prev = xpts[1]
    
    for n = eachindex(xpts)
        xpt = xpts[n]
        if abs(xpt-xpt_prev) > tol_Δui
            nwrap -= sign(xpt-xpt_prev)
        end

        xpts_unwrap[n] = nwrap + xpt

        xpt_prev = xpt
    
    end
    return xpts_unwrap
end
```
"""

# ╔═╡ 610b7b53-bd09-4a02-8e60-f43209da4128
import .Util_JLSD: u_gen_ir_rc

# ╔═╡ e97cfdc6-c2d8-4f68-bb0d-39117b02d560
#run this cell if you don't see any waveforms below
begin
param = TrxStruct.Param(  
			data_rate = 10e9,
			pam = 2,
			osr = 32,
			blk_size = 2^14,
			subblk_size = 32, 
			nsym_total = Int(1e6));

bist = TrxStruct.Bist(  
			param = param,
			polynomial = [28,31]);

drv = TrxStruct.Drv(
		param = param,
		ir = u_gen_ir_rc(param.dt, param.fbaud/4, 20*param.tui), #1st order ir
		fir = [1.0, -0.2],
		swing = 0.8,
		jitter_en = true,
		dcd = 0.03,
		rj_s = 300e-15,
		sj_amp_ui = 0.0,
		sj_freq = 10e6);
end;

# ╔═╡ 7644dc88-7529-4a0e-8646-b024dad540b3
Vosr = kron(bist.So, ones(param.osr)); 
#oh yeah, oversampling is this simple too

# ╔═╡ bb94d0a1-48f5-401e-b12a-2378435d4422
begin
f = Figure()
lines!(Axis(f[1,1]),Vosr, label="Vosr")
axislegend()
f
end

# ╔═╡ 7b13b871-3b0b-42c3-b9a7-00cf246c86b1
#call our convolution function; let's keep the input memory zero for now
#change the drv parameters in the struct definition to see the waveform/eye change
u_conv!(drv.Vo_conv, Vosr, drv.ir, Vi_mem=zeros(1), gain = drv.swing * param.dt);

# ╔═╡ e653c634-02db-44b3-8501-e2614ee2c07c
begin
f1 = Figure()
lines!(Axis(f1[1,1]),drv.Vo, label="drv output")
axislegend()
f1
end

# ╔═╡ 2a9645ac-a66e-4133-8dc7-2b455e0473f5
s_filt = u_filt(bist.So, fir_coeffs_norm); 
#change the FIR filter here to see the eye change reactively!

# ╔═╡ 0984e74f-c480-4aac-b429-fae9d10afd6d
Vosr_filt = kron(s_filt, ones(param.osr));

# ╔═╡ d3f0e565-51c3-497d-bd59-b0acb842e402
Δtt = zeros(param.blk_size); #clean slate no jitter for all symbols

# ╔═╡ e7e4e516-f575-4403-ad0a-cf0511d549b0
Δtt[1:2:end] .+= drv.dcd/2*param.osr;

# ╔═╡ f549409a-4ac7-4244-8d12-44d1eabe566b
Δtt[2:2:end] .+= -drv.dcd/2*param.osr;

# ╔═╡ b0c39645-22ec-49fb-ae84-900b169fdde7
rj_osr = drv.rj_s/param.tui * param.osr; #conver RJ unit from seconds to osr

# ╔═╡ b1ef7a3d-1eb3-4703-9827-8ede9165aab9
Δtt .+= rj_osr .* randn(param.blk_size);

# ╔═╡ a2c309cc-9fea-4ab0-909e-85f17ffa4e10
sj_amp_osr = drv.sj_amp_ui * param.osr; 

# ╔═╡ 8a9eefca-158b-400f-a7fa-2fcf18c02324
sj_freq_norm = drv.sj_freq * param.tui; #normalize SJ frequency to symbol rate

# ╔═╡ a5ec92b5-a0ee-45b4-b017-4822863b2055
begin
phi_sj = (drv.last_sj_phi .+ (2π*sj_freq_norm) * (1:param.blk_size)) .% (2π);
drv.last_sj_phi = phi_sj[end]; #store the last phase for next block
end;

# ╔═╡ 585edc1c-3349-4888-a8da-48d1868d2b36
Δtt .+= sj_amp_osr .* sin.(phi_sj); 

# ╔═╡ 7a2b8a2d-096e-4fea-b622-88e746981442
begin
drv.Δtt_ext[eachindex(drv.Δtt_prev_nui)] .= drv.Δtt_prev_nui
drv.Δtt_ext[lastindex(drv.Δtt_prev_nui)+1:end] .= Δtt
end;

# ╔═╡ b28ec8ad-a5df-4684-af4e-78eaac7fc506
begin
drv.Vext[eachindex(drv.V_prev_nui)] .= drv.V_prev_nui
drv.Vext[lastindex(drv.V_prev_nui)+1:end] .= Vosr 
	#we will just reuse the bit sequence above
end;

# ╔═╡ 22007bcd-0098-4c70-86d2-4b3d0ddb54b4
drv_jitter_tvec!(drv.tt_Vext, drv.Δtt_ext, param.osr);

# ╔═╡ 371dbb85-68d7-4ddc-bce5-4036d5c32314
itp = linear_interpolation(drv.tt_Vext, drv.Vext); #itp is a function object

# ╔═╡ 0b73cb4e-fe64-45df-9baf-13513c616bec
tt_uniform = (0:param.blk_size_osr-1) .+ drv.prev_nui/2*param.osr;
#note here tt_uniform is shifted by prev_nui/2 to give wiggle room for sampling "before" and "after" the current block. This is necessary for sinusoidal jitter

# ╔═╡ 55ab6dd0-4680-404d-9857-a757f15c6c5b
Vosr_jittered = itp.(tt_uniform); 
#To interpolate, use the itp object like a function and broadcast to a vector

# ╔═╡ a6919863-c211-4524-8cad-3e0845764625
drv_ir = u_gen_ir_rc(param.dt, drv_bw, 100*param.tui);

# ╔═╡ 5f137628-6545-4a08-a6ae-6d5f1bdd1241
Vo_filt_conv = u_conv(Vosr_filt, drv_ir, gain = drv.swing * param.dt);

# ╔═╡ d2a72861-e320-4547-885e-5ebe88473f9e
import .BlkBIST: pam_gen_top!

# ╔═╡ aa59ed17-207c-4a4a-80fd-93ca5999c207
pam_gen_top!(bist) #generate some PRBS bits

# ╔═╡ 3b04da34-0070-4c07-bd86-b7542372aafc
import .WvfmGen: w_gen_eye_simple

# ╔═╡ bcada016-27b2-4103-84bb-7a8015eae67c
eye_tx = w_gen_eye_simple(drv.Vo, 
							x_npts_ui, x_npts, 
							y_range, y_npts, 
							osr= param.osr, x_ofst=0);

# ╔═╡ 3c23ca8d-1202-4fb5-a082-3c1a76f1793e
begin
feye = Figure()
heatmap!(Axis(feye[1,1]), x_grid, y_grid, eye_tx, 
            colormap=:turbo, #try :inferno, :hot, :viridis
        )
feye
end


# ╔═╡ 6010b512-1b80-492a-a721-b6a1f2687a34
eye_tx_filt = w_gen_eye_simple(Vo_filt_conv[1:param.blk_size_osr], 
								x_npts_ui, x_npts, 
								y_range, y_npts, 
								osr= param.osr, x_ofst=0);

# ╔═╡ fe879f02-ce99-4745-be08-4a152e2a1d86
begin
feye_filt = Figure()
heatmap!(Axis(feye_filt[1,1]), x_grid, y_grid, eye_tx_filt, 
            colormap=:turbo, #try :inferno, :hot, :viridis
        )
feye_filt
end

# ╔═╡ cd5e4ece-3614-469c-9193-72b2d0590444
begin
#unit conversion
rj_osr1 = rj_s/param.tui*param.osr
sj_amp_osr1 = sj_amp_ui*param.osr
sj_freq_norm1 = sj_freq*param.tui

Δtt1 = zeros(param.blk_size);
Δtt1[1:2:end] .+= dcd/2*param.osr;
Δtt1[2:2:end] .-= dcd/2*param.osr;  #add dcd
Δtt1 .+= rj_osr1 .* randn(param.blk_size); #add rj
phi_sj1 = (0.0 .+ (2π*sj_freq_norm1) * (1:param.blk_size)) .% (2π);
Δtt1 .+= sj_amp_osr1 .* sin.(phi_sj1); #add sj

#gen jittered time
drv.Δtt_ext[eachindex(drv.Δtt_prev_nui)] .= drv.Δtt_prev_nui
drv.Δtt_ext[lastindex(drv.Δtt_prev_nui)+1:end] .= Δtt1
drv_jitter_tvec!(drv.tt_Vext, drv.Δtt_ext, param.osr);
#voltage
drv.Vext[eachindex(drv.V_prev_nui)] .= drv.V_prev_nui
drv.Vext[lastindex(drv.V_prev_nui)+1:end] .= Vosr 

#interpolate and remap to simulation grid. 
#compare and see the speed difference between the custom interpolation and the built-in one when you know the problem
Vosr_jit1 = Vector{Float64}(undef, length(tt_uniform))
# @time begin
itp1 = linear_interpolation(drv.tt_Vext, drv.Vext); #itp is a function object
Vosr_jit1 = itp1.(tt_uniform); 
# end
# @time drv_interp_jitter!(Vosr_jit1, drv.tt_Vext, drv.Vext, tt_uniform)
	
#convolve with ir
ir_high_bw = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui);
Vo_jit1 = u_conv(Vosr_jit1, ir_high_bw, gain = drv.swing * param.dt);
Vo_jit1_trunc = @views Vo_jit1[100*param.osr:end-100*param.osr]; #trim off some garbarge for now

eye_tx_jit = w_gen_eye_simple(Vo_jit1_trunc,
						x_npts_ui*2, x_npts_ui*2, 
						y_range, y_npts, 
						osr= param.osr, x_ofst=round(2*x_npts_ui/3));
#To see duty cycle effect better, here the eye diagram is plotted over just 1UI
end;

# ╔═╡ 9c3a3512-3ce2-45a0-a645-899eeed433f2
begin
cm = cgrad([:white, :blue, :red], [0.00001, 0.05], scale=:exp) 
#custom colormap for better visualization of jitter
feye_jit = Figure()
heatmap!(Axis(feye_jit[1,1]), x_grid, y_grid, eye_tx_jit, 
            colormap=cm, #try :inferno, :hot, :viridis
        )
feye_jit
end

# ╔═╡ 72e099af-d0b6-45c0-8106-7746eb2a3a2d
jit_0x = mod.(u_find_0x(Vo_jit1_trunc), param.osr) ./ param.osr;
#find zero crossing and normalize to within 1UI

# ╔═╡ 34895f52-6a85-4f03-9177-ad8432df125f
begin
f0x = Figure()
lines!(Axis(f0x[1,1]), jit_0x)
f0x
end

# ╔═╡ 09de39c9-7290-4cbc-ab8b-1343ab68abe1
jit_0x_unwrap = u_unwrap_0x(jit_0x);

# ╔═╡ 8b91cf27-e7e8-48a9-9104-3370d1095a4c
begin
f0x_uw = Figure()
lines!(Axis(f0x_uw[1,1]), jit_0x_unwrap)
f0x_uw
end

# ╔═╡ f8abc19c-ce00-43b6-b1b9-b34151c27d0b
begin
jitter_bnd = (-0.05,0.05) .+ (-6, 6).*(rj_s/param.tui) .+ (-1.2, 1.2).*sj_amp_ui
fd = Figure()
	density!(Axis(fd[1,1]), jit_0x_unwrap .- mean(jit_0x_unwrap), boundary=jitter_bnd, npoints=200)
fd
end

# ╔═╡ 53870af8-fd88-4f1a-b053-c409f92dfe98
md"""
## Appendix
"""

# ╔═╡ 2f90f15c-8bfd-454c-b90a-143bf133bf27
md"""
```Julia
@kwdef mutable struct Param
    const data_rate::Float64
    const osr::Int64
    const pam::Int8 = 2
    const bits_per_sym::Int8 = Int(log2(pam))
    const fbaud = data_rate/bits_per_sym
    const fnyq= fbaud/2
    const tui = 1/fbaud
    const dt= tui/osr

    const blk_size::Int64
    const blk_size_osr::Int64 = blk_size*osr
    const subblk_size::Int64 = blk_size
    const nsubblk::Int64 = Int(blk_size/subblk_size)
    const nsym_total::Int64
    const nblk = Int(round(nsym_total/blk_size))

    const rand_seed = 300  

    cur_blk = 0
    cur_subblk = 0

end
```
"""

# ╔═╡ 32513a0a-92cb-4c57-8a82-5463c0f37234
md"""
```Julia
@kwdef mutable struct Bist
    const param::Param
    const polynomial::Vector{UInt8}
    const order::UInt8 = maximum(polynomial)
    const inv = false

    gen_seed = ones(Bool,order)
    gen_gray_map::Vector{UInt8} = []
    gen_en_precode = false
    gen_precode_prev_sym = 0

    chk_start_blk = 100
    chk_seed = zeros(Bool,order)
    chk_precode_prev_sym = 0
    chk_lock_status = false
    chk_lock_cnt = 0
    chk_lock_cnt_threshold = 128 
    ber_err_cnt = 0
    ber_bit_cnt = 0

    So_bits::Vector = zeros(Bool, param.bits_per_sym*param.blk_size)
    So::Vector = zeros(param.blk_size)
    Si = CircularBuffer{UInt8}(param.blk_size)
    Si_bits::Vector = zeros(Bool, param.bits_per_sym*param.blk_size)
end
```
"""

# ╔═╡ bc91a803-b5ac-4dc2-aadd-f65e865c26ae
md"""
```Julia
function u_gen_ir_rc(dt,bw,t_len)
    tt = [0:dt:t_len-dt;]

    ω = (2*π*bw)
    ir = ω*exp.(-tt*ω)
    ir = ir/sum(ir*dt)

    return ir
end
```
"""

# ╔═╡ feca2e5a-5839-44e5-bf53-740d586bb8d0
function drv_interp_jitter!(vo, tt_jitter, vi, tt_uniform)
    last_idx = 1
    for n = eachindex(tt_uniform)
        t = tt_uniform[n]
        for m = last_idx:lastindex(tt_jitter)-1
            if (t >= tt_jitter[m]) && (t < tt_jitter[m+1])
                k = (vi[m+1]-vi[m])/(tt_jitter[m+1]-tt_jitter[m])
                vo[n] = vi[m] + k*(t-tt_jitter[m])
                last_idx = m
                break
            end
        end
    end

    return nothing
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
BenchmarkTools = "~1.5.0"
DSP = "~0.7.9"
DataStructures = "~0.18.20"
Distributions = "~0.25.109"
FFTW = "~1.8.0"
GLMakie = "~0.10.0"
Interpolations = "~0.15.1"
MAT = "~0.10.6"
Makie = "~0.21.0"
Parameters = "~0.12.3"
PlutoUI = "~0.7.59"
StatsBase = "~0.34.3"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "b7da7bf1c081688b2d963514112c93746db8e83d"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

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

[[deps.BufferedStreams]]
git-tree-sha1 = "4ae47f9a4b1dc19897d3743ff13685925c5202ec"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.2.1"

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

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

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

[[deps.GLFW]]
deps = ["GLFW_jll"]
git-tree-sha1 = "35dbc482f0967d8dceaa7ce007d16f9064072166"
uuid = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
version = "3.4.1"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GLMakie]]
deps = ["ColorTypes", "Colors", "FileIO", "FixedPointNumbers", "FreeTypeAbstraction", "GLFW", "GeometryBasics", "LinearAlgebra", "Makie", "Markdown", "MeshIO", "ModernGL", "Observables", "PrecompileTools", "Printf", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "21127b3106d8b97048975134008c143a8b8b33c5"
uuid = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
version = "0.10.0"

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

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "82a471768b513dc39e471540fdadc84ff80ff997"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.3+3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ca0f6bf568b4bfc807e7537f081c81e35ceca114"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.10.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

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

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

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

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "ed1cf0a322d78cee07718bed5fd945e2218c35a1"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.6"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "4099bb6809ac109bfc17d521dad33763bcf026b7"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.2.1+1"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "ce0ca3dd147c43de175c5aff161315a424f4b8ac"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.3.3+1"

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

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8c26ab950860dfca6767f2bbd90fdf1e8ddc678b"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.11"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f12a29c4400ba812841c6ace3f4efbb6dbb3ba01"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModernGL]]
deps = ["Libdl"]
git-tree-sha1 = "b76ea40b5c0f45790ae09492712dd326208c28b2"
uuid = "66fc600b-dfda-50eb-8b99-91cfa97b1301"
version = "1.1.7"

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

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e25c1778a98e34219a00455d6e4384e017ea9762"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.6+0"

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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

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

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

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

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "46bf7be2917b59b761247be3f317ddf75e50e997"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.2+0"

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
# ╟─248b7a30-25df-11ef-11c5-d33b88daa42b
# ╟─a934351d-c8c8-475e-a7de-5dee9e45d549
# ╟─938cc465-1b5a-4e02-9956-a57cb9e04ccb
# ╟─42718f5c-44a5-4d0f-b6f0-b6db8a499d14
# ╟─01820462-23f1-498e-95ac-3843101fe4eb
# ╟─ef180811-bb81-4b2d-955f-3785308a9bc9
# ╟─bdf79336-77db-491b-ba6d-049b33f088b6
# ╟─fcbd062f-d907-4691-9365-7fa3fe68790d
# ╟─b09b12d2-e7fb-4944-a680-a5f676176390
# ╟─f371095b-daa4-4b57-9c6d-274d32f801bd
# ╠═916ae0ce-04e6-4d78-9521-ff5501725bd2
# ╟─1ae615a4-7a81-4abb-b3e7-2435d8ea2616
# ╠═efa51aac-75c3-4621-bf04-ff42a87f5bd4
# ╟─5b16cc5d-6de8-4126-ab11-bfb0cabf44b7
# ╟─f0209a9a-ddfd-419f-85ae-fecf4882d7ab
# ╟─0f373aa6-55c9-4e48-abf0-019373f36a9b
# ╟─0e155843-85fe-4573-afce-0fc11f56d3e8
# ╟─e36746f1-19f4-4ecf-a6a0-29be5c307e26
# ╟─4e8daf25-af1c-4ac8-8436-7d971fd9c004
# ╟─cbcd9551-da4c-4c91-bec6-a698f5ec6fa5
# ╟─f0495ec1-04f5-4d72-82f9-35c768927e8e
# ╟─bb49149b-db9c-43d2-bf17-6862c26ebc11
# ╟─34319afe-4569-4cf1-9df4-21f976388d1c
# ╟─19c8f09b-e1f0-4117-8c6e-3623988b05fa
# ╟─c109ab44-1cc8-41f2-8f9b-2c04dbf824b7
# ╠═e97cfdc6-c2d8-4f68-bb0d-39117b02d560
# ╟─6a733af2-63f2-418c-8100-57195d9d790a
# ╟─eb15958d-acd1-484f-972c-b51f49352dfd
# ╟─111848d3-e97a-4869-807c-c9f927d2ca58
# ╟─4221ddc6-5958-4e75-9052-55929f7f6d74
# ╠═e103e442-69d3-44c6-8670-053534322412
# ╟─022cf285-2a25-407a-a322-5cbab294aeba
# ╟─6ac67d2d-79ae-4f9e-9328-c5a5981e7680
# ╠═3adc12e0-1ec7-4ea2-bc2a-dc35a987d96f
# ╟─f389a6fd-c388-4e9c-8d66-37c00305c1a6
# ╠═aa59ed17-207c-4a4a-80fd-93ca5999c207
# ╠═7644dc88-7529-4a0e-8646-b024dad540b3
# ╟─bb94d0a1-48f5-401e-b12a-2378435d4422
# ╠═7b13b871-3b0b-42c3-b9a7-00cf246c86b1
# ╟─e653c634-02db-44b3-8501-e2614ee2c07c
# ╟─2ab349bc-07af-4056-a729-37db20705389
# ╠═5062c917-f265-4e63-98a3-544bb9482bc1
# ╠═bcada016-27b2-4103-84bb-7a8015eae67c
# ╟─3c23ca8d-1202-4fb5-a082-3c1a76f1793e
# ╟─85f612fa-1889-4631-a334-27530881d7b8
# ╠═91469577-cb17-4505-bdcd-bafce396d564
# ╠═4651d7ef-3d52-4275-a0ce-13051d1a1cdf
# ╟─03ec64bb-7bb7-4e7f-9b4e-00ec713d739c
# ╟─f96baca0-637c-4b02-855d-d758121fa96e
# ╠═a6919863-c211-4524-8cad-3e0845764625
# ╠═01549352-1220-4e5b-ae1f-143a5bcada0c
# ╠═bfee37c6-fa4c-4e08-9b76-cfa298fff745
# ╠═2a9645ac-a66e-4133-8dc7-2b455e0473f5
# ╠═0984e74f-c480-4aac-b429-fae9d10afd6d
# ╠═5f137628-6545-4a08-a6ae-6d5f1bdd1241
# ╠═6010b512-1b80-492a-a721-b6a1f2687a34
# ╟─fe879f02-ce99-4745-be08-4a152e2a1d86
# ╟─c67fbbff-19e3-4ffd-bafd-d934e4882268
# ╟─dc08b7fa-f12a-4aa1-ac9d-049583cb91c1
# ╟─98fb5603-fe2e-455d-8d28-f96a8113121d
# ╟─55a76e1a-b520-4863-8864-ce2e772742ee
# ╟─24f1ef71-2d2f-42ff-a018-e2922a185de4
# ╠═d3f0e565-51c3-497d-bd59-b0acb842e402
# ╟─93413781-bc43-45d2-8d35-0017b948b189
# ╠═e7e4e516-f575-4403-ad0a-cf0511d549b0
# ╠═f549409a-4ac7-4244-8d12-44d1eabe566b
# ╟─74ffbe27-ae1b-4dac-84b8-83f074038a2e
# ╠═b0c39645-22ec-49fb-ae84-900b169fdde7
# ╠═b1ef7a3d-1eb3-4703-9827-8ede9165aab9
# ╟─22373810-8aa9-4950-9c41-e1cf08372421
# ╠═a2c309cc-9fea-4ab0-909e-85f17ffa4e10
# ╠═8a9eefca-158b-400f-a7fa-2fcf18c02324
# ╠═a5ec92b5-a0ee-45b4-b017-4822863b2055
# ╠═585edc1c-3349-4888-a8da-48d1868d2b36
# ╟─4c28d8ef-9169-4c9a-8a03-bc0f67e52991
# ╟─6be5384f-b06e-494e-a659-f772438e25d6
# ╠═7a2b8a2d-096e-4fea-b622-88e746981442
# ╟─510dec4a-07eb-4fa7-addb-876e3c171587
# ╠═b28ec8ad-a5df-4684-af4e-78eaac7fc506
# ╟─c66d9219-4f70-4917-855f-b3ff1c1c73f4
# ╠═63b400e2-c7d2-4873-ad93-6a23d14b760b
# ╟─b94044fa-193a-4bb9-998b-02e74bbee98c
# ╠═22007bcd-0098-4c70-86d2-4b3d0ddb54b4
# ╟─f063e527-2363-4692-b79a-9646a9b68723
# ╟─14ce9a21-09de-4ec6-9511-5a4281459032
# ╠═371dbb85-68d7-4ddc-bce5-4036d5c32314
# ╟─cebb0066-edf7-41c2-99ea-bd5c966fe04b
# ╠═0b73cb4e-fe64-45df-9baf-13513c616bec
# ╠═55ab6dd0-4680-404d-9857-a757f15c6c5b
# ╟─a9785c06-83d2-46a4-aec1-63d0644c19df
# ╟─d8da233c-81ef-4916-88e3-2b9b0fac5a57
# ╠═cd5e4ece-3614-469c-9193-72b2d0590444
# ╟─690ea232-0ef4-45ab-a158-6ff355d8101f
# ╟─9c3a3512-3ce2-45a0-a645-899eeed433f2
# ╟─0c05c7de-ff91-4927-8f8e-44340f014d55
# ╟─6d5ddcaf-bab8-44d8-a2f9-9a90bbe25fdf
# ╟─a296059d-7228-48da-b538-cac3e0c06ab7
# ╠═72e099af-d0b6-45c0-8106-7746eb2a3a2d
# ╟─34895f52-6a85-4f03-9177-ad8432df125f
# ╟─d130e739-b6d2-445e-ad27-510fb34da8bc
# ╠═09de39c9-7290-4cbc-ab8b-1343ab68abe1
# ╟─8b91cf27-e7e8-48a9-9104-3370d1095a4c
# ╟─24be3bdb-c3db-42ab-b2ef-c313b7108ff8
# ╠═f8abc19c-ce00-43b6-b1b9-b34151c27d0b
# ╟─78de21b5-984b-4e1f-9759-dfef88885c48
# ╟─1dfc2b52-1749-4782-9c5b-7c3289069803
# ╠═bb39273b-820f-472a-a6ae-c6ed9094044f
# ╟─bb408a87-8ded-403e-87a8-2a9cf42e79fb
# ╟─555d708b-ab85-4345-9e90-d6d5f480cbb1
# ╟─4deb9728-b83e-4f7e-8fd9-2cdf215c6408
# ╠═3bf00ba7-1593-4aad-bd92-49095c1f65e5
# ╟─8f763793-3e0c-4d72-988c-e5fe61585135
# ╠═41cb22b3-f000-4326-91ce-e52e95048cf0
# ╟─724709c6-1deb-4dd9-8d83-6db16a9bd2b6
# ╠═0f491d87-d14a-4bfe-8a2d-821797926590
# ╠═031f8558-b375-469f-8c62-c7d92a75441f
# ╠═610b7b53-bd09-4a02-8e60-f43209da4128
# ╠═ae7a47ae-2f31-426f-88ff-1cce13ff9622
# ╠═d2a72861-e320-4547-885e-5ebe88473f9e
# ╠═9c951989-97db-4ebd-9cd6-233beec23664
# ╠═3b04da34-0070-4c07-bd86-b7542372aafc
# ╠═907c712d-8fc0-49ac-9c72-ae53f8c78794
# ╠═1fcec875-7dda-4afe-bd9e-a9fea816fc32
# ╠═1417102e-62a3-4999-a278-9314e24a5edb
# ╠═a2334355-f50c-456f-8a51-578885f07528
# ╟─53870af8-fd88-4f1a-b053-c409f92dfe98
# ╟─2f90f15c-8bfd-454c-b90a-143bf133bf27
# ╟─32513a0a-92cb-4c57-8a82-5463c0f37234
# ╟─bc91a803-b5ac-4dc2-aadd-f65e865c26ae
# ╠═feca2e5a-5839-44e5-bf53-740d586bb8d0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
