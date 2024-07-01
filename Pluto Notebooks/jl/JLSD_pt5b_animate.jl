### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 8e12dc88-453c-4789-9079-3e3f7c789b9f
using Parameters, UnPack, DSP, FFTW, Random, Interpolations, DataStructures, Distributions, StatsBase, MAT

# ╔═╡ 33e5faec-0407-489e-a9e7-d90ce3ee9368
using GLMakie, Makie, ColorSchemes, Colors

# ╔═╡ 81b6950c-9eb2-4a72-820c-e8a246776c56
include("../../src/structs/TrxStruct.jl");

# ╔═╡ e032a9c9-d5dc-42eb-a321-0b70526e8e67
begin
include("../../src/blks/BlkBIST.jl");
import .BlkBIST: pam_gen_top!
end

# ╔═╡ 3d901245-bae6-4664-bfbe-cc48c9ed0004
begin
include("../../src/blks/BlkTX.jl");
import .BlkTX: dac_drv_top!
end

# ╔═╡ 656909c7-8a50-416c-a444-881d469c2a08
begin
include("../../src/blks/BlkCH.jl");
import .BlkCH: ch_top!
end

# ╔═╡ e35022f4-ee62-46c1-9e22-1171429c2572
begin
include("../../src/util/Util_JLSD.jl");
import .Util_JLSD: u_gen_ir_rc, u_fr_to_imp
end

# ╔═╡ 8dc41adf-b61d-49da-a04d-54442a051a67
begin
include("../../src/blks/WvfmGen.jl");
import .WvfmGen: w_gen_eye_simple, w_gen_eye_simple!
end

# ╔═╡ 69e46020-2f48-11ef-0e3a-23d26a279b9d
md"""
# Building SerDes in Julia pt. 5B - Plotting w/ Animation
"""

# ╔═╡ 1ce6e5d5-7deb-49a2-945a-875619ccb86c
md"""
![pam4eye](https://circuit-artists.com/wp-content/uploads/2024/06/pam4eye_animation_new.gif)
"""

# ╔═╡ bf3c796b-bcd8-493d-b036-de79991f6271
md"""
Finally we are here! This notebook will walk through how to build the animation and widget window shown above using Makie alone! The whole concept rests upon the idea of ```Observables``` - pointer like data structures that allow event listeners. Actually, you have already experienced the power of observables. This very Pluto notebook has observables in the background to achieve reactivity! Makie's own internal figure objects use observables too. Without further ado, let's begin building
"""

# ╔═╡ f0a8f5d2-572c-4e0e-b863-b99848cc1a1a
md"""
## So what are observables?
"""

# ╔═╡ 4decb1c3-7c59-4321-91eb-39d2395a093c
md"""
The idea of ```Observables``` is not that complicated. Just look at how simple the [official documentation](https://juliagizmos.github.io/Observables.jl/dev/) is compared to other Julia packages. The developers have put in quite some efforts to make Julia's observables easy to use and understand.
"""

# ╔═╡ e9edf872-886f-4937-8f00-6ed203630a06
md"""
At a high level, an observable is like a pointer (or ```Ref``` in Julia) that contains the data of interest but can be watched with handlers. We create an observable variable like below
"""

# ╔═╡ 5b78e003-edbd-4791-9fd1-48d0b4523bd4
ob_name = Observable("Iron Man");

# ╔═╡ b5d11536-5508-40b5-886c-8bbd6cbc73cc
md"""
The argument to the ```Observable()``` constructor provides both an initial value and determines the type of the observable variable. So if I want to assign an interger to ```ob_name```, it throw an error. There are two ways to access/update an observable's value - w/ ```.val``` or ```[]```.
"""

# ╔═╡ 7973422d-fbe9-4049-8325-a620520c9be8
begin #these are two ways to access the data in an observable
println(ob_name)
println(ob_name.val); 
println(ob_name[]);
end

# ╔═╡ d13dedd3-84c6-484e-b77f-98601d908474
md"""
The difference is that only using ```[]``` will trigger the listner event when updating the observable value.
"""

# ╔═╡ 9fd804d6-1bbd-4697-9cef-c3f1113a1967
begin
ob_name.val = "Hulk" #no event is triggered, but value still updated
println(ob_name[])
end

# ╔═╡ 303b9ac4-a7de-4d4d-89a7-aab55697cad0
begin
ob_name[] = "Thor" #an update event will triggered and handler function will execute if defined
println(ob_name[])
end

# ╔═╡ 069c9012-f240-464e-95ed-9e1338b2d1eb
ob_name[] = 1 #ERROR, ob_name can only take String now

# ╔═╡ 64a9d1c8-ae20-4731-b339-503ba5bbee60
md"""
Event handlers are defined with the ```on``` function. Here is a simple example
"""

# ╔═╡ 688b1610-ff94-4f9f-9d9e-49df00f9977e
begin

food_ob = Observable("pizza")
eat_func = on(food_ob) do food #every time the observable changes do the following
	println(uppercasefirst(food) * " time!")
end
	
end;

# ╔═╡ 55836ef7-0309-4bc1-8f43-9e5c28fe7100
food_ob[] = "coffee"; #try changing it here and see the handler function work

# ╔═╡ 10bfdc04-14b7-423a-8191-2a3324b752db
food_ob.val = "Fried chicken"; #no printout if using .val to update, no matter how bad you want fried chicken

# ╔═╡ 43a01dc0-cb9c-457f-ba87-2b8aaf81723a
md"""
So how does this translate to Makie plots? Using observables, we only need to update the data and Makie will update the plot automatically! Let's build a simple plot
"""

# ╔═╡ d3f53e27-5958-48d0-82ad-884500c60f15
# ╠═╡ disabled = true
#=╠═╡
begin
ftest = Figure();
axtest = Axis(ftest[1,1])
test_data = Observable(randn(10));
lines!(axtest, test_data);
ylims!(axtest,(-3.5, 3.5))
test_lbl = Label(ftest[2,1], "Rerun me!", tellwidth = false, fontsize=32);
display(GLMakie.Screen(), ftest); #run this cell to get a pop-up window
end
  ╠═╡ =#

# ╔═╡ 1422f9a9-024e-4715-b3cc-6f84961a8eac
# ╠═╡ disabled = true
#=╠═╡
begin #Rerun this cell to see the plot change by simply update test_data[]
test_lbl.text = "Running"
for n = 1:5
	test_data[] = randn(10)
	sleep(1)
end
test_lbl.text = "Rerun me!"
end
  ╠═╡ =#

# ╔═╡ d0fcb7ce-2d6b-453d-b9c7-aed81595aba3
md"""
Congratualtions! We just created our first animation with Makie! Using observable completely separates the figure configuration and the actual data pipeline. Here we simply kept updating the observable data, and Makie replots immediately (Makie figure is essentially an event handler on the data itself).
"""

# ╔═╡ 37a3ef57-eace-48e7-b170-e8bebad356b7
md"""
A little rant before I move on: this "paradigm" is way better than what MATLAB has, which is a ```drawnow``` command after redrawing all the plots over and over again. I find myself spending more time on the plotting code than actually working on the model sometimes. Makie + Observable should save people a lot of time when data has to be constantly updated and visualized.
"""

# ╔═╡ 1985a856-9479-445b-a0c1-d2c7117e6e98
md"""
### Chaining observables
"""

# ╔═╡ 6099bac2-6aa1-40b8-8f82-555d08ca1a84
md"""
Observables can also be chained together and each event can trigger more complicated processing down the line. We use ```lift``` to achieve this. 
"""

# ╔═╡ 7987b2c5-cc99-4676-ac22-ca7a4b96258a
ob1 = Observable(2.0);

# ╔═╡ 986a900f-738e-43df-bf10-bcf292660ca0
ob2 = lift((x)->x^2, ob1);

# ╔═╡ 55367577-43fa-4e72-849a-21a584590808
md"""
```ob2``` is now a new observable that's chained with ```ob1```, and we can create a event listener on ```ob2```. Then, if we update ```ob1```, the event handler on ```ob2``` will be triggered as well.
"""

# ╔═╡ 64c399e2-d200-49df-b633-ae74a0b23883
ob2_func = on(ob2) do val
	println(val-10);
end;

# ╔═╡ b396f7fe-e1b0-4d46-8acb-a4fd0e1eda0b
md"""
For more complicated chaining functions, we can either use ```lift() do``` or the ```@lift``` macro. Note w/ ```@lift```, we can directly refer another observable with the ```$``` prefix.
"""

# ╔═╡ ea1bfbb5-471a-410f-89b3-72c031705f8f
ob2_v2 = lift(ob1) do val
	if val > 4
		return 4*val;
	else
		return val^2;
	end
end;

# ╔═╡ 50026a47-7066-4e65-9eaf-40c53c59ca28
ob2_v3 = @lift begin
	return $ob1 > 4 ? 4*$ob1 : $ob1^2
end;

# ╔═╡ 4434df8c-fc34-4ffd-9c53-2b45df07f84a
ob1[] = 2.0; 

# ╔═╡ 3923a75e-8c3a-4ab0-8508-1fdc17d5a853
md"""
## Makie GridLayout and UI blocks
"""

# ╔═╡ 93b1f0ed-9de3-499e-99ad-87cb0f777c9d
md"""
Remember in [part 5a](https://github.com/kevjzheng/JLSD/tree/main/Pluto%20Notebooks/jl) I complained about Makie's syntax not being suited for quick plotting and debug, and created my own ```Figz``` wrapper? Well, as it turned out Makie has some good reasons to have separate ```GridLayout``` and ```Axis``` systems. Makie supports many interactive UI elements directly in the figure window, which means the user might want to create complex layout like a GUI window. That's why not everything has an ```Axis```, but they will all have associated grids.
"""

# ╔═╡ 20de9f86-3fec-4944-9bc3-72adfd47dd60
md"""
You can arbitrarily create grids and subgrids with ```[]``` indexing on the figure object. This grid system allows us to do very flexible UI interface when it comes to both plots and controls. It's super easy in Makie to create grids with different sizes.
"""

# ╔═╡ 6d447a7f-6816-4164-8eb2-2e7e85f7ab00
begin
fgrid = Figure();
#upper left 
Box(fgrid[1,1])
#lower left
Box(fgrid[2,1])
#lower right into left and right
gridlr = GridLayout(fgrid[2,2]) #instead of keep indexing, a Gridlayout object can be created
Box(gridlr[1,1])
Box(gridlr[1,2])
#upper right into even smaller grids
gridur = GridLayout(fgrid[1,2])
sub_gridur = GridLayout(gridur[1,1])
Box(fgrid[1,2][1,1][1,1]) #either keep indexing
Box(gridur[1,1][1,2]) #or use grids/subgrids
Box(gridur[1,1][2,1])
Box(sub_gridur[2,2])
Box(gridur[1,2])
Box(gridur[2,1])
Box(gridur[2,2])
fgrid
end

# ╔═╡ 72af7e60-3398-44e3-ba85-28376c50dd3b
md"""
Makie also provide built in basic interactive UI blocks, including ```Button```, ```Slider```, ```Menu```, ```Textbox```, and ```Toggle```.
"""

# ╔═╡ 0ee12efa-a36f-4638-88ab-2ed6094c6a81
begin
fui = Figure(size=(500,100));

Label(fui[1,1], "Some text")
Button(fui[1,2], label="Click me")
Toggle(fui[1,3], active = true)
Slider(fui[1,4], range = 0:10, startvalue = 5, width = 100)
Menu(fui[1,5], options =["coffee", "tea", "more coffee"], default="coffee")

fui
end

# ╔═╡ 11c65bc7-9ca9-4dc8-9ad8-0f630fc7ed1e
md"""
These grant us enough flexibility and options to begin building the eye diagram widget! To learn more about Makie, go through their well written [tutorials](https://docs.makie.org/stable/tutorials/getting-started), and find more examples at [BeatifulMakie](https://beautiful.makie.org/dev/), where you can copy/paste code and generate things like this right away:
"""

# ╔═╡ 96d45520-b93e-4864-a950-8002afe214a5
# ╠═╡ disabled = true
#=╠═╡
begin
function peaks(; n = 49)
    x = LinRange(-3, 3, n)
    y = LinRange(-3, 3, n)
    a = 3 * (1 .- x').^2 .* exp.(-(x'.^2) .- (y .+ 1).^2)
    b = 10 * (x' / 5 .- x'.^3 .- y.^5) .* exp.(-x'.^2 .- y.^2)
    c = 1 / 3 * exp.(-(x' .+ 1).^2 .- y.^2)
    return (x, y, a .- b .- c)
end
x, y, z = peaks()
with_theme(theme_dark()) do
    fig = Figure(size = (1200,800))
    ax1 = Axis(fig[1,1], aspect = 1)
    ax2 = Axis3(fig[1,2]; aspect = (1,1,0.7), perspectiveness = 0.5)
    axs = [ax1, ax2]
    cmap =  :diverging_bkr_55_10_c35_n256
    contourf!(axs[1], x, y, z; levels = 0.1:0.08:0.9, mode = :relative,
        colormap = cmap)
    # bug, colormap cannot be transparent
    contourf!(axs[2], x, y, z; levels = 14, colormap = cmap)
    contour3d!(axs[2], x, y, z; levels = 14, colormap = cmap,
        transparency = true, linewidth = 5)
    limits!(axs[1], -3,3,-3,3)
    hidedecorations!.(axs; grid = false)
    fig
end
end
  ╠═╡ =#

# ╔═╡ 22e655f4-d8e9-4900-affa-546ab5026e81
md"""
## It's alive!
"""

# ╔═╡ 44516906-4475-42bb-a3db-c1de5a2ea04c
md"""
Let's recap our simulation framework and generate some signals first before getting carried away with pretty plots.
"""

# ╔═╡ 8f93c1fb-c5bf-4ba7-8ae7-4923d10c57b4
md"""
Our simulation framework uses structs and mutating functions to model each circuit block's parameters and operations. In this tutorial, we will use the TX model (details in [part 4](https://github.com/kevjzheng/JLSD/tree/main/Pluto%20Notebooks/jl)) and add a toy channel. The channel will include both a PCB trace (translated to impulse response from .s4p frequency response) and another 1st order pad response. The additive Gaussian noise level is specified in dBm/Hz.
"""

# ╔═╡ 70d30065-2dde-446f-93f3-e86bc81d438f
begin
param = TrxStruct.Param(  
			data_rate = 56e9,
			pam = 4,
			osr = 20,
			blk_size = 2^10,
			subblk_size = 32, 
			nsym_total = Int(1e6));

bist = TrxStruct.Bist(  
			param = param,
			polynomial = [28,31]);

drv = TrxStruct.Drv(
		param = param,
		ir = u_gen_ir_rc(param.dt, param.fbaud/2, 20*param.tui), #1st order ir
		fir = [1.0, -0.0],
		swing = 0.8,
		jitter_en = true,
		dcd = 0.00,
		rj_s = 000e-15,
		sj_amp_ui = .0,
		sj_freq = 2e5);
	
ch = TrxStruct.Ch(
			param = param,
			ch_en = false,
			ir_ch = u_fr_to_imp("../../channel_data/TF_data/channel_4inch.mat", 
					param.tui, param.osr, npre = 20, npost= 79),
			ir_pad = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui),
			noise_en = true,
			noise_dbm_hz = -140)

trx = (;param, bist, drv, ch);
end;

# ╔═╡ 2f36a08e-3aeb-4dbd-a341-98e8f418ccb2
md"""
Now that our TX + channel isntances are instantiated, we will create a ```step_sim_blk``` function, which basically runs one block of the simulation and we will use the generated waveform to create one frame of the eye diagram. Note that this ```step_sim_blk``` function is taking an ```ch_en``` argument which allows us to toggle between channel included or excluded.

At a high level, if block size is too small, the frame rate will naturally be faster but noisier per frame. If the block size is too big, frame rate drops but each frame gets more data for plot. ~1000 symbols per block seem to be a good tradeoff. Feel free to play around with these paramters. To make it more fun, we will do a PAM4 eye this time. 
"""

# ╔═╡ f18af4fc-1bfc-45f3-88af-5feb27303155
function step_sim_blk(trx)
	@unpack bist, drv, ch = trx 
	
	pam_gen_top!(bist)
	
	dac_drv_top!(drv, bist.So)

	ch_top!(ch, drv.Vo)

	return copy(ch.Vo)
	#we are returning copies of drv/ch.Vo because remember they are views/subarrays
	
end

# ╔═╡ c1d734d9-5487-4a77-9bb8-30fdb905798d
md"""
If you followed my previous notebooks, the following code blocks should look familiar and you tweak the eye diagram parameters below.
"""

# ╔═╡ 696c8ef8-863c-4b7e-84e8-3e579cb9255e
begin
x_nui = 2
x_npts_ui = 256
x_npts = x_nui*x_npts_ui
y_range = 1 #+/-0.5
y_npts = 512
x_grid = -x_nui/2: 1/x_npts_ui:x_nui/2-1/x_npts_ui
y_grid = -y_range*(0.5-.5/y_npts):y_range/y_npts:y_range*(0.5-.5/y_npts)
cm_turbo = copy(colorschemes[:turbo].colors)
cm_turbo[1] = RGB(.97,.97,.97)
end;

# ╔═╡ 47803c96-38b7-4554-a5e9-031966756b5f
md"""
Now we will use observable to define the eye diagram to be plotted. ```eye_buffer``` will observe the time domain waveform of either driver or channel output (again, the initial value doesn't really matter, but the type needs to be defined). ```last_heatmap``` and ```shadow_weight``` are used for creating "after images" of eye diagrams - think of it as a shadowing that slowly decays.
"""

# ╔═╡ 4a1aa2b3-8fdd-4bc9-bfed-8bd546a24bef
eye_buffer = Observable(ch.Vo);

# ╔═╡ ac5ad027-a8cf-44a7-86ed-c19a589019ce
cur_heatmap = zeros(x_npts,y_npts);

# ╔═╡ 2cadbaef-eebb-4943-a08d-d70f62dc4c72
shadow_weight = Observable(0.5);

# ╔═╡ c51fbd6a-b07c-4ddb-927a-fa19186db679
xofst = Observable(0);

# ╔═╡ b4434801-fabf-4152-a136-06c66a7abe9c
md"""
Here we create a new observable, chained from ```eye_buffer```. So eveytime ```eye_buffer``` is updated, this new observable ```eye_heatmap``` will be update automatically. We will be fancy here and put the eye diagram through a leaky integrator (i.e. ```eye_out[n] = eye_in[n] + α × eye_out[n-1], α∈[0,1)```) to create the after image effect.
"""

# ╔═╡ 78b50e9d-e650-459b-a291-e33235974534
eye_heatmap = @lift begin
	w_gen_eye_simple!(cur_heatmap, $eye_buffer, x_npts_ui, x_npts, y_range, y_npts; osr = trx.param.osr, x_ofst=xofst[], shadow = shadow_weight[])
	cur_heatmap #returns the current heatmap
end;

# ╔═╡ 4dd5297b-0dd5-481c-b851-d9b4e009271c
md"""
Time to make the widget UI! The entire figure creation is capture in the following code block. The widget includes a run button to start/stop the simulation, a button for enabling/disabling channel (i.e. ISI).  A set of sliders for tuning TX/channel parameters real time. And a small slider for setting the after image weight.
"""

# ╔═╡ f269e7e6-42f0-4f49-a28a-252d1f572f8f
begin #run this block to see the widget!
f = Figure(size=(700,700))
#plot the eye diagram, pass the observable directly in
sl_xofst = Slider(f[1,1:3], range = -x_npts/2:x_npts/2-1, startvalue = 0)
heatmap!(Axis(f[2:3, 1:3]), x_grid, y_grid, eye_heatmap, colormap=cm_turbo)

gctrl = GridLayout(f[4,1:3])
#run button
btn_run = Button(gctrl[1,1], label ="Run", tellwidth=false, tellheight=false, width = 140, height=60, fontsize=28);
#enable channel button
btn_ch_en = Button(gctrl[2,1], label="Channel\nDisabled", tellwidth = false, tellheight=false, width=100, height=50, fontsize=18);
#slider for after image weight
sl_shadow = SliderGrid(gctrl[3,1], (label="After Image",  range=0:10, startvalue=shadow_weight[]), tellwidth = false, tellheight=false, width=200);
#sliders for TX/Channel parameters
sl = SliderGrid(gctrl[:,2:3],
	(label = "FIR post-tap", range=0:-0.05:-0.5, startvalue = trx.drv.fir[2]),
	(label = "DCD", range=0:1:20, format = "{:d} %", startvalue = trx.drv.dcd*100),
    (label = "RJ", range=0:0.1:2, format = "{:.1f} ps", startvalue = trx.drv.rj_s/1e-12),
	(label = "SJ freq", range=0:50:1000, format = "{:d} kHz", startvalue = trx.drv.sj_freq/1e3),
    (label = "SJ amp", range = 0:0.05:1, format = "{:.2f} UI", startvalue = trx.drv.sj_amp_ui),
	(label = "Channel noise", range=-150:.5:-130, format = "{:.1f} dBm/Hz", startvalue = trx.ch.noise_dbm_hz),
	tellheight=false, width = 300)
#open a new figure window
display(GLMakie.Screen(), f); #uncomment and see the widget in a pop-up window
end;

# ╔═╡ 7254e4f4-31f1-44ab-a9bf-c7d41828b7e2
md"""
When you run the code block above, a window should open up. After you hit run, you should see sometime like this

![widget_default](https://circuit-artists.com/wp-content/uploads/2024/06/eye_widget_default.png)
"""

# ╔═╡ b4eb62d7-3942-4987-a106-ffe305c2a997
md"""
The section below contains the code related to how each UI element works. Each one is simple enough and should be self explanatory
"""

# ╔═╡ 3b9c6cf2-0667-480b-9d50-3a7d9952cc91
begin
ch_en = Observable(false); #this is passed to the step_sim_blk function above
ch_en_func = on(btn_ch_en.clicks) do clicks
	trx.ch.ch_en = ~trx.ch.ch_en;
	btn_ch_en.label = trx.ch.ch_en ? "Channel\nEnabled" : "Channel\nDisabled"
end
end;

# ╔═╡ 8b44925f-d944-4498-a926-11ebb56db7d8
begin
eye_shadow = lift(sl_shadow.sliders[1].value) do val
	shadow_weight[] = 1.0 - .5^val 
end
end;

# ╔═╡ 5a04a0b6-1fdb-49cf-8f49-fb12db7f8597
begin
eye_ofst = lift(sl_xofst.value) do val
	xofst[] = round(Int, val)
end
end;

# ╔═╡ d10036aa-eae6-43bb-893c-c11380f61d7f
begin
sliderobservables = [s.value for s in sl.sliders]
trx_settings = lift(sliderobservables...) do slvalues...
	trx.drv.fir[2] = slvalues[1]
    trx.drv.dcd = slvalues[2]/100
	trx.drv.rj_s = slvalues[3]*1e-12
	trx.drv.sj_freq = slvalues[4]*1e3
	trx.drv.sj_amp_ui = slvalues[5]
	trx.ch.noise_dbm_hz = slvalues[6]

	trx.drv.fir_norm = trx.drv.fir/sum(abs.(trx.drv.fir))
	trx.ch.noise_rms = sqrt(0.5/trx.param.dt*10^((trx.ch.noise_dbm_hz-30)/10)*trx.ch.noise_Z) #explicitly recalculate noise rms because it's not an observable
end
end;

# ╔═╡ 4c9f4620-8df6-4b3e-b836-9989ac5bb4cf
md"""
The most important UI element is the run button and this is similar to the one shown in the tutorial [here](https://www.youtube.com/watch?v=L-gyDvhjzGQ&t=612s&ab_channel=JuliaDynamics). Our animation step/frame is simply done by updating the ```eye_buffer``` observable, which in turn triggers the ```eye_heatmap``` generation, then Makie will update the plot. I wish I can say more about this, but it's really just this simple to make it work!
"""

# ╔═╡ 445f6bea-97b7-43ce-a51f-1a4ae15968db
begin
isrunning = Observable(false);
	
run_func1 = on(btn_run.clicks) do clicks
	isrunning[] = ~isrunning[]; 
	btn_run.label = isrunning[] ? "Stop" : "Run"
end;
	
run_func2 = on(btn_run.clicks) do clicks
	@async while isrunning[]
		isopen(f.scene) || break
		step_sim_blk(trx)
		eye_buffer[] = ch.Vo #this is the animation step
		yield()
		# sleep(0.001) #use sleep to slow down the frame rate if desired
	end
end
end;

# ╔═╡ 399da64a-2eff-4ae4-9e23-a498b6860a68
md"""
Now that we have reached the end, I will conclude not with too many words, but GIFs of what this widget can generate. **Enjoy!** 👏
"""

# ╔═╡ 3e9cc2f3-89b9-4a4d-887c-7111bdc051c5
md"""
### ISI

![demo_isi](https://circuit-artists.com/wp-content/uploads/2024/06/demo_isi.gif)
"""

# ╔═╡ e8308591-935e-4ec6-84c4-432a4a5ec8b1
md"""
### DCD
![demo_dcd](https://circuit-artists.com/wp-content/uploads/2024/06/demo_dcd.gif)
"""

# ╔═╡ 81fe8fe7-f7be-4b87-872c-8fa79ea4fe53
md"""
### Noise
![demo_noise](https://circuit-artists.com/wp-content/uploads/2024/06/demo_noise.gif)
"""

# ╔═╡ f9e43f6f-c725-4476-8fee-ac7c4c566296
md"""
### RJ
![demo_rj](https://circuit-artists.com/wp-content/uploads/2024/06/demo_rj.gif)
"""

# ╔═╡ bb8382df-3408-4f09-970b-e08aecf384c2
md"""
### SJ
![demo_sj](https://circuit-artists.com/wp-content/uploads/2024/06/demo_sj.gif)
"""

# ╔═╡ 63728684-c42e-40c3-aa58-ef23eb9701c7
md"""
## Helper and import section
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
ColorSchemes = "~3.25.0"
Colors = "~0.12.11"
DSP = "~0.7.9"
DataStructures = "~0.18.20"
Distributions = "~0.25.109"
FFTW = "~1.8.0"
GLMakie = "~0.10.0"
Interpolations = "~0.15.1"
MAT = "~0.10.6"
Makie = "~0.21.0"
Parameters = "~0.12.3"
StatsBase = "~0.34.3"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "c58efdc54e1225f014605a89aa56163c91a8fdcb"

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
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
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
git-tree-sha1 = "94997910aca72897524d2237c41eb852153b0f65"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.3"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ab3f7e1819dba9434a3a5126510c8fda3a4e7000"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.1+0"

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
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

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
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

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
git-tree-sha1 = "6e00379a24597be4ae1ee6b2d882e15392040132"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.5"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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
git-tree-sha1 = "a947ea21087caba0a798c5e494d0bb78e3a1a3a0"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.9"
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
# ╟─69e46020-2f48-11ef-0e3a-23d26a279b9d
# ╟─1ce6e5d5-7deb-49a2-945a-875619ccb86c
# ╟─bf3c796b-bcd8-493d-b036-de79991f6271
# ╟─f0a8f5d2-572c-4e0e-b863-b99848cc1a1a
# ╟─4decb1c3-7c59-4321-91eb-39d2395a093c
# ╟─e9edf872-886f-4937-8f00-6ed203630a06
# ╠═5b78e003-edbd-4791-9fd1-48d0b4523bd4
# ╟─b5d11536-5508-40b5-886c-8bbd6cbc73cc
# ╠═7973422d-fbe9-4049-8325-a620520c9be8
# ╟─d13dedd3-84c6-484e-b77f-98601d908474
# ╠═9fd804d6-1bbd-4697-9cef-c3f1113a1967
# ╠═303b9ac4-a7de-4d4d-89a7-aab55697cad0
# ╠═069c9012-f240-464e-95ed-9e1338b2d1eb
# ╟─64a9d1c8-ae20-4731-b339-503ba5bbee60
# ╠═688b1610-ff94-4f9f-9d9e-49df00f9977e
# ╠═55836ef7-0309-4bc1-8f43-9e5c28fe7100
# ╠═10bfdc04-14b7-423a-8191-2a3324b752db
# ╟─43a01dc0-cb9c-457f-ba87-2b8aaf81723a
# ╠═d3f53e27-5958-48d0-82ad-884500c60f15
# ╠═1422f9a9-024e-4715-b3cc-6f84961a8eac
# ╟─d0fcb7ce-2d6b-453d-b9c7-aed81595aba3
# ╟─37a3ef57-eace-48e7-b170-e8bebad356b7
# ╟─1985a856-9479-445b-a0c1-d2c7117e6e98
# ╟─6099bac2-6aa1-40b8-8f82-555d08ca1a84
# ╠═7987b2c5-cc99-4676-ac22-ca7a4b96258a
# ╠═986a900f-738e-43df-bf10-bcf292660ca0
# ╟─55367577-43fa-4e72-849a-21a584590808
# ╠═64c399e2-d200-49df-b633-ae74a0b23883
# ╟─b396f7fe-e1b0-4d46-8acb-a4fd0e1eda0b
# ╠═ea1bfbb5-471a-410f-89b3-72c031705f8f
# ╠═50026a47-7066-4e65-9eaf-40c53c59ca28
# ╠═4434df8c-fc34-4ffd-9c53-2b45df07f84a
# ╟─3923a75e-8c3a-4ab0-8508-1fdc17d5a853
# ╟─93b1f0ed-9de3-499e-99ad-87cb0f777c9d
# ╟─20de9f86-3fec-4944-9bc3-72adfd47dd60
# ╠═6d447a7f-6816-4164-8eb2-2e7e85f7ab00
# ╟─72af7e60-3398-44e3-ba85-28376c50dd3b
# ╠═0ee12efa-a36f-4638-88ab-2ed6094c6a81
# ╟─11c65bc7-9ca9-4dc8-9ad8-0f630fc7ed1e
# ╟─96d45520-b93e-4864-a950-8002afe214a5
# ╟─22e655f4-d8e9-4900-affa-546ab5026e81
# ╟─44516906-4475-42bb-a3db-c1de5a2ea04c
# ╟─8f93c1fb-c5bf-4ba7-8ae7-4923d10c57b4
# ╠═70d30065-2dde-446f-93f3-e86bc81d438f
# ╟─2f36a08e-3aeb-4dbd-a341-98e8f418ccb2
# ╠═f18af4fc-1bfc-45f3-88af-5feb27303155
# ╟─c1d734d9-5487-4a77-9bb8-30fdb905798d
# ╠═696c8ef8-863c-4b7e-84e8-3e579cb9255e
# ╟─47803c96-38b7-4554-a5e9-031966756b5f
# ╠═4a1aa2b3-8fdd-4bc9-bfed-8bd546a24bef
# ╠═ac5ad027-a8cf-44a7-86ed-c19a589019ce
# ╠═2cadbaef-eebb-4943-a08d-d70f62dc4c72
# ╠═c51fbd6a-b07c-4ddb-927a-fa19186db679
# ╟─b4434801-fabf-4152-a136-06c66a7abe9c
# ╠═78b50e9d-e650-459b-a291-e33235974534
# ╟─4dd5297b-0dd5-481c-b851-d9b4e009271c
# ╠═f269e7e6-42f0-4f49-a28a-252d1f572f8f
# ╟─7254e4f4-31f1-44ab-a9bf-c7d41828b7e2
# ╟─b4eb62d7-3942-4987-a106-ffe305c2a997
# ╠═3b9c6cf2-0667-480b-9d50-3a7d9952cc91
# ╠═8b44925f-d944-4498-a926-11ebb56db7d8
# ╠═5a04a0b6-1fdb-49cf-8f49-fb12db7f8597
# ╠═d10036aa-eae6-43bb-893c-c11380f61d7f
# ╟─4c9f4620-8df6-4b3e-b836-9989ac5bb4cf
# ╠═445f6bea-97b7-43ce-a51f-1a4ae15968db
# ╟─399da64a-2eff-4ae4-9e23-a498b6860a68
# ╟─3e9cc2f3-89b9-4a4d-887c-7111bdc051c5
# ╟─e8308591-935e-4ec6-84c4-432a4a5ec8b1
# ╟─81fe8fe7-f7be-4b87-872c-8fa79ea4fe53
# ╟─f9e43f6f-c725-4476-8fee-ac7c4c566296
# ╟─bb8382df-3408-4f09-970b-e08aecf384c2
# ╟─63728684-c42e-40c3-aa58-ef23eb9701c7
# ╠═8e12dc88-453c-4789-9079-3e3f7c789b9f
# ╠═33e5faec-0407-489e-a9e7-d90ce3ee9368
# ╠═81b6950c-9eb2-4a72-820c-e8a246776c56
# ╠═e032a9c9-d5dc-42eb-a321-0b70526e8e67
# ╠═3d901245-bae6-4664-bfbe-cc48c9ed0004
# ╠═656909c7-8a50-416c-a444-881d469c2a08
# ╠═e35022f4-ee62-46c1-9e22-1171429c2572
# ╠═8dc41adf-b61d-49da-a04d-54442a051a67
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
