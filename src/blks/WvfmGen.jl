module WvfmGen
using UnPack, Interpolations
using GLMakie, Makie, Printf
include("../util/Util_JLSD.jl"); using .Util_JLSD


export init_plot, w_plot_test, w_plot_test2, w_eye_gen_heatmap, update_eye
export w_gen_eye_simple, w_gen_eye_simple!
export w_newfig

function init_plot(wvfm)
    @unpack nrow, ncol, sizex, sizey = wvfm
    @unpack V11_x, V11_y, V21_x, V21_y, V31_x, V31_y = wvfm
    @unpack V12_x, V12_y, V22_x, V22_y = wvfm
    @unpack eye1, eslc_ref_ob = wvfm

    display(wvfm.screen, wvfm.fig)


    for c = 1:ncol
        for r = 1:nrow
            wvfm.axes[r,c] = Axis(wvfm.fig[r,c])
        end
    end

    ax11 = wvfm.axes[1,1]
    ax11.title = "Driver output"
    lines!(ax11, V11_x, V11_y)
    
    ax21 = wvfm.axes[2,1]
    ax21.title = "RX input"
    lines!(ax21, V21_x, V21_y)

    ax31 = wvfm.axes[3,1]
    ax31.title = "Sampled voltage"
    ylims!(ax31, -.4, .4)
    scatter!(ax31, V31_x, V31_y, alpha=0.2, markersize = 6)
    hlines!(ax31, eslc_ref_ob, color=:red, linestyle=:dash, linewidth=3)


    ax12 = wvfm.axes[1,2]
    ax12.title = "Error Slicer DAC code"
    lines!(ax12, V12_x, V12_y)

    ax22 = wvfm.axes[2,2]
    ax22.title = "PI code"
    lines!(ax22, V22_x, V22_y)

    ax_eye = wvfm.axes[3,2]
    ax_eye.title = "Eye"
    append!(eye1.buffer,zeros(eye1.buffer_size))
    eye1.heatmap_ob = lift(w_eye_gen_heatmap, eye1.heatmap_ob_trig, eye1)
    heatmap!(ax_eye, eye1.x_grid, eye1.y_grid, eye1.heatmap_ob, 
            colormap=eye1.colormap,
            inspector_label = (self, i, p) -> @sprintf("x = %.3f, y = %.3f", p[1], p[2]))
    hlines!(ax_eye, eslc_ref_ob, color=:red, linestyle=:dash, linewidth=3)

    DataInspector(wvfm.fig)

    display(wvfm.fig)
end



function w_plot_test(wvfm; cond = true)
    @unpack cur_blk, blk_size = wvfm.param

    if wvfm.en_plot & cond
            wvfm.V11_x.val = eachindex(wvfm.buffer11)
            wvfm.V11_y[] = wvfm.buffer11
    
            wvfm.V21_x.val = eachindex(wvfm.buffer21)
            wvfm.V21_y[] = wvfm.buffer21

            wvfm.V31_x.val = eachindex(wvfm.buffer31)
            wvfm.V31_y[] = wvfm.buffer31
            wvfm.eslc_ref_ob[] = wvfm.eslc_ref_ob.val

            wvfm.V12_x.val = eachindex(wvfm.buffer12)
            wvfm.V12_y[] = wvfm.buffer12

            wvfm.V22_x.val = eachindex(wvfm.buffer22)
            wvfm.V22_y[] = wvfm.buffer22

            wvfm.eye1.heatmap_ob_trig[] = true

            
            display(wvfm.screen, wvfm.fig)
            sleep(.0001)
            # yield()
    end

end

function w_gen_eye_simple(input,x_npts_ui, x_npts, y_range, y_npts; osr, x_ofst=0)
    heatmap = zeros(x_npts, y_npts)

    input_x = 0:1/osr:(lastindex(input)-1)/osr
    itp_resample = linear_interpolation(input_x, input)
    idx_itp = 0:1/x_npts_ui:input_x[end]
    input_itp = itp_resample.(idx_itp)

    for n = 1:x_npts
        heatmap[n,:] = u_hist(input_itp[n:x_npts:end], -y_range/2, y_range/2, y_npts)
    end
   
    return circshift(heatmap, (Int(x_ofst), 0))

end

function w_gen_eye_simple!(heatmap, input,x_npts_ui, x_npts, y_range, y_npts; osr, x_ofst=0, shadow=0.0)

    heatmap .*= shadow

    input_x = 0:1/osr:(lastindex(input)-1)/osr
    itp_resample = linear_interpolation(input_x, input)
    idx_itp = 0:1/x_npts_ui:input_x[end]
    input_itp = itp_resample.(idx_itp)

    for n = 1:x_npts
        k = Int(mod(n-1+x_ofst, x_npts)+1)
        heatmap[k,:] .+= (1.0-shadow) .* u_hist(input_itp[n:x_npts:end], -y_range/2, y_range/2, y_npts)
    end
   
    return nothing
end

function w_eye_gen_heatmap(heatmap_ob_trig, eye)
    if heatmap_ob_trig
        @unpack osr, tui = eye.param
        @unpack x_npts_ui, x_npts, y_npts, y_range = eye
        @unpack buffer, buffer_plt_len = eye
        @unpack clk_skews, clk_rj, noise_rms = eye

        x_npts_ui_val = x_npts_ui.val
        x_npts_val = x_npts.val
        y_npts_val = y_npts.val
        y_range_val = y_range.val

        if ~isempty(clk_skews) | (clk_rj > 0.0)
            #add 1UI in front and back
            buffer_x = 0:buffer_plt_len+2*osr-1
            itp_buffer = linear_interpolation(buffer_x, view(buffer, 1:buffer_plt_len+2*osr)) 
            
            n_resample = Int(buffer_plt_len/osr)
            buffer_resample = zeros(buffer_plt_len) 
            
            Φskew = zeros(n_resample)
            if ~isempty(clk_skews)
                Φskew .= kron(ones(Int(n_resample/length(clk_skews))), clk_skews/tui*osr)              
            end
            for n = 1:osr
                Φrj = clk_rj/tui*osr*randn(n_resample)
                Φnom = n:osr:n+n_resample*osr-1
                Φ =  osr .+ Φnom .+ Φskew .+ Φrj 

                buffer_resample[n:osr:end] .= itp_buffer.(Φ) 
            end  
        else
            buffer_resample = buffer[1:buffer_plt_len] #leaves 1UI in front and back
        end

        return w_gen_eye_simple(buffer_resample, x_npts_ui_val, x_npts_val, y_range_val, y_npts_val, osr=osr, x_ofst=eye.x_ofst)
    end
end

function update_eye(eye; x_nui, x_npts_ui, y_npts, x_ofst_ui)
    eye.x_ofst = round(-x_ofst_ui*x_npts_ui+x_nui*x_npts_ui/2)
    eye.x_nui.val = x_nui
    eye.x_npts_ui[] = x_npts_ui
    eye.y_npts[] = y_npts
    eye.heatmap_ob_trig[] = true
end



function w_newfig(nrow=1, ncol=1; kwargs...)
	sc = GLMakie.Screen();
	f = Figure(; kwargs...);
	display(sc,f)
	axes = Array{Axis}(undef,nrow,ncol)

	for r=1:nrow
		for c = 1:ncol
			axes[r,c] = Axis(f[r,c])
		end
	end
	
	return f, axes
end


end