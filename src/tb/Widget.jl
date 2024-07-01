using GLMakie, Makie
using UnPack, Random, Interpolations, ColorSchemes, Colors
include("../structs/TrxStruct.jl"); #using .TrxStruct
includet("../util/Util_JLSD.jl"); using .Util_JLSD
includet("../blks/BlkBIST.jl"); using .BlkBIST
includet("../blks/BlkTX.jl"); using .BlkTX
includet("../blks/BlkCH.jl"); using .BlkCH
includet("../blks/BlkRX.jl"); using .BlkRX
includet("../blks/WvfmGen.jl"); using .WvfmGen

function make_widget(trx)
    @unpack param, drv, ch, splr, clkgen, cdr = trx
    @unpack fig, screen, eye1 = trx.wvfm




    eye1.colormap[][1] = RGB(1.0,1.0,1.0)

    geye = GridLayout(fig[1:3,1:3])

    #plot the eye diagram, pass the observable directly in
    sl_xofst = Slider(geye[1,1:2], range = -eye1.x_npts[]/2:eye1.x_npts[]/2-1, startvalue = 0)
    ax_eye = Axis(geye[2, 1:2])
    heatmap!(ax_eye, eye1.x_grid, eye1.y_grid, eye1.heatmap_ob, 
            colormap=eye1.colormap)
    xlims!(eye1.x_grid[][1], eye1.x_grid[][end])

    Φo_hist = Observable(zeros(eye1.x_npts[]))
    ax_jitter = Axis(geye[3, 1:2], height = 60)
    band!(ax_jitter, eye1.x_grid, 0, Φo_hist)
    xlims!(eye1.x_grid[][1], eye1.x_grid[][end])
    hidespines!(ax_jitter, :t, :l, :r)
    hideydecorations!(ax_jitter)

    yk_hist = Observable(zeros(eye1.y_npts[]))
    ax_yk = Axis(geye[2, 3], width = 60)
    lower = Point2f.(0, eye1.y_grid[])
    upper = lift(yk_hist) do yk_h
        Point2f.(yk_h, eye1.y_grid[])
    end
    band!(ax_yk, lower, upper)
    ylims!(eye1.y_grid[][1], eye1.y_grid[][end])
    xlims!()
    hidespines!(ax_yk, :t, :b, :r)
    hidexdecorations!(ax_yk)

    display(screen, fig)
    # DataInspector(fig)

    gctrl = GridLayout(fig[1:3,4:5])
    #run button
    btn_run = Button(gctrl[1,1], label ="Run", tellwidth=false, tellheight=false, width = 140, height=60, fontsize=28);
    
    ch_lbl = ch.ch_en ? "Channel\nEnabled" : "Channel\nDisabled"
    btn_ch_en = Button(gctrl[1,2], label=ch_lbl, tellwidth = false, tellheight=false, width=100, height=60, fontsize=18);
 
    #slider for after image weight
    sl_config = SliderGrid(gctrl[2:3,:], 
        (label="Avg factor",  range=0:10, startvalue=eye1.shadow_weight), 
        (label="BG color",  range=0:0.01:1, startvalue=1.0), 
        (label="Frame wait",  range=0:100, format = "{:d} ms", startvalue=0.0), 
        tellheight=false, value_column_width=100)

    #sliders for TX/Channel parameters

    sl_tx = SliderGrid(gctrl[4:5,:],
        (label = "TX swing", range=0.4:0.02:1, startvalue = drv.swing),
        (label = "TX FIR h1", range=0:-0.05:-0.5, startvalue = drv.fir[2]),
        (label = "TX DCD", range=-10:10, format = "{:d} %", startvalue = drv.dcd*100),
        (label = "TX RJ", range=0:0.1:2, format = "{:.1f} ps", startvalue = drv.rj_s/1e-12),
        (label = "TX SJ freq", range=[0:10:90; 100:100:900; 1000:1000:1e4], format = "{:d} kHz", startvalue = drv.sj_freq/1e3),
        (label = "TX SJ amp", range = 0:0.05:1, format = "{:.2f} UI", startvalue = drv.sj_amp_ui),
        (label = "CH noise", range=-150:.5:-130, format = "{:.1f} dBm/Hz", startvalue = ch.noise_dbm_hz),
        tellheight=false, value_column_width=100)

    #enable channel button
    sl_rx = SliderGrid(gctrl[6:7,:],
        (label = "RX CDR Kp", range= 10:-1:0, format = "1/2^{:d}", startvalue = round(Int,-log2(cdr.kp))),
        (label = "RX CDR Ki", range= [32; 20:-1:8], format = "1/2^{:d}", startvalue = round(Int,-log2(cdr.ki))),
        (label = "RX IQ skew", range= -3:0.1:3, format = "{:.1f} ps", startvalue = clkgen.skews[2]/1e-12),
        (label = "RX RJ", range= 0:0.1:2, format = "{:.1f} ps", startvalue = clkgen.rj/1e-12),
        tellheight=false, value_column_width=100)



    #listeners   
    eye_ofst = lift(sl_xofst.value) do val
        eye1.x_ofst = val
    end

    sleep_sec = Observable(0.0)
    eye_config = lift([s.value for s in sl_config.sliders]...) do vals...
        eye1.shadow_weight = 1.0-.5^vals[1] #when weight =1, we purposely make it very close 1 becausue a perfect integrator can cause Float overflow
        eye1.colormap.val[1] = RGB(vals[2],vals[2],vals[2])
        eye1.colormap[] = eye1.colormap[]
        sleep_sec[] = 1e-3*vals[3]
    end

    ch_en_func = on(btn_ch_en.clicks) do clicks
        ch.ch_en = ~ch.ch_en;
        btn_ch_en.label = ch.ch_en ? "Channel\nEnabled" : "Channel\nDisabled"
    end

    txobservables = [s.value for s in sl_tx.sliders]
    tx_settings = lift(txobservables...) do slvals...
        drv.swing = slvals[1]
        drv.fir[2] = slvals[2]
        drv.dcd = slvals[3]/100
        drv.rj_s = slvals[4]*1e-12
        drv.sj_freq = slvals[5]*1e3
        drv.sj_amp_ui = slvals[6]
        ch.noise_dbm_hz = slvals[7]

        drv.fir_norm = drv.fir/sum(abs.(drv.fir))
        ch.noise_rms = sqrt(0.5/param.dt*10^((ch.noise_dbm_hz-30)/10)*ch.noise_Z) #explicitly recalculate noise rms because it's not an observable
    end

    rxobservables = [s.value for s in sl_rx.sliders]
    rx_settings = lift(rxobservables...) do slvals...
        cdr.kp = 1.0/2^slvals[1]
        cdr.ki = 1.0/2^slvals[2]
        clkgen.skews[[2,4]] .= slvals[3]*1e-12
        clkgen.rj = slvals[4]*1e-12

    end


    eye_buffer = Observable(splr.Vo)
    on(eye_buffer) do buffer
        w_gen_eye_simple!(eye1.heatmap_ob[], buffer, 
                        eye1.x_npts_ui[], eye1.x_npts[], eye1.y_range[], eye1.y_npts[]; 
                        osr = trx.param.osr, x_ofst=eye1.x_ofst, shadow = eye1.shadow_weight)
        eye1.heatmap_ob[] = eye1.heatmap_ob[]
        # eye1.heatmap_ob[] = randn(eye1.x_npts[],eye1.y_npts[])
    end


    isrunning = Observable(false)

    run_func1 = on(btn_run.clicks) do clicks
        isrunning[] = ~isrunning[]; 
        btn_run.label = isrunning[] ? "Stop" : "Run"
    end
        
    run_func2 = on(btn_run.clicks) do clicks
        @async while isrunning[]
                isopen(wvfm.fig.scene) || break
                step_sim_blk(trx)
                eye_buffer[] = splr.Vo #this is the animation step
                Φo_wrap = mod.(clkgen.Φo./param.osr .+ (eye1.x_ofst/eye1.x_npts_ui[]), eye1.x_nui[]).- (eye1.x_nui[]/2)
                Φo_hist.val .*= eye1.shadow_weight
                Φo_hist.val .+= (1-eye1.shadow_weight) .* u_hist(Φo_wrap, -eye1.x_nui[]/2, eye1.x_nui[]/2, eye1.x_npts[])
                Φo_hist[] = Φo_hist[] 

                yk_hist.val .*= eye1.shadow_weight
                yk_hist.val .+= (1-eye1.shadow_weight) .* u_hist(splr.So, -eye1.y_range[]/2, eye1.y_range[]/2, eye1.y_npts[])
                yk_hist[] = yk_hist[]
                ylims!(ax_jitter)
                xlims!(ax_yk)
                sleep(sleep_sec[]) #use sleep to slow down the frame rate if desired
        end
    end
end


function step_sim_blk(trx)
    @unpack param, bist, drv, ch = trx 

    pam_gen_top!(bist)

    dac_drv_top!(drv, bist.So)
    
    ch_top!(ch, drv.Vo)

    sample_itp_top!(splr, ch.Vo)

    run_blk_iter(trx, 0, param.nsubblk, step_sim_subblk)


end

function step_sim_subblk(trx, blk_idx)
    @unpack clkgen, splr = trx 
    @unpack dslc, eslc, cdr, adpt = trx

    trx.param.cur_subblk = blk_idx

    clkgen_pi_itp_top!(clkgen, pi_code=cdr.pi_code)
        
    sample_phi_top!(splr, clkgen.Φo_subblk)

    slicers_top!(dslc, splr.So_subblk, ref_code=[[128],[128],[128],[128]])
    slicers_top!(eslc, splr.So_subblk, ref_code=adpt.eslc_ref_vec)

    cdr_top!(cdr, dslc.So, eslc.So)

    adpt_top!(adpt, dslc.So, eslc.So)

end
