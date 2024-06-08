using Revise, BenchmarkTools
using GLMakie, Makie
using UnPack, Random, Interpolations
include("../structs/TrxStruct.jl")
includet("../util/Util_JLSD.jl"); using .Util_JLSD
includet("../blks/BlkBIST.jl"); using .BlkBIST
includet("../blks/BlkTX.jl"); using .BlkTX
includet("../blks/BlkCH.jl"); using .BlkCH
includet("../blks/BlkRX.jl"); using .BlkRX
includet("../blks/WvfmGen.jl"); using .WvfmGen



function init_trx()
    #global params for simulation
    params = TrxStruct.Params(  
                data_rate = 56e9,
                pam = 2,
                osr = 32,
                blk_size = 2^10,
                subblk_size = 16, 
                nsym_total = Int(2e6))
    Random.seed!(params.rand_seed)

    #bist params
    bist = TrxStruct.Bist(  
                params = params,
                polynomial = [28,31])


    #TX driver params
    drv = TrxStruct.Drv(
                params = params,
                ir = u_gen_ir_rc(params.dt, params.fbaud, 20*params.tui),
                fir = [-0.0, 1., -0.2],
                swing = 0.8,
                jitter_en = true,
                dcd = 0.03,
                rj_s = 300e-15,
                sj_amp_ui = 0.2,
                sj_freq = 1e5)

    #AWGN ch params
    ch = TrxStruct.Ch(
                params = params,
                ir_ch = u_fr_to_imp("./channel_data/TF_data/channel_4inch.mat", 
                        params.tui, params.osr, npre = 20, npost= 79),
                ir_pad = u_gen_ir_rc(params.dt, params.fbaud, 20*params.tui),
                noise_en = true,
                noise_rms = 2e-3)

    #clkgen params
    clkgen = TrxStruct.Clkgen(
                params = params,
                nphases = 4,
                rj = .3e-12,
                skews = [0e-12, 0e-12, 0e-12, 0e-12])

    #sampler params
    splr = TrxStruct.Splr(
                params = params,
                ir = u_gen_ir_rc(params.dt, params.fbaud, 20*params.tui))

    #slicer params
    dslc = TrxStruct.Slicers(
                params = params,
                N_per_phi = ones(UInt8, clkgen.nphases),
                noise_rms = 1.5e-3,
                ofst_std = 7e-3)

    eslc = TrxStruct.Slicers(
                params = params,
                N_per_phi = UInt8.([1,0,0,0]),
                noise_rms = 1.5e-3,
                ofst_std = 7e-3)

    cdr = TrxStruct.Cdr(
                params = params,
                Neslc_per_phi = eslc.N_per_phi,
                kp = 1/2^6,
                ki = 1/2^14,
                pi_res = clkgen.pi_res)
    
    adpt = TrxStruct.Adpt(
                params = params,
                Neslc_per_phi = eslc.N_per_phi,
                mu_eslc = 1/64)

    #waveform plotting params
    wvfm = TrxStruct.Wvfm(
                params = params,
                en_plot = true,
                nrow = 3,
                ncol = 2)
    wvfm.eye1.clk_skews = clkgen.skews
    wvfm.eye1.clk_rj = clkgen.rj
    wvfm.eye1.noise_rms = dslc.noise_rms

    init_plot(wvfm)
    
    println("init done")

    return (;params, bist, drv, ch, clkgen, splr, dslc, eslc, cdr, adpt, wvfm)
end


function run_sim_blk(trx)
    @unpack params, bist, drv, ch, clkgen, splr = trx 
    @unpack dslc, eslc, cdr, adpt, wvfm = trx

    pam_gen_top!(bist)

    dac_drv_top!(drv, bist.So)
    
    ch_top!(ch, drv.Vo)

    sample_itp_top!(splr, ch.Vo)

    for m = 1:params.nsubblk
        params.cur_subblk = m
        
        clkgen_pi_itp_top!(clkgen, pi_code=cdr.pi_code)
        
        sample_phi_top!(splr, clkgen.Φo)

        slicers_top!(dslc, splr.So_subblk, ref_code=[[0],[0],[0],[0]])
        slicers_top!(eslc, splr.So_subblk, ref_code=adpt.eslc_ref_vec)

        cdr_top!(cdr, dslc.So, eslc.So)

        adpt_top!(adpt, dslc.So, eslc.So)

        append!(bist.Si, [sum(dvec) for dvec in dslc.So])

        push!(wvfm.buffer12, adpt.eslc_ref_code)
        push!(wvfm.buffer22, cdr.pi_code)
    end

    ber_checker_top!(bist)

    #record waveform here
    append!(wvfm.buffer11, drv.Vo)
    append!(wvfm.buffer21, ch.Vo)
    append!(wvfm.buffer31, splr.So)
    append!(wvfm.eye1.buffer, splr.Vo)
    wvfm.eslc_ref_ob.val = adpt.eslc_ref_code*eslc.dac_lsb+eslc.dac_min
    wvfm.eye1.x_ofst = round(-(cdr.pi_code/clkgen.pi_codes_per_ui)*wvfm.eye1.x_npts_ui[]+wvfm.eye1.x_npts[]/2)
    w_plot_test(wvfm, cond=(params.cur_blk%wvfm.plot_every_nblk==0))
end


function run_sim(trx, sim_blk)
    for n = 1:trx.params.nblk
        trx.params.cur_blk = n
        sim_blk(trx)
    end
end

    