using Revise, BenchmarkTools

include("./tb/Widget.jl")


GLMakie.closeall()

#global param for simulation
param = TrxStruct.Param(  
    data_rate = 56e9,
    pam = 2,
    osr = 20,
    blk_size = 512,
    subblk_size = 32, 
    nsym_total = Int(1e6))
Random.seed!(param.rand_seed)

#bist param
bist = TrxStruct.Bist(  
    param = param,
    polynomial = TrxStruct.PRBS31)


#TX driver param
drv = TrxStruct.Drv(
    param = param,
    ir = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui),
    fir = [1., -0.25],
    swing = 0.8,
    jitter_en = true,
    dcd = 0.03,
    rj_s = 300e-15,
    sj_amp_ui = 0.0,
    sj_freq = 2e5)

#AWGN ch param
ch = TrxStruct.Ch(
    param = param,
    ir_ch = u_fr_to_imp("./channel_data/TF_data/channel_4inch.mat", 
            param.tui, param.osr, npre = 20, npost= 79),
    ir_pad = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui),
    noise_en = true,
    noise_dbm_hz = -150 )

#clkgen param
clkgen = TrxStruct.Clkgen(
    param = param,
    nphases = 4,
    rj = .3e-12,
    skews = [0e-12, 0e-12, 0e-12, 0e-12])

#sampler param
splr = TrxStruct.Splr(
    param = param,
    ir = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui))

#slicer param
dslc = TrxStruct.Slicers(
    param = param,
    N_per_phi = ones(UInt8, clkgen.nphases),
    noise_rms = 1.5e-3,
    ofst_std = 7e-3,
    dac_min = -0.1,
    dac_max = 0.1)

eslc = TrxStruct.Slicers(
    param = param,
    N_per_phi = UInt8.([1,0,0,0]),
    noise_rms = 1.5e-3,
    ofst_std = 7e-3,
    dac_min = 0,
    dac_max = 0.5)

cdr = TrxStruct.Cdr(
    param = param,
    Neslc_per_phi = eslc.N_per_phi,
    kp = 1/2^5,
    ki = 1/2^16,
    pi_res = clkgen.pi_res)

adpt = TrxStruct.Adpt(
    param = param,
    Neslc_per_phi = eslc.N_per_phi,
    mu_eslc = 1/64)

#waveform plotting param
wvfm = TrxStruct.Wvfm(
    param = param,
    sizex = 1200,
    sizey = 800,
    en_plot = true)

wvfm.eye1.x_nui[] = 2
wvfm.eye1.x_npts_ui[] = 256
wvfm.eye1.y_npts[] = 256

println("init done")

trx =  (;param, bist, drv, ch, clkgen, splr, dslc, eslc, cdr, adpt, wvfm)


@time make_widget(trx);


