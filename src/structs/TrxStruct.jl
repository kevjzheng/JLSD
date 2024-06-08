module TrxStruct
using Parameters, DataStructures, Interpolations
using GLMakie, Makie


export Params, Bist, Drv, Ch, Clkgen, Splr, Slicers, Cdr, Adpt, Eye, Wvfm

@kwdef mutable struct Params
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
    const pr_sinc::Vector{Float64} = sinc.(-16:1/osr:16-1/osr)

    cur_blk = 0
    cur_subblk = 0

end

@kwdef mutable struct Bist
    const params::Params
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

    So_bits::Vector = zeros(Bool, params.bits_per_sym*params.blk_size)
    So::Vector = zeros(params.blk_size)
    Si = CircularBuffer{UInt8}(params.blk_size)
    Si_bits::Vector = zeros(Bool, params.bits_per_sym*params.blk_size)
end


@kwdef mutable struct Drv
    const params::Params
    ir::Vector{Float64}
    rlm_en = false
    rlm = 1.0

    swing = 0.7
    fir = [1,0]
    fir_norm = fir/sum(abs.(fir))
    quantize = false
    dac_res = 7 

    jitter_en = false
    dcd = 0.0
    rj_s = 0.0
    sj_amp_ui = 0.0
    sj_freq = 0.0
    last_rj_osr = 0.0
    last_sj_phi = 0.0

    Sfir::Vector = zeros(params.blk_size)
    Sfir_mem::Vector = zeros(length(fir)-1)

    prev_nui = 4
    V_prev_nui::Vector = zeros(prev_nui*params.osr)
    Vext::Vector = zeros(prev_nui*params.osr+params.blk_size_osr)
    tt_Vext = -prev_nui/2*params.osr:length(Vext)-prev_nui/2*params.osr-1

    Vo::Vector = zeros(params.blk_size_osr) 
    Vo_mem::Vector = zeros(lastindex(ir)-1)

end

@kwdef mutable struct Ch
    const params::Params
    noise_en = true
    noise_rms::Float64 = 0.0

    ir_ch::Vector{Float64}
    Vch::Vector = zeros(params.blk_size_osr) 
    Vch_mem::Vector = zeros(lastindex(ir_ch)-1) 
    
    No::Vector = zeros(params.blk_size_osr)
    No_mem::Vector = zeros(lastindex(params.pr_sinc)-1)

    ir_pad::Vector{Float64}
    Vo::Vector = zeros(params.blk_size_osr) 
    Vo_mem::Vector = zeros(lastindex(ir_pad)-1)
end

# @with_kw mutable struct Ctle
#     params::Params

#     bypass::Bool = false
#     nonlin_en::Bool = true
#     code
#     min_hf_code
#     max_hf_code

#     hf_ir_table
#     hf_dc_table

#     ir::AbstractArray
#     dc_tf::AbstractArray
#     ofst = 0.0
#     noise_rms = 0.0

#     oscal_en::Bool = false
#     oscal_dac_range = 60e-3
#     oscal_dac_res = 6
#     oscal_dac_lsb = oscal_dac_range/2^oscal_dac_res
#     oscal_code = 2^(oscal_dac_res-1)

#     No = zeros(params.blk_size_osr)
#     No_mem = zeros(1)

#     Vo
# end

@kwdef mutable struct Clkgen
    const params::Params

    nphases::Int8
    rj = 0
    skews = zeros(nphases)

    pi_res = 8
    pi_max_code = 2^pi_res-1
    pi_ui_cover = 4
    pi_codes_per_ui = 2^pi_res/pi_ui_cover
    pi_nonlin_lut = zeros(2^pi_res) #introduce INL here
    pi_code_prev = 0
    pi_wrap_ui = 0
    pi_wrap_ui_Δcode = pi_max_code-10

    Φo::Vector = zeros(params.subblk_size)

end

@kwdef mutable struct Splr
    const params::Params

    
    ir::Vector{Float64}

    Vo::Vector = zeros(params.blk_size_osr)
    Vo_mem::Vector = zeros(length(ir)-1)

    prev_nui = 16
    V_prev_nui::Vector = zeros(prev_nui*params.osr)
    Vext::Vector = zeros(prev_nui*params.osr+params.blk_size_osr)
    tt_Vext = -prev_nui/2*params.osr:length(Vext)-prev_nui/2*params.osr-1

    So::Vector = CircularBuffer{Float64}(params.blk_size)
    So_subblk::Vector = zeros(params.subblk_size)
    itp_Vext = nothing
end

@kwdef mutable struct Slicers
    const params::Params

    N_per_phi::Vector
    nphases = length(N_per_phi)
    noise_rms = 0.0
    ofst_std = 0.0
    
    ofsts = [ofst_std*randn(n) for n in N_per_phi]

    dac_res = 8
    dac_min = 0
    dac_max = 0.5
    dac_lsb = (dac_max-dac_min)/2^dac_res


    So = [zeros(Bool, Int(N_per_phi[n%nphases+1])) for n in 0:params.subblk_size-1]

end

@kwdef mutable struct Cdr
    const params::Params
    
    Neslc_per_phi::Vector
    nphases = length(Neslc_per_phi)

    Sd_prev = 0
    eslc_nvec = kron(ones(Int(params.subblk_size/nphases)),Neslc_per_phi)
    filt_patterns = [[0, 1, 1], [1, 1, 0]]

    kp = 1/32
    ki = 0/1024
    pd_gain = 1.0

    ki_accum = 0.0
    pd_accum = 128.0

    pi_res = 8
    pi_code = Int(floor(pd_accum))
end

@kwdef mutable struct Adpt 
    const params::Params

    Neslc_per_phi::Vector
    nphases = length(Neslc_per_phi)

    Sd_prev = 0
    eslc_nvec = kron(ones(Int(params.subblk_size/nphases)),Neslc_per_phi)
    eslc_filt_patterns = [[0, 1, 1], [1, 1, 0]]
    eslc_ref_accum = 128.0
    mu_eslc = 1/1024
    eslc_ref_max = 255
    eslc_ref_code = floor(eslc_ref_accum)
    eslc_ref_vec = [eslc_ref_code*ones(Int,n) for n in Neslc_per_phi]

    
end

@kwdef mutable struct Eye
    const params::Params
    x_npts_ui = Observable(Int(params.osr))
    x_nui = Observable(1)
    x_npts = @lift($x_npts_ui*$x_nui)
    x_grid = @lift(-$x_nui/2 : 1/$x_npts_ui: $x_nui/2-1/$x_npts_ui)
    x_ofst = 0
    y_npts = Observable(64)
    y_range = Observable(0.8)
    y_grid_edge = @lift(-$y_range/2: $y_range/$y_npts: $y_range/2)
    y_grid = @lift(($y_grid_edge[1:end-1] .+ $y_grid_edge[2:end])/2)
    buffer_size::Int = 2^15*params.osr
    buffer = CircularBuffer{Float64}(buffer_size)
    buffer_plt_len::Int = 2^14*params.osr
    heatmap_ob_trig = Observable{Bool}(true)
    heatmap_ob = Observable{Matrix{Float64}}(zeros(x_npts.val, y_npts.val))
    colormap = Observable(:turbo)

    clk_skews::Vector{Float64} = Float64[]
    clk_rj = 0.0
    noise_rms = 0.0
end

@kwdef mutable struct Wvfm
    const params::Params

    en_plot = true
    plot_every_nblk = Int(round(5e5/params.blk_size))

    sizex = 800
    sizey = 840 
    fig = Figure(backgroundcolor = :white, size = (sizex, sizey))
    nrow::Int8 = 3
    ncol::Int8 = 2
    axes::Array = Array{Axis}(undef,nrow,ncol)

    buffer11 = CircularBuffer{Float64}(8192)
    # buffer12 = CircularBuffer{Float64}(2^16)
    buffer12 = Float64[]

    buffer21 = CircularBuffer{Float64}(8192)
    # buffer22 = CircularBuffer{Float64}(2^16)
    buffer22 = Float64[]

    buffer31 = CircularBuffer{Float64}(8192)
    
    V11_x = Observable(zeros(1))
    V11_y = Observable(zeros(1))

    V21_x = Observable(zeros(1))
    V21_y = Observable(zeros(1))

    V31_x = Observable(zeros(1))
    V31_y = Observable(zeros(1))

    V12_x = Observable(zeros(1))
    V12_y = Observable(zeros(1))

    V22_x = Observable(zeros(1))
    V22_y = Observable(zeros(1))

    eye1 = Eye(params=params)
    
    eslc_ref_ob = Observable(0.0)
end


end