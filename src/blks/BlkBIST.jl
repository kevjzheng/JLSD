module BlkBIST
using UnPack, Random, Distributions

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

function bist_prbs_gen!(seq; poly, inv, Nsym, seed)
    for n = 1:Nsym
        seq[n] = inv
        for p in poly
            seq[n] ⊻= seed[p]
        end
        seed .= [seq[n]; seed[1:end-1]]
    end
    return nothing
end


function pam_gen_top!(bist)
    @unpack pam, bits_per_sym, blk_size = bist.param
    @unpack polynomial, inv, gen_seed = bist
    @unpack gen_gray_map, gen_en_precode, gen_precode_prev_sym, So_bits, So = bist

    bist_prbs_gen!(So_bits, poly=polynomial, inv=inv,
                    Nsym=bits_per_sym*blk_size, seed=gen_seed)


    fill!(So, zero(Float64))
    for n = 1:bits_per_sym
        @. So = So + 2^(bits_per_sym-n)*So_bits[n:bits_per_sym:end]
    end

    #gray encoding
    if ~isempty(gen_gray_map)
        for n in 1:blk_size
            So[n] = gen_gray_map[So[n] + 1]
        end
    end

    if gen_en_precode
        for n = 1:blk_size
            So[n] = mod(So[n]-gen_precode_prev_sym , pam)
            gen_precode_prev_sym = So[n]
        end
        bist.gen_precode_prev_sym = gen_precode_prev_sym #need to write back
    end
    @. So = 2/(pam-1)*So - 1

    return nothing
end

function ber_checker_top!(bist)
    @unpack cur_blk, pam, bits_per_sym = bist.param
    @unpack gen_gray_map, chk_precode_prev_sym, chk_start_blk, Si, Si_bits = bist 


    if cur_blk >= chk_start_blk #make start blk a parameter later
        if bist.gen_en_precode
            bist.chk_precode_prev_sym = Si[end]
            Si .= mod.([chk_precode_prev_sym; Si[1:end-1]].+ Si , pam)
        end
    
        if ~isempty(gen_gray_map)
            for n in 1:blk_size
                Si[n] = gen_gray_map[Si[n] + 1]
            end
        end
    
        Si_bits .= vec(stack(int2bits.(Si, bits_per_sym)))

        ber_check_prbs!(bist)
    end
end

function ber_check_prbs!(bist)
    @unpack polynomial, inv, chk_seed, ref_bits, Si_bits = bist
    nbits_rcvd = lastindex(Si_bits)

    # err_loc = rand(Uniform(0,1.0), nbits_rcvd).< 1e-4;
    # Si_bits .= Si_bits .⊻ err_loc

    if bist.chk_lock_status
        bist_prbs_gen!(ref_bits, poly=polynomial, inv=inv,
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
                ~, chk_seed = bist_prbs_gen(poly=polynomial, inv=inv,
                                            Nsym=nbits_rcvd-n, seed=chk_seed)
                #run prbs towards the end of the block to get the right seed
                break
            end
        end
    end
end

function int2bits(num, nbit)
    return [Bool((num>>k)%2) for k in nbit-1:-1:0]
end



end