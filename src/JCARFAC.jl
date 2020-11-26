"""
    JCARFAC

Wrapper material for calling a simple C library implementation of Richard Lyon's
CARFAC-AGC (Cascade of Asymmetric Resonators with Fast-Acting Compression and
Automatic Gain Control) model of the auditory periphery.

SAI (Stabilized Auditory Image) functionality forthcoming.  
"""
module JCARFAC

using Match

export calcnap, defaultparams, BMParams, IHCParams, OHCParams

################################################################################
# Data Types
################################################################################

struct BMParams
    num_sections::Cint
    x_lo::Cdouble
    x_hi::Cdouble
    damping::Cdouble
end

struct IHCParams
    hpf_cutoff::Cdouble
    tau_in::Cdouble
    tau_out::Cdouble
    tau_ihc::Cdouble
end

struct OHCParams
    scale::Cdouble
    offset::Cdouble
    b::Ptr{Cdouble} 
end

struct SAIParams
    trig_win_t::Cdouble
    adv_t::Cdouble
    num_trig_win::Cint
    num_sections::Cint
    num_samples::Cint
end

struct SAI
    num_frames::Cint
    num_sections::Cint
    frame_len_n::Cint
    images::Ptr{Cdouble} 
    delay_t::Ptr{Cdouble} 
    times::Ptr{Cdouble} 
end

struct CARFACAGCState
    num_sections::Cint
    block_size::Cint

    f::Ptr{Cdouble}
    a0::Ptr{Cdouble}
    c0::Ptr{Cdouble}
    r::Ptr{Cdouble}
    r1::Ptr{Cdouble}
    h::Ptr{Cdouble}
    g::Ptr{Cdouble}

    q::Cdouble
    c_in::Cdouble
    c_out::Cdouble
    c_ihc::Cdouble

    scale::Cdouble
    offset::Cdouble
    b::Ptr{Cdouble}
    d_rz::Ptr{Cdouble}

    c_agc::Ptr{Cdouble}
    sa::Ptr{Cdouble}
    sb::Ptr{Cdouble}
    sc::Ptr{Cdouble}

    bm::Ptr{Cdouble}
    bm_hpf::Ptr{Cdouble}
    ihc_out::Ptr{Cdouble}
    ihc_state::Ptr{Cdouble}
    w0::Ptr{Cdouble}
    w1::Ptr{Cdouble}
    w1_old::Ptr{Cdouble}
    trans::Ptr{Cdouble}
    acc8::Ptr{Cdouble}
    acc16::Ptr{Cdouble}
    acc32::Ptr{Cdouble}
    acc64::Ptr{Cdouble}
    agc::Ptr{Cdouble}
    agc0::Ptr{Cdouble}
    agc1::Ptr{Cdouble}
    agc2::Ptr{Cdouble}
    agc3::Ptr{Cdouble}

    ihc_new::Ptr{Cdouble}
    z::Ptr{Cdouble}
    v_mem::Ptr{Cdouble}
    v_ohc::Ptr{Cdouble}
    sqrd::Ptr{Cdouble}
    nlf::Ptr{Cdouble}
    prev::Cdouble
    w0_new::Cdouble
end


################################################################################
# Utility Functions (for internal use)
################################################################################
"""
    init(bmp::Ref{BMParams}, ihcp::Ref{IHCParams}, ohcp::Ref{OHCParams},

Initialize the CARFAC state structure.
"""
function init(bmp::Ref{BMParams}, ihcp::Ref{IHCParams}, ohcp::Ref{OHCParams},
              blocksize::Int, fs::Float64)
    blocksize = Int32(blocksize)
    carfacstateptr = ccall(
                        (:carfacagc_init, "../ccarfac/src/libcarfac.so"),
                        Ptr{CARFACAGCState},
                        (Ref{BMParams}, Ref{IHCParams}, Ref{OHCParams}, Int32, Float64),
                        bmp, ihcp, ohcp, blocksize, fs)
    if carfacstateptr == C_NULL
        error("CARFAC initialization process returned null.")
    end
    return carfacstateptr
end

init(bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams, blocksize::Int, fs::Float64) = 
    init(Ref(bmp), Ref(ihcp), Ref(ohcp), blocksize, fs)


"""
    processblock(carfacstate::Ptr{CARFACAGCState}, signal::Array{Float64,1})

Process a block of input.
"""
function processblock(carfacstate::Ptr{CARFACAGCState}, signal::Array{Float64,1})
    # @assert carfacstate.num_sections == length(signal) "Block size doesn't match signal length"
    signal_ref = Base.unsafe_convert(Ptr{Cdouble}, Base.cconvert(Ptr{Cdouble}, signal))
    success = ccall(
                    (:carfacagc_process_block, "../ccarfac/src/libcarfac.so"),
                    Cint,
                    (Ptr{CARFACAGCState},Ptr{Cdouble}),
                    carfacstate, signal_ref)
    if success != 0
        error("CARFAC failed to process the block.")
    end
    return carfacstate
end


"""
    copynap(carfacstate::Ptr{CARFACAGCState})

Pull the output of the inner hair cells from the CARFAC state data so it is
available in Julia.

(NOTE: transferring memory ownership through own=false in unsafe_unwrap
seemed to cause problems when using carfacagc_free_except_signal...

Here we're just copying the data over so the whole structure can be freed at once.)
"""
function copynap(carfacstate::Ptr{CARFACAGCState})
    carfacstate_local = Base.unsafe_load(carfacstate)
    numsections = carfacstate_local.num_sections
    numpoints = carfacstate_local.block_size
    fcs = unsafe_wrap(Array, carfacstate_local.f, numsections, own=false) |> copy
    nap = unsafe_wrap(Array, carfacstate_local.ihc_out, numsections*numpoints, own=false) |> copy
    fcs = reverse(fcs)
    nap = nap |> x -> reshape(x, (numpoints, numsections)) |> transpose |> x -> reverse(x, dims=1)
    # success = ccall(
    #               (:carfacagc_free_except_signal,
    #               "/home/dahlbom/research/ccarfac/libcarfac.so"),
    #               Cint,
    #               (Ptr{CARFACAGCState},),
    #               carfacstate)
    # if success != 0
    #   error("CARFAC failed to free state data.")
    # end
    return nap, fcs 
end


"""
    free(carfacstate::Ptr{CARFACAGCState})

Free the CARFAC state data.
"""
function free(carfacstate::Ptr{CARFACAGCState})
    success = ccall(
          (:carfacagc_free,"../ccarfac/src/libcarfac.so"),
          Cint,
          (Ptr{CARFACAGCState},),
          carfacstate)
    if success != 0
        error("Failed to free CARFAC state.")
    end
    return nothing
end

################################################################################
# Public Interface
################################################################################
"""
    defaultparams(numsections=72)

Generate default parameters for the CARFAC model. Returns parameter structures
for the basilar membrane (bmp), inner hair cells (ihcp), and outer hair cells
(ohcp).
"""
function defaultparams(numsections=72)
    bvals = ones(Float64,numsections) 
    bvals_ref = Base.unsafe_convert(Ptr{Cdouble}, Base.cconvert(Ptr{Cdouble}, bvals)) 
    ohcp  = OHCParams(0.1,0.04,bvals_ref)
    bmp   = BMParams(numsections,0.1,0.9,0.2)
    ihcp  = IHCParams(20.0,10.0e-3,0.5e-3,80.0e-6)
    return bmp, ihcp, ohcp
end


"""
    calcnap(bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams,
            signal::Array{Float64,1}, fs::T) where T <: Number

Calculate the neural activity pattern (NAP) generated by a signal. This
is the main interface to the module. Returns a 2D array containing
the NAP as well as a set of center frequencies corresponding to each channel.
"""
function calcnap(bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams,
                signal::Array{Float64,1}, fs::T) where T <: Number
    state = init(bmp, ihcp, ohcp, length(signal), Float64(fs))
    processblock(state, signal)
    nap, fcs = copynap(state)
    free(state)
    return nap, fcs
end

calcnap(signal::Array{Float64,1}, fs::T) where T <: Number = calcnap(defaultparams()..., signal, fs)


end
