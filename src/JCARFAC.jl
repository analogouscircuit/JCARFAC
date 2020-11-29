"""
    JCARFAC

Wrapper material for calling a simple C library implementation of Richard Lyon's
CARFAC-AGC (Cascade of Asymmetric Resonators with Fast-Acting Compression and
Automatic Gain Control) model of the auditory periphery.

SAI (Stabilized Auditory Image) functionality forthcoming.  
"""
module JCARFAC

using Match

export calcnap, calcsai, defaultparams, BMParams, IHCParams, OHCParams, JSAI

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

struct JSAI
    numframes::Int64
    numsections::Int64
    framelen::Int64
    images::Array{Float64,3}
    delays::Array{Float64,1}
    times::Array{Float64,1}
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
    initcarfac(bmp::Ref{BMParams}, ihcp::Ref{IHCParams}, ohcp::Ref{OHCParams},

Initialize the CARFAC state structure.
"""
function initcarfac(bmp::Ref{BMParams}, ihcp::Ref{IHCParams}, ohcp::Ref{OHCParams},
              blocksize::Int, fs::Float64)
    blocksize = Int32(blocksize)
    carfacstateptr = ccall(
                        (:carfacagc_init, "/home/dahlbom/research/carfac/ccarfac/src/libcarfac.so"),
                        Ptr{CARFACAGCState},
                        (Ref{BMParams}, Ref{IHCParams}, Ref{OHCParams},
                         Int32, Float64),
                        bmp, ihcp, ohcp, blocksize, fs)
    if carfacstateptr == C_NULL
        error("CARFAC initialization process returned null.")
    end
    return carfacstateptr
end

initcarfac(bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams, blocksize::Int, fs::Float64) = 
    initcarfac(Ref(bmp), Ref(ihcp), Ref(ohcp), blocksize, fs)


"""
    processblock(carfacstate::Ptr{CARFACAGCState}, signal::Array{Float64,1})

Process a block of input.
"""
function processblock(carfacstate::Ptr{CARFACAGCState}, signal::Array{Float64,1})
    # @assert carfacstate.num_sections == length(signal) "Block size doesn't match signal length"
    signal_ref = Base.unsafe_convert(Ptr{Cdouble}, Base.cconvert(Ptr{Cdouble}, signal))
    success = ccall(
                    (:carfacagc_process_block, "/home/dahlbom/research/carfac/ccarfac/src/libcarfac.so"),
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
"""
function copynap(carfacstate::Ptr{CARFACAGCState}; julialayout=true)
    carfacstate_local = Base.unsafe_load(carfacstate)
    numsections = carfacstate_local.num_sections
    numpoints = carfacstate_local.block_size
    fcs = unsafe_wrap(Array, carfacstate_local.f,
                      numsections, own=false) |> copy
    nap = unsafe_wrap(Array, carfacstate_local.ihc_out,
                      numsections*numpoints, own=false) |> copy
    fcs = reverse(fcs)
    if julialayout
        nap = nap |>
              x -> reshape(x, (numpoints, numsections)) |>
              transpose |>
              x -> reverse(x, dims=1)
    end
    return nap, fcs 
end


"""
    copysai(sai::Ptr{SAI})

Copy SAI data into Julia structure.  
"""
function copysai(sai::Ptr{SAI})
    # Copy data locally
    sai_local = Base.unsafe_load(sai)
    numframes = Int64(sai_local.num_frames)
    numsections = Int64(sai_local.num_sections)
    framelen = Int64(sai_local.frame_len_n)
    delays = unsafe_wrap(Array, sai_local.delay_t, framelen, own=false) |> copy
    times = unsafe_wrap(Array, sai_local.times, numframes, own=false) |> copy
    images = unsafe_wrap(Array, sai_local.images,
                         numframes*numsections*framelen, own=false) |> copy
    # Massage data layout
    offset = numsections*framelen 
    imagesnew = Array{Float64,3}(undef,numframes,numsections,framelen)
    for k ∈ 1:numframes
        imagesnew[k,:,:] = images[(k-1)*offset+1:k*offset] |>
                       x -> reshape(x, (framelen, numsections)) |>
                       transpose |>
                       x -> reverse(x, dims=1)
    end
    JSAI(numframes, numsections, framelen, imagesnew, delays, times) 
end


"""
    free(carfacstate::Ptr{CARFACAGCState})

Free the CARFAC state data.
"""
function free(carfacstate::Ptr{CARFACAGCState})
    success = ccall((:carfacagc_free,"/home/dahlbom/research/carfac/ccarfac/src/libcarfac.so"),
                    Cint,
                    (Ptr{CARFACAGCState},),
                    carfacstate)
    if success != 0
        error("Failed to free CARFAC state.")
    end
    return nothing
end


"""
    free(sai::Ptr{SAI})

Free the SAI data.
"""
function free(sai::Ptr{SAI})
    success = ccall((:sai_free, "/home/dahlbom/research/carfac/ccarfac/src/libcarfac.so"),
                    Cint,
                    (Ptr{SAI},),
                    sai)
    if success != 0
        error("Failed to free SAI state.")
    end
    return nothing
end


################################################################################
# Public Interface
################################################################################
OHCParams() = OHCParams(0.1,0.04,bvals_ref)
BMParams() = BMParams(numsections,0.1,0.9,0.2)
IHCParams() = IHCParams(20.0,10.0e-3,0.5e-3,80.0e-6)


"""
    defaultparams(numsections=40)

Generate default parameters for the CARFAC model. Returns parameter structures
for the basilar membrane (bmp), inner hair cells (ihcp), and outer hair cells
(ohcp).
"""
function defaultparams(numsections=40)
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
function calcnap(signal::Array{Float64,1}, fs::T,
                 bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams) where T <: Number
    state = initcarfac(bmp, ihcp, ohcp, length(signal), Float64(fs))
    processblock(state, signal)
    nap, fcs = copynap(state)
    free(state)
    return nap, fcs
end

calcnap(signal::Array{Float64,1}, fs::T) where T <: Number =
        calcnap(signal, fs, defaultparams()...)


"""
    calcsai(signal::Array{Float64,1}, fs::T,
            bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams;
            trigwindowtime = 0.010, advancestep = 0.005,
            numtrigwindows = 4) where T <: Number

Calculate frames of a movie consisting of stabilized auditory
images using the CARFAC model of the auditory periphery. Returns
a three dimensional array, the first index of which corresponds
to frame number, the second to channel, and the third to sample.
"""
function calcsai(signal::Array{Float64,1}, fs::T,
                 bmp::BMParams, ihcp::IHCParams, ohcp::OHCParams;
                 trigwindowtime = 0.010, advancestep = 0.005,
                 numtrigwindows = 4) where T <: Number
    fs = Float64(fs)
    state_ptr = initcarfac(bmp, ihcp, ohcp, length(signal), fs)
    processblock(state_ptr, signal)
    nap, fcs = copynap(state_ptr; julialayout=false)
    carfacstate_local = Base.unsafe_load(state_ptr)
    numsections = carfacstate_local.num_sections
    numpoints = carfacstate_local.block_size
    saip = SAIParams(trigwindowtime, advancestep, numtrigwindows,
                     numsections, numpoints)
    #saip_ref = Base.unsafe_convert(Ptr{SAI}, Base.cconvert(Ptr{SAI}, saip))
    sai_ptr = ccall((:sai_generate, "/home/dahlbom/research/carfac/ccarfac/src/libcarfac.so"),
                    Ptr{SAI},
                    (Ptr{Cdouble}, Cdouble, Ptr{SAIParams}),
                    nap, fs, Ref(saip))
    sai = copysai(sai_ptr)
    free(state_ptr)
    free(sai_ptr)
    return sai, fcs
end

calcsai(signal::Array{Float64,1}, fs::T;
        trigwindowtime = 0.010, advancestep = 0.005,
        numtrigwindows = 4) where T <: Number = 
    calcsai(signal, fs, defaultparams()...; trigwindowtime=trigwindowtime,
            advancestep=advancestep, numtrigwindows=numtrigwindows)

end
