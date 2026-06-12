# n13_gks.jl — O-N13-1 (iter 47): LF/GKS analysis of the Z4c constraint-row
# boundary realizations on the linearized constraint subsystem.
#
# Subsystem (1010.0523v2 eq:sys-theta-Z; linear, flat, normal incidence;
# Z4c damping kappa1):
#   dt Theta = -dx Z - k1 (2 + k2) Theta
#   dt Z     = -dx Theta - k1 Z
# Characteristics: w+- = Theta +- Z advect at speed +-1 (w+ toward the
# boundary at x = 0 is OUTGOING).
#
# Boundary-row candidates at x = 0 (continuum form; sigma/r are the
# Sommerfeld falloff weights; the AthenaK kernel applies rows at the face
# node with centered stencils reading order-2 extrapolated ghosts):
#   :somm   dt q = -dx q - sigma_q q/r on BOTH q = Theta, Z (sigma 1, 1)
#   :v1     same + the constraint sink: sigma_Z = 3 (cpbc=1)
#   :adv2   Z-row with ONE-SIDED inward dx and sigma_Z = 5 (cpbc=2 model:
#           -dZ - 2Z/r on the Z-channel TWICE the falloff), Theta somm
#   :bjorhus  face node: dt w- = 0 (incoming-mode replacement),
#           dt w+ = -dx w+ - w+/r (one-sided)  [the textbook maximally
#           dissipative realization]
#
# Tools: (1) frozen semi-discrete operator spectrum (GKS instability shows
# as Re(lambda) growth ~ 1/h under refinement); (2) packet-absorption: a
# clean outgoing constraint packet hits the boundary; absorption =
# 1 - (post-passage energy)/(initial energy); (3) continuum reflection
# |R(omega)| for the row family by 2x2 mode matching.
using LinearAlgebra

export gks_operator, gks_spectrum_maxre, gks_packet_absorption, gks_reflection

# 6th-order centered first derivative + order-2 extrapolated ghosts at the
# right edge (the AthenaK face structure); LEFT edge: clean characteristic
# inflow (w- = 0 there is enforced by the row; packet starts interior)
function _D6_rows(n)
    D = zeros(n, n)
    c = [-1, 9, -45, 0, 45, -9, 1] ./ 60.0
    for i in 4:n-3
        for (k, off) in enumerate(-3:3)
            D[i, i+off] += c[k]
        end
    end
    # near-edge rows: 2nd-order centered (the AthenaK boundary kernel uses
    # order-2 stencils in the boundary region)
    for i in (1, 2, 3, n-2, n-1, n)
        if i == 1
            D[i, 1:3] .= [-1.5, 2.0, -0.5]          # one-sided (left inflow)
        elseif i == n
            D[i, n-2:n] .= [0.5, -2.0, 1.5]          # one-sided fallback
        else
            D[i, i-1] = -0.5; D[i, i+1] = 0.5
        end
    end
    D
end

"""assemble the frozen semi-discrete operator L (2n x 2n, state [Theta; Z])
for a boundary treatment; h = grid spacing, r = boundary radius (falloff),
k1/k2 = Z4c damping. Face = LAST node (x = 0). `bc` in
(:somm, :v1, :adv2, :bjorhus). Ghost-extrapolation effects are modeled by
the one-sided/centered row choices above (the kernel's order-2 region)."""
function gks_operator(n::Int, h::Float64, bc::Symbol;
                      r::Float64 = 41.0, k1::Float64 = 0.02,
                      k2::Float64 = 0.0)
    D = _D6_rows(n) ./ h
    L = zeros(2n, 2n)
    T = 1:n; Z = n+1:2n
    # interior: dt Theta = -D Z - k1(2+k2) Theta ; dt Z = -D Theta - k1 Z
    L[T, Z] .= -D
    L[T, T] .-= k1*(2 + k2)*I(n)
    L[Z, T] .= -D
    L[Z, Z] .-= k1*I(n)
    # left edge: incoming-mode kill (clean inflow): dt w-|_1 = 0 keeps the
    # left boundary passive: rows for node 1 replaced by
    # dt Theta_1 = dt Z_1 = (1/2) dt w+_1, dt w+_1 = +D1(one-sided) w+
    d1 = zeros(n); d1[1:3] .= [-1.5, 2.0, -0.5] ./ h   # toward interior
    wp_row_T = zeros(2n); wp_row_T[T] .= -d1; wp_row_T[Z] .= -d1
    # w+ = Theta + Z moves RIGHT: at the LEFT edge it is incoming-from-
    # outside; w- = Theta - Z is outgoing there. Keep the left edge
    # transparent for w- and kill incoming w+: dt w+|_1 = 0,
    # dt w-|_1 = +D w- (one-sided; speed -1 means dt w- = +dx w-)
    wm_row = zeros(2n); wm_row[T] .= d1; wm_row[Z] .= -d1
    L[1, :] .= 0.5 .* wm_row                 # dt Theta_1
    L[n+1, :] .= -0.5 .* wm_row              # dt Z_1 (w+ contributes 0)
    # right (outer) face rows by candidate
    dn = zeros(n); dn[n-2:n] .= [0.5, -2.0, 1.5] ./ h      # one-sided inward
    dc = zeros(n); dc[n-1] = -0.5/h; dc[n] = 0.5/h         # centered w/ghost
    # ghost value g_{n+1} = 2 q_n - q_{n-1} (order-2 extrapolation):
    # centered derivative with extrapolated ghost == one-sided 2nd-order?
    # (q_{n+1} - q_{n-1})/2h with q_{n+1} = 2q_n - q_{n-1} ->
    # (2q_n - 2q_{n-1})/2h = (q_n - q_{n-1})/h  (FIRST-order one-sided)
    dg = zeros(n); dg[n-1] = -1.0/h; dg[n] = 1.0/h
    if bc == :somm || bc == :v1
        sZ = bc == :v1 ? 3.0 : 1.0
        L[n, :] .= 0; L[2n, :] .= 0
        L[n, T] .= -dg; L[n, n] -= 1.0/r                    # Theta row
        L[2n, Z] .= -dg; L[2n, 2n] -= sZ/r                  # Z row
    elseif bc == :adv2
        L[n, :] .= 0; L[2n, :] .= 0
        L[n, T] .= -dg; L[n, n] -= 1.0/r
        L[2n, Z] .= -dn; L[2n, 2n] -= 2.0/r                 # one-sided, 2Z/r
    elseif bc == :bjorhus
        # dt w-|_n = 0 ; dt w+|_n = -D(one-sided) w+ - w+/r
        wp = zeros(2n); wp[T] .= -dn; wp[Z] .= -dn
        wp[n] -= 1.0/r; wp[2n] -= 1.0/r
        L[n, :] .= 0.5 .* wp                                 # dt Theta_n
        L[2n, :] .= 0.5 .* wp                                # dt Z_n
    else
        error("unknown bc $bc")
    end
    L
end

"""max real part of the spectrum (frozen stability indicator)."""
gks_spectrum_maxre(n, h, bc; kw...) =
    maximum(real, eigvals(gks_operator(n, h, bc; kw...)))

"""packet test: outgoing w+ Gaussian hits the boundary; returns
(absorption, refl) where refl = max |w-| generated / peak |w+| in, and
absorption = 1 - E_after/E_before (E = sum Theta^2 + Z^2)."""
function gks_packet_absorption(n::Int, h::Float64, bc::Symbol;
                               r::Float64 = 41.0, k1::Float64 = 0.0,
                               k2::Float64 = 0.0, cfl::Float64 = 0.2)
    L = gks_operator(n, h, bc; r=r, k1=k1, k2=k2)
    x = collect(range(-(n-1)*h, 0.0; length=n))
    x0, sg = -(n-1)*h*0.6, (n-1)*h*0.06
    wp0 = exp.(-((x .- x0)./sg).^2)
    u = vcat(0.5 .* wp0, 0.5 .* wp0)         # Theta = Z = w+/2, w- = 0
    dt = cfl*h
    nt = ceil(Int, ((abs(x0) + 8sg)/1.0)/dt) # travel + clearance
    E0 = sum(abs2, u)
    reflmax = 0.0
    for _ in 1:nt
        k1v = L*u; k2v = L*(u .+ dt/2 .* k1v)
        k3v = L*(u .+ dt/2 .* k2v); k4v = L*(u .+ dt .* k3v)
        u .+= dt/6 .* (k1v .+ 2k2v .+ 2k3v .+ k4v)
        wm = u[1:n] .- u[n+1:2n]
        reflmax = max(reflmax, maximum(abs.(wm)))
    end
    (1 - sum(abs2, u)/E0, reflmax)
end

"""continuum reflection |R| at s = i*omega for the (sigma_T, sigma_Z)
advection-row family (kappa1 = 0, normal incidence): interior modes
Theta = A e^{s x} + B e^{-s x}, Z = -A e^{s x} + B e^{-s x} (A = incoming
w- amplitude, B = outgoing w+); both rows hold in least-squares sense at
the node — return the two single-row solutions' |A/B| (their spread is the
scheme ambiguity) as a (min, max) tuple."""
function gks_reflection(omega::Float64, sigma_T::Float64, sigma_Z::Float64;
                        r::Float64 = 41.0)
    s = im*omega
    # row Theta: s(A+B) = -(sA - sB) - sigma_T (A+B)/r -> A(2s + sigma_T/r)
    #            = -B sigma_T/r
    RT = abs(-sigma_T/r / (2s + sigma_T/r))
    # row Z: s(-A+B) = -(-sA - sB)... -> -2sA = sigma_Z (A - B)/r
    RZ = abs(sigma_Z/r / (2s + sigma_Z/r))
    (min(RT, RZ), max(RT, RZ))
end
