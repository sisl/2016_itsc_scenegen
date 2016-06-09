######################################

function are_collinear(a::AbstractVec, b::AbstractVec, c::AbstractVec, tol::Float64=1e-8)
    # http://mathworld.wolfram.com/Collinear.html
    # if val = 0 then they are collinear
    val = a.x*(b.y-c.y) + b.x*(c.y-a.y)+c.x*(a.y-b.y)
    abs(val) < tol
end

function left_rotate!(arr::AbstractVector, d::Int, n::Int=length(a))
    #=
    Perform a cyclic rotation of the elements in the array, in place
        d = amount to shift
        n = length of array (if we only want to work with first n elements)
    =#

    for i in 1 : gcd(d, n)
    # move i-th values of blocks

        temp = arr[i]
        j = i
        while true
            k = j + d
            if k > n
                k = k - n
            end
            if k == i
                break
            end
            arr[j] = arr[k]
            j = k
        end
        arr[j] = temp
    end

    arr
end
@assert(left_rotate!([1,2,3,4], 1, 4) == [2,3,4,1])
@assert(left_rotate!([1,2,3,4], 2, 4) == [3,4,1,2])
@assert(left_rotate!([1,2,3,4,5,6,7], 1, 4) == [2,3,4,1,5,6,7])
@assert(left_rotate!([1,2,3,4,5,6,7], 2, 4) == [3,4,1,2,5,6,7])

######################################

immutable LineSegment
    a::VecE2
    b::VecE2
end
Base.(:(+))(seg::LineSegment, v::VecE2) = LineSegment(seg.a + v, seg.b + v)
Base.(:(-))(seg::LineSegment, v::VecE2) = LineSegment(seg.a - v, seg.b - v)

function get_polar_angle(seg::LineSegment)
    θ = atan2(seg.b.y - seg.a.y, seg.b.x - seg.a.x)
    if θ < 0.0
        θ += 2π
    end
    θ
end

######################################

function get_signed_area(pts::Vector{VecE2}, npts::Int = length(pts))

    # https://en.wikipedia.org/wiki/Shoelace_formula
    # sign of -1 means clockwise, sign of 1 means counterclockwise

    retval = pts[end].x*pts[1].y - pts[1].x*pts[end].y
    for i in 1 : npts-1
        retval += pts[i].x * pts[i+1].y
        retval -= pts[i+1].x * pts[i].y
    end

    retval / 2
end
function get_edge(pts::Vector{VecE2}, i::Int, npts::Int=length(pts))
    a = pts[i]
    b = i+1 ≤ npts ? pts[i+1] : pts[1]
    LineSegment(a,b)
end

######################################

type ConvexPolygon
    pts::Vector{VecE2} # ordered counterclockwise along polygon boundary s.t. first edge has minimum polar angle in [0,2π]
    npts::Int # number of pts that are currently used (since we preallocate a longer array)

    ConvexPolygon(npts::Int) = new(Array(VecE2, npts), 0)
end

Base.start(poly::ConvexPolygon) = 1
Base.done(poly::ConvexPolygon, i::Int) = i > length(poly)
Base.next(poly::ConvexPolygon, i::Int) = (poly.pts[i], i+1)

Base.length(poly::ConvexPolygon) = poly.npts
Base.isempty(poly::ConvexPolygon) = poly.npts == 0
get_edge(P::ConvexPolygon, i::Int) = get_edge(P.pts, i, P.npts)
get_signed_area(poly::ConvexPolygon) = get_signed_area(poly.pts, poly.npts)
function Base.empty!(poly::ConvexPolygon)
    poly.npts = 0
    poly
end
function Base.copy!(dest::ConvexPolygon, src::ConvexPolygon)
    dest.npts = src.npts
    copy!(dest.pts, 1, src.pts, 1, src.npts)
    dest
end
function Base.push!(poly::ConvexPolygon, v::VecE2)
    poly.pts[poly.npts+1] = v
    poly.npts += 1
    poly
end
function Base.contains(poly::ConvexPolygon, v::VecE2)

    # does NOT include pts on the physical boundary

    previous_side = 0

    for i in 1 : length(poly)
        seg = get_edge(poly, i)
        affine_segment = seg.b - seg.a
        affine_point = v - seg.a
        current_side = round(Int, sign(cross(affine_segment, affine_point))) # sign indicates side
        if current_side == 0
            # outside or over an edge
            return false
        elseif previous_side == 0
            # first segment
            previous_side = current_side
        elseif previous_side != current_side
            # all most be on the same side
            return false
        end
    end
    true
end
Base.contains(poly::ConvexPolygon, v::VecSE2) = contains(poly, convert(VecE2, v))
function Base.show(io::IO, poly::ConvexPolygon)
    @printf(io, "ConvexPolygon: len %d (max %d pts)\n", length(poly), length(poly.pts))
    for i in 1 : length(poly)
        print(io, "\t"); show(io, poly.pts[i])
        print(io, "\n")
    end
end

function ensure_pts_sorted_by_min_polar_angle!(poly::ConvexPolygon, npts::Int=poly.npts)

    @assert(npts ≥ 3)
    @assert(sign(get_signed_area(poly)) == 1) # must be counter-clockwise

    # ensure that edges are sorted by minimum polar angle in [0,2π]

    angle_start = Inf
    index_start = -1
    for i in 1 : npts
        seg = get_edge(poly.pts, i)

        θ = atan2(seg.b.y - seg.a.y, seg.b.x - seg.a.x)
        if θ < 0.0
            θ += 2π
        end

        if θ < angle_start
            angle_start = θ
            index_start = i
        end
    end

    if index_start != 1
        left_rotate!(poly.pts, index_start-1, npts)
    end
    poly
end
function shift!(poly::ConvexPolygon, v::VecE2)
    for i in 1 : length(poly)
        poly.pts[i] += v
    end
    ensure_pts_sorted_by_min_polar_angle!(poly)
    poly
end
function rotate!(poly::ConvexPolygon, θ::Float64)
    for i in 1 : length(poly)
        poly.pts[i] = Vec.rot(poly.pts[i], θ)
    end
    ensure_pts_sorted_by_min_polar_angle!(poly)
    poly
end
function mirror!(poly::ConvexPolygon)
    for i in 1 : length(poly)
        poly.pts[i] = -poly.pts[i]
    end
    ensure_pts_sorted_by_min_polar_angle!(poly)
    poly
end

function minkowksi_sum!(retval::ConvexPolygon, P::ConvexPolygon, Q::ConvexPolygon)

    #=
    For two convex polygons P and Q in the plane with m and n vertices, their Minkowski sum is a
    convex polygon with at most m + n vertices and may be computed in time O (m + n) by a very simple procedure,
    which may be informally described as follows.

    Assume that the edges of a polygon are given and the direction, say, counterclockwise, along the polygon boundary.
    Then it is easily seen that these edges of the convex polygon are ordered by polar angle.
    Let us merge the ordered sequences of the directed edges from P and Q into a single ordered sequence S.
    Imagine that these edges are solid arrows which can be moved freely while keeping them parallel to their original direction.
    Assemble these arrows in the order of the sequence S by attaching the tail of the next arrow to the head of the previous arrow.
    It turns out that the resulting polygonal chain will in fact be a convex polygon which is the Minkowski sum of P and Q.
    =#

    empty!(retval)

    index_P = 1
    index_Q = 1
    θp = get_polar_angle(get_edge(P, index_P))
    θq = get_polar_angle(get_edge(Q, index_Q))

    while index_P ≤ length(P) || index_Q ≤ length(Q)
        # select next edge with minimum polar angle

        if θp == θq
            seg_p = get_edge(P, index_P)
            seg_q = get_edge(Q, index_Q)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg_p.b - seg_p.a + seg_q.b - seg_q.a)
            index_P += 1
            θp = index_P ≤ length(P) ? get_polar_angle(get_edge(P, index_P)) : Inf
            index_Q += 1
            θq = index_Q ≤ length(Q) ? get_polar_angle(get_edge(Q, index_Q)) : Inf
        elseif θp ≤ θq
            seg = get_edge(P, index_P)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg.b - seg.a)
            index_P += 1
            θp = index_P ≤ length(P) ? get_polar_angle(get_edge(P, index_P)) : Inf
        else
            seg = get_edge(Q, index_Q)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg.b - seg.a)
            index_Q += 1
            θq = index_Q ≤ length(Q) ? get_polar_angle(get_edge(Q, index_Q)) : Inf
        end
    end

    ensure_pts_sorted_by_min_polar_angle!(retval)

    retval
end
function minkowksi_difference!(retval::ConvexPolygon, P::ConvexPolygon, Q::ConvexPolygon)

    #=
    The minkowski difference is what you get by taking the minkowski sum of a shape and the mirror of another shape.
    So, your second shape gets flipped about the origin (all of its points are negated).

    The idea is that you do a binary operation on two shapes to get a new shape,
    and if the origin (the zero vector) is inside that shape, then they are colliding.
    =#

    minkowksi_sum!(retval, P, mirror!(Q))
end

function is_colliding(P::ConvexPolygon, Q::ConvexPolygon, temp::ConvexPolygon=ConvexPolygon(length(P) + length(Q)))
    minkowksi_difference!(temp, P, Q)
    contains(temp, VecE2(0,0))
end


function set_to_oriented_bounding_box!(retval::ConvexPolygon, veh::Vehicle)

    # get a bounding box relative to the center of the vehicle

    x = polar(veh.length/2, veh.state.posG.θ)
    y = polar(veh.width/2, veh.state.posG.θ+π/2)

    retval.pts[1] =  x - y
    retval.pts[2] =  x + y
    retval.pts[3] = -x + y
    retval.pts[4] = -x - y
    retval.npts = 4

    ensure_pts_sorted_by_min_polar_angle!(retval)

    retval
end
function set_to_positioned_oriented_bounding_box!(retval::ConvexPolygon, veh::Vehicle)

    # get a bounding box relative to the center of the vehicle

    c = convert(VecE2, get_center(veh))
    x = polar(veh.length/2, veh.state.posG.θ)
    y = polar(veh.width/2, veh.state.posG.θ+π/2)

    retval.pts[1] =  x - y + c
    retval.pts[2] =  x + y + c
    retval.pts[3] = -x + y + c
    retval.pts[4] = -x - y + c
    retval.npts = 4

    ensure_pts_sorted_by_min_polar_angle!(retval)

    retval
end

######################################

function get_future_position(ray::VecSE2, speed::Float64, t::Float64)
    convert(VecE2, ray) + polar(speed*t, ray.θ)
end

function is_colliding(ray::VecSE2, seg::LineSegment)
    v₁ = convert(VecE2, ray) - seg.a
    v₂ = seg.b - seg.a
    v₃ = polar(1.0, ray.θ + π/2)

    denom = dot(v₂, v₃)

    if abs(denom) > 0.0
        t₁ = cross(v₂, v₁) / denom # time for ray (0 ≤ t₁)
        t₂ = dot(v₁, v₃) / denom # time for segment (0 ≤ t₂ ≤ 1)
        return 0.0 ≤ t₁ && 0.0 ≤ t₂ ≤ 1.0
    else
        # denom is zero if the segment and the ray are parallel
        # only collide if they are perfectly aligned
        return are_collinear(ray, seg.a, seg.b)
    end
end
function get_collision_time(ray::VecSE2, seg::LineSegment, ray_speed::Float64)

    o = convert(VecE2, ray)
    v₁ = o - seg.a
    v₂ = seg.b - seg.a
    v₃ = polar(1.0, ray.θ + π/2)

    denom = dot(v₂, v₃)

    if abs(denom) > 0.0
        d₁ = cross(v₂, v₁) / denom # time for ray (0 ≤ t₁)
        t₂ = dot(v₁, v₃) / denom # time for segment (0 ≤ t₂ ≤ 1)
        if 0.0 ≤ d₁ && 0.0 ≤ t₂ ≤ 1.0
            return d₁/ray_speed
        end
    else
        # denom is zero if the segment and the ray are parallel
        # only collide if they are perfectly aligned
        if are_collinear(ray, seg.a, seg.b)
            dist_a = abs2(seg.a - o)
            dist_b = abs2(seg.b - o)
            return sqrt(min(dist_a, dist_b)) / ray_speed
        end
    end

    NaN
end
function get_time_and_dist_of_closest_approach(ray::VecSE2, seg::LineSegment, ray_speed::Float64, if_no_col_skip_eval::Bool=false)

    o = convert(VecE2, ray)
    v₁ = o - seg.a
    v₂ = seg.b - seg.a
    v₃ = polar(1.0, ray.θ + π/2)

    denom = dot(v₂, v₃)

    if abs(denom) > 0.0

        d₁ = cross(v₂, v₁) / denom # time for ray (0 ≤ d₁)
        t₂ = dot(v₁, v₃) / denom # time for segment (0 ≤ t₂ ≤ 1)

        if 0.0 ≤ d₁ && 0.0 ≤ t₂ ≤ 1.0
            t = d₁/ray_speed
            return (t, 0.0)
        else
            # no collision, get time of closest approach to each endpoint
            if !if_no_col_skip_eval
                ray_pt = polar(1.0, ray.θ)
                projA = proj(seg.a - o, ray_pt, VecE2)
                projB = proj(seg.b - o, ray_pt, VecE2)
                tA = max(abs(projA)/ray_speed * sign(dot(ray_pt, projA)), 0.0)
                tB = max(abs(projB)/ray_speed * sign(dot(ray_pt, projB)), 0.0)
                pA = get_future_position(ray, ray_speed, tA)
                pB = get_future_position(ray, ray_speed, tB)
                distA = abs2(seg.a - pA)
                distB = abs2(seg.b - pB)
                if distA < distB
                    return (tA, sqrt(distA))
                else
                    return (tB, sqrt(distB))
                end
            else
                return (Inf, Inf)
            end
        end
    else
        # denom is zero if the segment and the ray are parallel
        # only collide if they are perfectly aligned
        if are_collinear(ray, seg.a, seg.b)
            dist_a = abs2(seg.a - o)
            dist_b = abs2(seg.b - o)
            t = sqrt(min(dist_a, dist_b)) / ray_speed
            return (t, 0.0)
        else
            if !if_no_col_skip_eval
                ray_pt = polar(1.0, ray.θ)
                projA = proj(seg.a - o, ray_pt, VecE2)
                projB = proj(seg.b - o, ray_pt, VecE2)
                tA = max(abs(projA)/ray_speed * sign(dot(ray_pt, projA)), 0.0)
                tB = max(abs(projB)/ray_speed * sign(dot(ray_pt, projB)), 0.0)
                pA = get_future_position(ray, ray_speed, tA)
                pB = get_future_position(ray, ray_speed, tB)
                distA = abs2(seg.a - pA)
                distB = abs2(seg.b - pB)
                if distA < distB
                    return (tA, distA)
                else
                    return (tB, distB)
                end
            else
                return (Inf,Inf)
            end
        end
    end
end

function is_colliding(ray::VecSE2, poly::ConvexPolygon)
    # collides if at least one of the segments collides with the ray

    for i in 1 : length(poly)
        seg = get_edge(poly, i)
        if is_colliding(ray, seg)
            return true
        end
    end
    false
end

######################################


immutable CPAMemory
    vehA::ConvexPolygon # bounding box for vehicle A
    vehB::ConvexPolygon # bounding box for vehicle B
    mink::ConvexPolygon # minkowski bounding box

    CPAMemory() = new(ConvexPolygon(4), ConvexPolygon(4), ConvexPolygon(8))
end
function get_time_and_dist_of_closest_approach(a::Vehicle, b::Vehicle, mem::CPAMemory=CPAMemory())

    set_to_oriented_bounding_box!(mem.vehA, a)
    set_to_oriented_bounding_box!(mem.vehB, b)
    minkowksi_sum!(mem.mink, mem.vehA, mem.vehB)

    rel_pos = convert(VecE2, b.state.posG) - a.state.posG
    rel_velocity = polar(b.state.v, b.state.posG.θ) - polar(a.state.v, a.state.posG.θ)
    ray_speed = hypot(rel_velocity)
    ray = VecSE2(rel_pos, atan2(rel_velocity))

    if contains(mem.mink, convert(VecE2, ray))
        return (0.0, 0.0)
    end

    best_t_CPA = NaN
    best_d_CPA = Inf
    if_no_col_skip_eval = false

    for i in 1 : length(mem.mink)
        seg = get_edge(mem.mink, i)
        t_CPA, d_CPA = get_time_and_dist_of_closest_approach(ray, seg, ray_speed, if_no_col_skip_eval)
        if d_CPA < best_d_CPA
            best_t_CPA = t_CPA
            best_d_CPA = d_CPA
            if d_CPA == 0.0
                if_no_col_skip_eval = true
            end
        end
    end

    (best_t_CPA, best_d_CPA)
end