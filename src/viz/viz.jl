abstract Overlay
function render!(rendermodel::RenderModel, overlay::Overlay, trajdata::Trajdata, frame::Int)
    # does nothing
    rendermodel
end

type LineToCenterline <: Overlay
    target_id::Int # if -1 does it for all
end
function render!(rendermodel::RenderModel, overlay::LineToCenterline, trajdata::Trajdata, frame::Int,
    color::Colorant = colorant"white",
    )

    line_width = 0.5

    if overlay.target_id < 0
        target_ids = carsinframe(trajdata, frame)
    else
        target_ids = [overlay.target_id]
    end

    for carid in target_ids
        veh = get_vehicle(trajdata, carid, frame)
        footpoint = get_footpoint(veh)
        Renderer.add_instruction!(rendermodel, render_line_segment, (veh.state.posG.x, veh.state.posG.y, footpoint.x, footpoint.y, color, line_width))
    end

    rendermodel
end

type FrenetDisplay <: Overlay
    target_id::Int # if -1 does it for all
end
function render!(rendermodel::RenderModel, overlay::FrenetDisplay, trajdata::Trajdata, frame::Int,
    color::Colorant = colorant"white",
    )

    render!(rendermodel, LineToCenterline(overlay.target_id), trajdata, frame)

    if overlay.target_id > 0
        veh = get_vehicle(trajdata, overlay.target_id, frame)
        posF, laneid = project_to_closest_lane(veh.state.posG, trajdata.roadway)

        text_y = -5
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("extind: %.6f", posF.x), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("laneid: %d", laneid), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("d_cl: %.6f", posF.y), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("theta: %.5f", rad2deg(posF.θ)), 10, text_y+=20, 15, color), incameraframe=false)
        text_y+=20
        posF = veh.state.posF
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("extind: %.6f", posF.extind), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("laneid: %d", posF.laneid), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("d_cl: %.6f", posF.t), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("theta: %.5f", rad2deg(posF.ϕ)), 10, text_y+=20, 15, color), incameraframe=false)

        posG = veh.state.posG
        text_y+=20
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("x: %.5f", posG.x), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("y: %.5f", posG.y), 10, text_y+=20, 15, color), incameraframe=false)
        Renderer.add_instruction!(rendermodel, render_text, (@sprintf("t: %.5f", rad2deg(posG.θ)), 10, text_y+=20, 15, color), incameraframe=false)
    end

    rendermodel
end

type SceneValidityDisplay <: Overlay
    # scene_params::SceneExtractParams
    mem::CPAMemory
end
function render!(rendermodel::RenderModel, overlay::SceneValidityDisplay, scene::Scene, structure::SceneStructure,
    color::Colorant = colorant"white",
    )

    # _is_in_bounds = is_in_bounds(overlay.scene, overlay.scene_params)
    _is_collision_free = is_collision_free(scene, overlay.mem)
    _is_there_longitudinal_room = is_there_longitudinal_room(scene)

    text_y = -5
    # Renderer.add_instruction!(rendermodel, render_text, ("is_in_bounds: " * string(_is_in_bounds), 10, text_y+=20, 15, color), incameraframe=false)
    Renderer.add_instruction!(rendermodel, render_text, ("is_collision_free: " * string(_is_collision_free), 10, text_y+=20, 15, color), incameraframe=false)
    Renderer.add_instruction!(rendermodel, render_text, ("is_there_longitudinal_room: " * string(_is_there_longitudinal_room), 10, text_y+=20, 15, color), incameraframe=false)

    rendermodel
end

type SceneBox <: Overlay
    scene_extract::SceneExtractParams
end
render!(rendermodel::RenderModel, overlay::SceneBox, trajdata::Trajdata, frame::Int) = render!(rendermodel, overlay.scene_extract)

type SceneStructureOverlay <: Overlay
    scene_extract::SceneExtractParams
    extractor::Scene
    factors::Vector{Factor}
end
function render!(rendermodel::RenderModel, overlay::SceneStructureOverlay, trajdata::Trajdata, frame::Int)
    scene = pull_subscene(trajdata, frame, overlay.scene_extract, overlay.extractor)
    structure = gen_scene_structure(scene, factors)
    render!(rendermodel, structure, scene)
end

type MinTargetOverlay <: Overlay
end
function render!(rendermodel::RenderModel, overlay::MinTargetOverlay, scene::Scene, structure::SceneStructure)

    min_target = Inf
    target_index = -1
    front_index = -1

    for (factor_index, vehicle_indeces) in structure.factor_assignments
        if factor_index == 2
            rear = scene.vehicles[vehicle_indeces[1]]
            fore = scene.vehicles[vehicle_indeces[2]]
            a_req = get_required_acceleration(rear, fore)
            if a_req < min_target
                min_target = a_req
                target_index = vehicle_indeces[1]
                front_index = vehicle_indeces[2]
            end
        end
    end

    fore = scene.vehicles[front_index]
    rear = scene.vehicles[target_index]
    t_CPA, d_CPA = closest_time_of_approach_and_distance(fore.state.posG, fore.state.v,
                                                         rear.state.posG, rear.state.v)

    color = colorant"white"
    Renderer.add_instruction!(rendermodel, render_text, (@sprintf("value: %.6f", -min_target), 10, 20, 15, color), incameraframe=false)
    Renderer.add_instruction!(rendermodel, render_text, (@sprintf("t_CPA: %.6f", t_CPA), 10, 40, 15, color), incameraframe=false)
    Renderer.add_instruction!(rendermodel, render_text, (@sprintf("d_CPA: %.6f", d_CPA), 10, 60, 15, color), incameraframe=false)

    render!(rendermodel, rear, colorant"red")
    render!(rendermodel, fore, colorant"blue")
    rendermodel
end

type CollisionOverlay <: Overlay
end
function render!(rendermodel::RenderModel, overlay::CollisionOverlay, scene::Scene)

    P = ConvexPolygon(4)
    Q = ConvexPolygon(4)
    M = ConvexPolygon(8)

    has_col = falses(length(scene))

    for (a, vehA) in enumerate(scene)

        set_to_positioned_oriented_bounding_box!(P, vehA)

        for b in a + 1 : length(scene)

            vehB = scene[b]
            set_to_positioned_oriented_bounding_box!(Q, vehB)

            if is_colliding(P, Q, M)
                has_col[a] = has_col[b] = true
            end
        end
    end

    for (a, veh) in enumerate(scene)
        color = has_col[a] ? RGB(1.0,0.0,0.0) : RGB(0.0,0.0,1.0)
        add_instruction!(rendermodel, render_vehicle, (veh.state.posG.x, veh.state.posG.y, veh.state.posG.θ, veh.length, veh.width, color))
    end
end

type Overwash <: Overlay
    color::Colorant
end
function render!(rendermodel::RenderModel, overlay::Overwash, trajdata::Trajdata, frame::Int)
    Renderer.add_instruction!(rendermodel, render_paint, (overlay.color,))
end

###############

function plot_lines(lines::AbstractVector{Vector{VecE2}};
    rand_color::Bool=false,
    highlight_target::Int=0,
    )

    plots = Array(Plots.Linear, length(lines))

    for (i,line) in enumerate(lines)
        arr_x = map(pt->-pt.y, line)
        arr_y = map(pt->-pt.x, line)

        bluePercent = rand_color ? rand(20:100) : 80
        greenPercent = rand_color ? rand(20:100) : 80

        if i == highlight_target
            plots[i] = Plots.Linear(arr_x, arr_y, style=@sprintf("thick, red, mark=none, solid"))
        else
            plots[i] = Plots.Linear(arr_x, arr_y, style=@sprintf("thick, blue!%d!green!%d, mark=none, solid", bluePercent, greenPercent))
        end
    end

    Axis(plots, xlabel=L"-y", ylabel=L"-x", height="8cm", width="25cm", axisEqual=true)
end
function plot_boundaries_and_centerlines(
    boundaries::AbstractVector{Vector{VecE2}},
    centerlines::Dict{AbstractString, Vector{VecE2}};
    rand_color::Bool=false,
    )

    plots = Array(Plots.Linear, length(boundaries))

    for (i,line) in enumerate(boundaries)
        arr_x = map(pt->-pt.y, line)
        arr_y = map(pt->-pt.x, line)

        bluePercent = rand_color ? rand(20:100) : 80
        greenPercent = rand_color ? rand(20:100) : 80

        plots[i] = Plots.Linear(arr_x, arr_y, style=@sprintf("thick, blue!%d!green!%d, mark=none, solid", bluePercent, greenPercent))
    end

    for (name,line) in centerlines
        arr_x = map(pt->-pt.y, line)
        arr_y = map(pt->-pt.x, line)
        push!(plots, Plots.Linear(arr_x, arr_y, style="red!80, mark=none, solid"))
    end

    Axis(plots, xlabel=L"-y", ylabel=L"-x", height="18cm", width="100cm", axisEqual=true)
end

###############

center_camera_on!(rendermodel::RenderModel, veh::Vehicle) = camera_set!(rendermodel, convert(VecE2, veh.state.posG), 3.0)
center_camera_on!(rendermodel::RenderModel, trajdata::Trajdata, carid::Int, frame::Int) = center_camera_on!(rendermodel, get_vehicle(trajdata, carid, frame))
function scene_center(scene::Scene)
    center = VecE2(0.0,0.0)
    for veh in scene
        center += convert(VecE2, veh.state.posG)
    end
    center / length(scene)
end
function scene_center(rec::SceneRecord, past_frame::Int)
    center = VecE2(0.0,0.0)
    n_vehicles = rec.n_vehicles[1-past_frame]
    for i in 1 : n_vehicles
        center += convert(VecE2, rec[i,past_frame].state.posG)
    end
    center / n_vehicles
end

function render!(rendermodel::RenderModel, roadway::Roadway;
    color::Colorant = RGB(0x88, 0x88, 0xFF),
    line_width::Real=1.5 # [ft]
    )

    for line in roadway.boundaries
        pts = Array(Float64, 2, length(line))
        for (i,v) in enumerate(line)
            pts[1,i] = v.x
            pts[2,i] = v.y
        end
        Renderer.add_instruction!(rendermodel, render_line, (pts, color, line_width))
    end

    # color_centerline = RGB(0.1, 0.9, 0.9)
    # for line in roadway.centerlines
    #     pts = Array(Float64, 2, length(line))
    #     for (i,v) in enumerate(line)
    #         pts[1,i] = v.pos.x
    #         pts[2,i] = v.pos.y
    #     end
    #     Renderer.add_instruction!(rendermodel, render_line, (pts, color_centerline, line_width))
    # end

    rendermodel
end
function render!(rendermodel::RenderModel, veh::Vehicle,
    color::Colorant = RGB(rand(), rand(), rand())
    )


    v = clamp(veh.state.v, BOUNDS_V_LO, BOUNDS_V_HI)
    speed_t = 0.5*(v - BOUNDS_V_LO) / (BOUNDS_V_HI - BOUNDS_V_LO) + 0.4
    # speed_t = 0.9
    color = RGB(speed_t,speed_t,speed_t)
    color_arrow = RGB(0.0,0.0,0.0)

    # println("speed_t: ", speed_t)
    # println("length:  ", veh.length)
    # println("width:   ", veh.width)
    # println("color:   ", color)
    # println("color_arrow:   ", color_arrow)


    add_instruction!(rendermodel, render_vehicle, (veh.state.posG.x, veh.state.posG.y, veh.state.posG.θ, veh.length, veh.width, color, color, color_arrow))
    rendermodel
end
function render!(rendermodel::RenderModel, poly::ConvexPolygon, color::Colorant, line_width::Real)

    pts = Array(Float64, 2, length(poly))
    for (i,p) in enumerate(poly)
        pts[1,i] = p.x
        pts[2,i] = p.y
    end

    Renderer.add_instruction!(rendermodel, render_closed_line, (pts, color, line_width))

    rendermodel
end
render!(rendermodel::RenderModel, scene_extract::SceneExtractParams) =
    render!(rendermodel, scene_extract.box, colorant"red", 1.0)

function render!(rendermodel::RenderModel, trajdata::Trajdata, frame::Int, carid_target::Int=0)

    for carid in carsinframe(trajdata, frame)
        veh = get_vehicle(trajdata, carid, frame)
        if carid == carid_target
            render!(rendermodel, veh, COLOR_CAR_EGO)
        else
            render!(rendermodel, veh, COLOR_CAR_OTHER)
        end
    end

    rendermodel
end
function render!(rendermodel::RenderModel, trajdata::Trajdata, frame::Int, color::Colorant)

    for carid in carsinframe(trajdata, frame)
        veh = get_vehicle(trajdata, carid, frame)
        render!(rendermodel, veh, color)
    end

    rendermodel
end
function render!(rendermodel::RenderModel, scene::Scene)

    for i in 1 : length(scene)
        render!(rendermodel, scene.vehicles[i], COLOR_CAR_OTHER)
    end

    rendermodel
end
function render!(rendermodel::RenderModel, structure::SceneStructure, scene::Scene)

    line_len = 8.0 # ft
    line_width = 1.5 # ft
    circle_radius = 2.0 # ft

    color_line = colorant"white"
    color_road = RGBA(0x52/0xFF,0xE3/0xFF,0xF6/0xFF)
    color_following = RGBA(0xA7/0xFF,0xEC/0xFF,0x21/0xFF)
    color_neighbor = RGBA(0xFF/0xFF,0x00/0xFF,0x7F/0xFF)

    # render all links first
    for (factor_index, vehicle_indeces) in structure.factor_assignments
        @assert(length(vehicle_indeces) > 0)
        if length(vehicle_indeces) == 2
            # draw a line between the vehicle centers

            pos1 = get_center(scene.vehicles[vehicle_indeces[1]])
            pos2 = get_center(scene.vehicles[vehicle_indeces[2]])
            center = lerp(pos1, pos2, 0.5)

            Renderer.add_instruction!(rendermodel, render_line_segment, (pos1.x, pos1.y, pos2.x, pos2.y, color_line, line_width))
        end
    end

    for (factor_index, vehicle_indeces) in structure.factor_assignments
        @assert(length(vehicle_indeces) > 0)
        if length(vehicle_indeces) == 1
            # draw a line from the vehicle center to the bottom right

            veh = scene.vehicles[vehicle_indeces[1]]
            pos = get_center(veh)

            Renderer.add_instruction!(rendermodel, render_circle, (pos.x, pos.y, circle_radius, color_road))
        elseif length(vehicle_indeces) == 2
            # draw a line between the vehicle centers

            pos1 = get_center(scene.vehicles[vehicle_indeces[1]])
            pos2 = get_center(scene.vehicles[vehicle_indeces[2]])
            center = lerp(pos1, pos2, 0.5)
            color = factor_index == 2 ? color_following : color_neighbor

            Renderer.add_instruction!(rendermodel, render_circle, (center.x, center.y, circle_radius, color))
        else
            error("not implemented")
        end
    end
end
function render!(rendermodel::RenderModel, rec::SceneRecord, past_frame::Int)

    for i in 1 : rec.n_vehicles[1-past_frame]
        render!(rendermodel, rec[i,past_frame], COLOR_CAR_OTHER)
    end

    rendermodel
end

function render_scene!(
    trajdata::Trajdata,
    frame::Int,
    carid_target::Int=0;

    camerazoom::Real=3.0, # [pix/m]
    canvas_width::Integer=1100,
    canvas_height::Integer=800,
    rendermodel::RenderModel=RenderModel(),
    overlays::Vector{Overlay} = Overlay[],
    )

    s = CairoRGBSurface(canvas_width, canvas_height)
    ctx = creategc(s)
    clear_setup!(rendermodel)

    render!(rendermodel, trajdata.roadway)
    render!(rendermodel, trajdata, frame, carid_target)

    for overlay in overlays
        render!(rendermodel, overlay, trajdata, frame)
    end

    if carid_target != 0
        center_camera_on!(rendermodel, trajdata, carid_target, frame)
    else
        camera_set!(rendermodel, trajdata.roadway.boundaries[1][5], 3.0)
    end
    camera_setzoom!(rendermodel, camerazoom)

    render(rendermodel, ctx, canvas_width, canvas_height)
    s
end
function render_scene!(
    trajdata::Trajdata,
    frame::Int,
    scene_extract::SceneExtractParams,
    extractor::Scene;

    camerazoom::Real=3.0, # [pix/m]
    canvas_width::Integer=1100,
    canvas_height::Integer=800,
    rendermodel::RenderModel=RenderModel(),
    overlays::Vector{Overlay} = Overlay[],
    )

    s = CairoRGBSurface(canvas_width, canvas_height)
    ctx = creategc(s)
    clear_setup!(rendermodel)

    render!(rendermodel, trajdata.roadway)
    render!(rendermodel, trajdata, frame, RGBA(0.8,0.8,0.8,0.5))
    render!(rendermodel, pull_subscene(trajdata, frame, scene_extract, extractor))

    for overlay in overlays
        render!(rendermodel, overlay, trajdata, frame)
    end

    camera_set!(rendermodel, scene_extract.center, camerazoom)

    render(rendermodel, ctx, canvas_width, canvas_height)
    s
end
function render_scene!(
    scene::Scene,
    structure::Union{Void, SceneStructure} = nothing;

    camerazoom::Real=5.0, # [pix/m]
    canvas_width::Integer=600,
    canvas_height::Integer=600,
    rendermodel::RenderModel=RenderModel(),
    overlays::Vector{Overlay} = Overlay[],
    s::CairoSurface = CairoRGBSurface(canvas_width, canvas_height),
    )

    ctx = creategc(s)
    clear_setup!(rendermodel)

    render!(rendermodel, get_roadway(scene))
    render!(rendermodel, scene)
    # render!(rendermodel, CollisionOverlay(), scene)

    if isa(structure, SceneStructure)
        Renderer.add_instruction!(rendermodel, render_paint, (RGBA(0.0,0.0,0.0,0.7),)) # overwash
        render!(rendermodel, structure, scene)
        for overlay in overlays
            render!(rendermodel, overlay, scene, structure::SceneStructure)
        end
    end

    camera_set!(rendermodel, scene_center(scene), camerazoom)

    render(rendermodel, ctx, canvas_width, canvas_height)
    s
end
function render_scene_overlay!(
    trajdata::Trajdata,
    roadway::Roadway,
    frame::Int,
    carid_target::Int=0;

    camerazoom::Real=3.0, # [pix/m]
    canvas_width::Integer=1100,
    canvas_height::Integer=800,
    rendermodel::RenderModel=RenderModel(),
    )

    scene = get_scene(trajdata, roadway, frame, carid_target, -100.0, 100.0)

    s = CairoRGBSurface(canvas_width, canvas_height)
    ctx = creategc(s)
    clear_setup!(rendermodel)

    render!(rendermodel, trajdata, frame, RGBA(0.8,0.8,0.8,0.5))
    render!(rendermodel, scene)

    center_camera_on!(rendermodel, scene.vehicles[1])
    camera_setzoom!(rendermodel, camerazoom)

    render(rendermodel, ctx, canvas_width, canvas_height)
    s
end

function render_rec!(
    rec::SceneRecord,
    past_frame::Int=0;

    camerazoom::Real=5.0, # [pix/m]
    canvas_width::Integer=600,
    canvas_height::Integer=600,
    rendermodel::RenderModel=RenderModel(),
    overlays::Vector{Overlay} = Overlay[],
    s::CairoSurface = CairoRGBSurface(canvas_width, canvas_height),
    )

    ctx = creategc(s)
    clear_setup!(rendermodel)

    render!(rendermodel, get_roadway(rec))
    render!(rendermodel, rec, past_frame)

    color = RGB(1.0,1.0,0.0)
    line_width = 2.0
    for i in 1 : rec.n_vehicles[1-past_frame]
        veh = rec[i, past_frame]
        footpoint = get_footpoint(veh)
        Renderer.add_instruction!(rendermodel, render_line_segment, (veh.state.posG.x, veh.state.posG.y, footpoint.x, footpoint.y, color, line_width))
    end

    camera_set!(rendermodel, scene_center(rec, past_frame), camerazoom)

    render(rendermodel, ctx, canvas_width, canvas_height)
    s
end

