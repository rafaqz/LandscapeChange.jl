using LandscapeChange, Test
NV = NamedVector

k = (:native, :cleared, :abandoned, :urban, :forestry, :water)

# to = from
transitions = NV(
    native    = NV(native=true,  cleared=false, abandoned=false, urban=false, forestry=false,  water=false),
    cleared   = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    abandoned = NV(native=false, cleared=true,  abandoned=true,  urban=false,  forestry=false,  water=false),
    urban     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=false),
    forestry  = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=true,   water=false),
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=true),
)
# from = to
reversed = LandscapeChange.reverse_transitions(transitions) 
indirect = LandscapeChange.indirect_transitions(transitions)
reversed_indirect = LandscapeChange.reverse_transitions(indirect) 
force = map(_ -> false, transitions)

@testset "generate_bitmasks" begin
    @test LandscapeChange.generate_bitmasks(transitions) == NV(
        native    = NV(native=true , cleared=false, abandoned=false, urban=false, forestry=false, water=false),
        cleared   = NV(native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false),
        abandoned = NV(native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false),
        urban     = NV(native=false, cleared=false, abandoned=false, urban=true , forestry=false, water=false),
        forestry  = NV(native=false, cleared=false, abandoned=false, urban=false, forestry=true , water=false),
        water     = NV(native=false, cleared=false, abandoned=false, urban=false, forestry=false, water=true),
    )
end

@testset "reverse_transitions" begin
    reversed = LandscapeChange.reverse_transitions(transitions) 
    @test reversed == NV(;
        native    = NV(native = true , cleared = true , abandoned = false, urban = true , forestry = true, water = true),
        cleared   = NV(native = false, cleared = true , abandoned = true , urban = true , forestry = true, water = true),
        abandoned = NV(native = false, cleared = true , abandoned = true , urban = true , forestry = true, water = true),
        urban     = NV(native = false, cleared = false, abandoned = false, urban = true , forestry = false, water = true),
        forest    = NV(native = false, cleared = false, abandoned = false, urban = false, forestry = true, water = false),
        water     = NV(native = false, cleared = false, abandoned = false, urban = false, forestry = false, water = true),
    )
end

@testset "merge_all" begin
    a = NV(native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
    b = NV(native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=true)
    @test LandscapeChange.merge_all(a, b, reversed, reversed_indirect) ==
        NV(native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=true)

    a = NV(native=true, cleared=false, abandoned=true , urban=false, forestry=true, water=false)
    b = NV(native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=true)
    @test LandscapeChange.merge_all(a, b, reversed, reversed_indirect) ==
        NV(native=false, cleared=true , abandoned=false, urban=false, forestry=true, water=true)
end

@testset "simplify_end!" begin
    @testset "Contiguous category is kept in end uncertainty" begin
        end_timeline = NV.([
            (native=true , cleared=false, abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=true, urban=false)
            (native=true , cleared=true,  abandoned=true, urban=false)
            (native=false, cleared=true , abandoned=true, urban=true)
        ])
        @test LandscapeChange.simplify_end!(end_timeline) == NV.([
            (native=true , cleared=false, abandoned=false, urban=false)
            (native=false, cleared=true , abandoned=false, urban=false)
            (native=false, cleared=true , abandoned=false, urban=false)
            (native=false, cleared=true , abandoned=false, urban=false)
            (native=false, cleared=true , abandoned=false, urban=false)
        ])
    end
    @testset "No change without contiguity" begin
        end_timeline = NV.([
            (native=true , cleared=false, abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=true, urban=false)
            (native=true , cleared=false, abandoned=true, urban=false)
            (native=false, cleared=true , abandoned=true, urban=true)
        ])
        @test LandscapeChange.simplify_end!(end_timeline) == NV.([
            (native=true , cleared=false, abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=false, urban=false)
            (native=false, cleared=true,  abandoned=true, urban=false)
            (native=true , cleared=false, abandoned=true, urban=false)
            (native=false, cleared=true , abandoned=true, urban=true)
        ])
    end

    timeline = NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (true, true, false, false, false, false)
        (true, true, false, true, false, false)
        (true, true, false, true, false, false)
        (true, true, true, true, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, true, true, false, false, false)
        (false, true, true, true, false, false)
        (false, true, true, true, true, true)
    ])
    simplified = LandscapeChange.simplify_end!(timeline)
    @test simplified == NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (true, true, false, false, false, false)
        (true, true, false, false, false, false)
        (true, true, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
    ])

    timeline = LandscapeChange.LOG[end]
    simplified = LandscapeChange.simplify_end!(timeline)
end


@testset "apply_transitions!" begin
    timeline = NV.([
        (native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=true , cleared=false, abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=true , forestry=false, water=false)
        (native=true , cleared=false, abandoned=false, urban=false, forestry=false, water=false)
        (native=true , cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=true , forestry=false, water=false)
        (native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=true , cleared=false, abandoned=false, urban=true , forestry=false, water=false)
        (native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=false, cleared=false, abandoned=false, urban=false, forestry=false, water=false)
    ])
    errors = [NV(native=false, cleared=false, abandoned=false, urban=false, forestry=false, water=false) for _ in timeline]

    timeline1 = copy(timeline); LandscapeChange.apply_transitions!(timeline1, reversed, reversed_indirect)
    timeline1 .== NV.([
        (native=false, cleared=false, abandoned=false, urban=false, forestry=false, water=false)
        (native=true , cleared=false, abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=false, forestry=false, water=false)
        (native=false, cleared=true , abandoned=false, urban=true , forestry=false, water=false)
        (native=false, cleared=false, abandoned=true , urban=false, forestry=false, water=false)
        (native=false, cleared=false, abandoned=false, urban=true , forestry=false, water=false)
        (native=false, cleared=false, abandoned=false, urban=true , forestry=false, water=false)
        (native=false, cleared=false, abandoned=false, urban=true , forestry=false, water=false)
    ])

    # rev = lastindex(timeline):-1:1
    # timeline2 = copy(timeline); LandscapeChange._apply_transitions!(view(timeline2, rev), transitions, indirect, force)
    # timeline2 .== NV.([
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=false, cleared=false, abandoned=true , urban=false)
    #     (native=false, cleared=false, abandoned=false, urban=false)
    # ])

    # timeline_merged = map(.|, timeline1, timeline2)
    # timeline1 = copy(timeline); errors1 = copy(errors); LandscapeChange.apply_transitions!(timeline1, errors1, reversed, reversed_indirect)
    # timeline2 = copy(timeline_merged); errors2 = copy(errors); LandscapeChange.apply_transitions!(view(timeline2, rev), view(errors2, rev), transitions, indirect)


    # timeline = NV.([
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=false, cleared=true,  abandoned=false, urban=false)
    #     (native=false, cleared=false,  abandoned=true, urban=false)
    #     (native=true , cleared=false, abandoned=false, urban=false)
    #     (native=false, cleared=true , abandoned=false, urban=false)
    # ])

    # timeline1 = copy(timeline); LandscapeChange.apply_transitions!(timeline1, reversed, reversed_indirect, force)
    # timeline2 = copy(timeline); LandscapeChange.apply_transitions!(view(timeline2, rev), transitions, indirect, force)
    # map(.|, timeline1, timeline2)
end


@testset "cross_validate_timeline!" begin
    timeline = NV{k}.([
        (true, false, false, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, true, false, false)
        (true, false, true, false, false, false)
        (true, false, true, false, false, false)
    ])
    simplified, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)
    @test simplified == NV{k}.([
        (true, false, false, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
    ])

    timeline = NV{k}.([
        (true, false, false, false, false, false)
        (false, true, false, false, false, false)
        (false, false, false, true, false, false)
        (false, true,  true, false, false, false)
        (false, false, false, true, false, false)
    ])

    # Simplify
    simplified, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=false)
    @test simplified == NV{k}.([
        (true, false, false, false, false, false)
        (false, true, false, false, false, false)
        (false, true, false, true, false, false)
        (false, true, false, true, false, false)
        (false, false, false, true, false, false)
    ])

    # Cull 
    culled, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)
    # cleared and urban cannot be resolved for two timesteps
    @test culled == NV{k}.([
        (true, false, false, false, false, false)
        (false, true, false, false, false, false)
        (false, true, false, true, false, false)
        (false, true, false, true, false, false)
        (false, false, false, true, false, false)
    ])

    # Force urban and native
    force = NV{k}(map(x -> x in (:native, :urban), propertynames(transitions)))
    forced, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true, force)
    @test forced == NV{k}.([
        (true, false, false, false, false, false)
        (false, true, false, false, false, false)
        (false, false, false, true, false, false)
        (false, false, false, true, false, false)
        (false, false, false, true, false, false)
    ])

    timeline = NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (false, true, true, true, false, false)
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (false, false, false, false, false, false)
        (true, false, true, false, false, false)
        (false, false, false, false, false, false)
        (false, false, false, false, false, false)
        (false, false, false, false, false, false)
        (true, false, true, false, false, false)
        (true, false, true, true, true, true)
        (true, true, true, false, false, false)
        (false, true, true, true, false, false)
        (false, false, false, false, false, false)
    ])
    
    expected, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)
    expected
    expected == NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (true,  true,  true, false, false, false)
        (true,  true,  true, false, false, false)
        (true,  true,  true, false, false, false)
        (true,  true,  true, false, false, false)
        (true,  true,  true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
        (false, false, true, false, false, false)
    ])
    map(expected) do +
        map(first(expected)) do _
            rand([true, false, false, false])
        end
    end

    timeline = NV{k}.([
        (false, true, false, false, true, false)
        (false, false, false, false, true, false)
        (false, false, true, true, false, true)
        (false, true, false, true, false, false)
        (false, false, false, false, false, true)
        (false, false, true, false, false, false)
        (false, true, false, false, true, true)
        (true, true, false, false, true, true)
        (true, true, true, false, false, false)
        (false, true, true, false, false, false)
        (false, true, false, true, true, false)
        (false, true, false, false, false, false)
        (false, false, true, false, false, false)
        (true, false, true, true, true, false)
        (true, false, false, false, true, false)
    ])
    transitions
    expected, _ = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)
    expected

    expected == NV{k}.([
     # (:native, :cleared, :abandoned, :urban, :forestry, :water)
        (false, false, true, false, true, false)  needed = (false, false, false, false, true, false)
        (false, false, true, false, true, false)  needed = (false, false, false, false, true, false)
        (false, false, true, false, true, false)    needed = (false, false, false, false, true, false)
        (false, true, false, false, true, false)   needed = (false, false, false, false, true, false)
        (false, true, true, false, true, true)  needed = (false, false, false, false, false, true)
        (false, false, true, false, true, false)  needed = (false, false, true, false, false, true)
        (false, true, false, false, true, true)    needed = (false, false, false, false, false, false)
        (false, true, false, false, true, true)
        (false, true, true, false, true, false)
        (false, true, true, false, true, false)
        (false, true, false, true, true, false)
        (false, true, false, false, true, false)
        (false, false, true, false, true, false)
        (false, false, true, false, true, false)
        (false, false, false, false, true, false)
    ])

    timeline = NV{k}.([
        (false, false, false, true, false, false)
        (false, false, false, true, false, false)
        (false, true, false, true, false, false)
        (false, false, false, false, true, false)
        (false, false, false, false, false, true)
        (false, false, true, true, false, false)
        (false, false, false, false, false, false)
        (false, false, false, false, false, true)
        (true, false, false, false, false, false)
        (false, false, true, false, false, true)
        (false, false, false, false, true, false)
        (false, false, false, false, false, false)
        (false, false, false, true, false, false)
        (true, false, false, false, true, false)
        (true, true, false, true, false, true)
    ])
    expected = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)[1]
    (:native, :cleared, :abandoned, :urban, :forestry, :water)
    expected == NV{k}.([
        (true, false, false, true, false, false)
        (true, false, false, true, false, false)
        (true, false, false, true, false, false)
        (true, false, false, true, true, false)
        (true, false, false, true, false, true)
        (true, false, false, true, false, true)
        (true, false, false, false, false, true)
        (true, false, false, false, false, true)
        (true, false, false, false, true, true)
        (false, false, false, false, true, true)
        (false, false, false, false, true, true)
        (false, false, false, false, true, true)
        (false, false, false, true, true, true)
        (false, false, false, true, true, true)
        (false, false, false, true, false, true)
    ])

    timeline = NV{k}.([
        (true, false, false, false, false, true)
        (false, false, false, false, true, false)
        (true, true, false, false, false, true)
        (false, false, false, false, true, false)
        (true, true, true, false, false, false)
        (false, true, true, false, true, false)
        (false, false, false, false, false, true)
        (true, false, true, false, true, false)
        (true, false, false, false, false, true)
        (false, false, false, false, false, false)
        (false, true, false, true, true, false)
        (false, false, true, false, true, false)
        (false, false, false, false, false, true)
        (false, false, false, false, true, true)
        (true, true, false, true, false, false)
    ])
    expected = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true)[1]
    expected == NV{k}.([
        (true, false, false, false, false, true)
        (false, false, false, false, true, false)
        (true, true, false, false, false, true)
        (false, false, false, false, true, false)
        (true, true, true, false, false, false)
        (false, true, true, false, true, false)
        (false, false, false, false, false, true)
        (true, false, true, false, true, false)
        (true, false, false, false, false, true)
        (false, false, false, false, false, false)
        (false, true, false, true, true, false)
        (false, false, true, false, true, false)
        (false, false, false, false, false, true)
        (false, false, false, false, true, true)
        (true, true, false, true, false, false)
    ])

    timeline = NV{k}.([
        (true, false, false, false, false, false) certain
        (true, false, false, true, false, false) uncertain
        (true, true, false, true, false, false) uncertain
        (false, true, true, false, true, false) uncertain
        (true, true, false, false, false, true) uncertain
    ])

    LandscapeChange.apply_transitions!(reverse(timeline), transitions, indirect) |> reverse
    LandscapeChange.apply_transitions!(copy(timeline), reversed, reversed_indirect)
    expected = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=false, cull=false)[1]
    expected == NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, true, false, false)
        (true, true, false, true, false, false)
        (true, true, true, true, true, false)
        (true, true, false, true, true, false)
    ])
    wanted = NV{k}.([
        (true, false, false, false, false, false)
        (true, false, false, false, false, false)
        (true, true, false, false, false, false)
        (false true, false, false, false, false)
        (false true, false, false, true, false)
    ])

        (true, false, false, true, false, false)
        (true, true, false, true, false, false)
end


nv_rasts.mus[Y=Contains(-20.202), X=Contains(57.500)] |> parent
compiled = map(nv_rasts) do nv_rast
    cross_validate_timeline(logic, nv_rast; simplify=true, cull=true)#, force)
end
compiled.mus.timeline[Y=Contains(-20.202), X=Contains(57.500)]
