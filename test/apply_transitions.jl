using LandscapeChange, Test
using LandscapeChange: apply_transitions, Forced, Uncertain, all_transitions

const NV = NamedVector

direct_transitions = NV(
    native    = NV(native=true,  cleared=false, abandoned=false, urban=false, forestry=false,  water=false),
    cleared   = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    abandoned = NV(native=false, cleared=true,  abandoned=true,  urban=false,  forestry=false,  water=false),
    urban     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=false),
    forestry  = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=true,   water=false),
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=true,  water=true),
)
transitions = all_transitions(direct_transitions)
reversed = all_transitions(LandscapeChange.reverse_transitions(direct_transitions))
logic = (; transitions, reversed)
k = (:native, :cleared, :abandoned, :urban, :forestry, :water)

@testset "merge forced" begin
    a = NV{k}((true, false, false, false, false, false))
    b = NV{k}((true, true,  false, false, false, true))
    # Dest is totally forced
    @test LandscapeChange._merge(Forced(), Forced(), a, b, reversed) ==
        NV{k}((true, true,  false, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), a, b, transitions) ==
        NV{k}((true, true, false, false, false, true))
    # Source is minimised where transitions are possible
    @test LandscapeChange._merge(Forced(), Forced(), b, a, reversed) ==
        NV{k}((true, false,  false, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), b, a, transitions) ==
        NV{k}((true, false, false, false, false, false))

    a = NV{k}((false, false, true, false, false, true))
    b = NV{k}((false, true,  false, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), a, b, reversed) ==
        NV{k}((false, true,  false, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), b, a, reversed) ==
        NV{k}((false, false, true, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), a, b, transitions) ==
        NV{k}((false, true,  false, false, false, true))
    @test LandscapeChange._merge(Forced(), Forced(), b, a, transitions) ==
        NV{k}((false, false, true, false, false, true))
end

@testset "merge uncertain" begin
    a = NV{k}((true, false, false, true, false, false))
    b = NV{k}((true, true,  false, false, false, true))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), a, b, reversed) ==
        NV{k}((true, false,  false, false, false, false))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), b, a, reversed) ==
        NV{k}((true, false,  false, false, false, false))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), a, b, transitions) ==
        NV{k}((true, false,  false, false, false, false))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), b, a, transitions) ==
        NV{k}((true, false,  false, false, false, false))

    # a = NV{k}((true, true, false, true, false, false))
    # b = NV{k}((false, false, true, false, false, false))
    # @test LandscapeChange._merge(Uncertain(), Uncertain(), a, b, reversed) ==
    #     NV{k}((false, false,  true, true, false, false))
    # @test LandscapeChange._merge(Uncertain(), Uncertain(), b, a, reversed) ==
    #     NV{k}((false, false,  true, true, false, false))
    # LandscapeChange._merge(Uncertain(), Uncertain(), a, b, transitions)
    # LandscapeChange._merge(Uncertain(), Uncertain(), b, a, transitions)

    a = NV{k}((false, true, false, true, false, true))
    b = NV{k}((true, false, true, false, false, false))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), a, b, reversed) ==
        NV{k}((false, false,  true, false, false, false))
    @test LandscapeChange._merge(Uncertain(), Uncertain(), b, a, reversed) ==
        NV{k}((false, true,  false, true, false, true))
    # @test LandscapeChange._merge(Uncertain(), Uncertain(), a, b, transitions, (reversed,)) ==
    #     NV{k}((true, false,  true, false, false, false))
    # @test LandscapeChange._merge(Uncertain(), Uncertain(), b, a, transitions, (reversed,)) ==
    #     NV{k}((false, true,  false, true, false, true))
end

@testset "apply_transitions" begin
    @testset "Certainty wins" begin
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (true, true,  false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, true, false, false, false, true),
            (true, false,  false, false, false, false),
        ))
        result = apply_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, true, false, false, false, true),
            (true, true,  false, true, true, false),
        ))
        result = apply_transitions(timeline, logic)
        @test result == NV{k}.((
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
    end

    @testset "Uncertainty simplifies to certainty" begin
        timeline = NV{k}.((
            (true, false, true, false, false, false),
            (true, true,  false, false, false, true),
        ))
        @test apply_transitions(timeline, logic) == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
    end

    @testset "Multiple possible paths are combined while uncertainty is removed" begin
        timeline = NV{k}.((
            (true, false, true, false, true, true),
            (true, true,  false, true, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, true, true, false, true, true),
            (true, true, true, false, true, true),
            (false, true,  true, true, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, true, true, false, false, true),
            (false, true, true, false, false, true),
            (false, true,  true, false, false, true),
        ))
        timeline = NV{k}.((
            (true, true, true, false, true, true),
            (true, true, true, false, true, true),
            (false, true,  true, true, false, true),
            (false, true,  true, true, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, true, true, false, false, true),
            (false, true, true, false, false, true),
            (false, true,  true, false, false, true),
            (false, true,  true, false, false, true),
        ))
    end

    @testset "Allow transitions in known values" begin
        timeline = NV{k}.((
            (false, false, false, true, false, false),
            (false, false, false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, false, true, false, false),
            (false, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, true, false),
            (false, false, false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, false, false, true, false),
            (false, false, false, false, false, true),
        ))
    end

    @testset "Keep both certainties through impossible transitions" begin
        timeline = NV{k}.((
            (false, true, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        results = apply_transitions(timeline, logic); collect(results)
        @test results == NV{k}.((
            (true, true, false, false, false, false),
            (true, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
        ))
        results = apply_transitions(timeline, logic); collect(results)
        @test results == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
        ))
    end
    a = NV{k}((false, false, false, false, false, true))
    b = NV{k}((true, false, false, true, false, false))
    @testset "Transitions take the most direct path through uncertainty" begin
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
    end

    @testset "Unresolvable uncertainty is kept" begin
        timeline = NV{k}.((
            (false, false, true, true, false, false),
            (false, true,  false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, true, true, false, false),
            (false, true,  false, false, false, true),
        ))
    end

    @testset "All the same" begin
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
            (true, true, true, true, true, false),
        ))
    end

    @testset "Multiple lines working" begin
        timeline = NV{k}.((
            (true, false, false, false, false, true),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, true, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (true,  false, false, false, false, false),
            (false, false, false, false, true, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (false, false, false, false, true, false),
        ))
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, true, false),
            (false, true, true, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, false, false, false),
            (false, false, true,  false, false, false),
            (true,  false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (true, false, true , false, false, true),
            (true, false, true , false, false, true),
            (true, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (false, true, false, false, false, true),
            (false, true, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, true, true, true, true, true),
            (false, true, true, true, true, false),
            (false, false, true, true, false, false),
            (false, false, false, true, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, false, true, false, false),
            (false, false, false, true, false, false),
            (false, false, false, true, false, false),
            (false, false, false, true, false, false),
        ))
        timeline = NV{k}.((
            (true, true, true, true, true, false),
            (false, true, true, true, true, false),
            (false, false, true, true, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, true, true, false, false),
            (false, false, true, true, false, false),
            (false, false, true, true, false, false),
            (false, false, true, true, false, false),
        ))
        timeline = NV{k}.((
            (false, false, true, true, true, true),
            (false, false, false, true, true, true),
            (false, false, false, false, true, true),
            (false, false, false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, false, false, false, false, true),
            (false, false, false, false, false, true),
            (false, false, false, false, false, true),
            (false, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, true, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (false, true, false, false, false, true),
            (false, true, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, true, true, true, false, false),
            (true, true, true, false, false, false),
            (true, true, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
    end

    @testset "No unnecessary forcing with multiple forced" begin
        timeline = NV{k}.((
            (true,  false, false, false, false, false),
            (true,  false, false, false, false, false),
            (true,  false, false, false, false, false),
            (true,  false, false, false, false, false),
            (true,  false, false, false, false, false),
            (true,  false, false, false, false, false),
            (true,  false, true , false, false, false),
            (true,  false, false, false, false, false),
            (false, true , false, false, false, false),
            (false, false, false, false, false, false),
            (false, true , false, true , false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, true , false, true , false, false),
            (false, false, false, false, false, false),
            (false, true , false, false, false, false),
            (false, false, true , true , false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)

        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, true,  false, false),
            (true,  false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (true, false, false, true, false, true),
            (true, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (false, false, false, false, false, false),
            (false, false, true,  false, false, false),
            (true,  false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, true , false, false, true),
            (true, false, true , false, false, true),
            (true, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, true,  false, false),
            (true,  false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true ),
            (true, false, false, true,  false, true ),
            (true, false, false, false,  false, true ),
            (true, false, false, false,  false, true ),
            (true, false, false, false,  false, true ),
        ))

        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, true,  false, false),
            (false, true,  false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (false, true,  false, false, false, true),
            (false, true,  false, false, false, true ),
            (false, true,  false, true,  false, true ),
            (false, true,  false, false, false, true ),
            (false, true, false, false, false, true ),
            (false, true, false, false, false, true ),
        ))

        timeline = NV{k}.((
            (false, false, false, true, false,  false),
            (false, false, false, false, true,  true),
            (true,  false, false, false, false, false),
            (false, false, true,  false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, true,  false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true,  false, false, true,  false, false),
            (true,  false, false, false, false, true ),
            (true,  false, false, false, false, true ),
            (false, false, true,  false, false, true ),
            (false, false, true,  false, true,  true ),
            (false, false, false, false, true,  true ),
        ))

        # timeline = NV{k}.((
        #     (false, false, false, true,  false,  false),
        #     (true,  false, true,  false, true,  true),
        #     (true,  false, false, false, false, false),
        #     (false, false, true,  false, false, false),
        #     (false, false, false, false, false, false),
        #     (false, false, false, false, true,  false),
        # ))
        # result = apply_transitions(timeline, logic); collect(result)

        timeline = NV{k}.((
            (false, true , false, false, false, false),
            (false, false, true , false, false, false),
            (true,  false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, true ),
            (true,  false, true , true , false, true ),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true,  true , false, false, false, false),
            (true,  false, true , false, false, false),
            (true,  false, true , false, false, false),
            (true,  false, true , false, false, true ),
            (false, false, false, false, false, true ),
            (false, false, false, false, false, true ),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, true , false, true , false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, true,  false, false, false, true ),
            (false, false, false, false, false, true ),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
             (false, false, true,  false, true , false),
             (false, false, true,  false, true , false),
             (false, false, true,  false, true , true ),
             (false, false, true , false, true , true ),
             (false, false, false, false, false, true ),
             (false, false, false, false, false, true ),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true ),
            (true , false, false, true , true , false),
            (false, true , false, true,  false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test_broken result == NV{k}.((
            (false, false, false, false, false, true),
            (false, false, false, true , false, true),
            (false, false, false, true , false, true),
            (false, false, false, false, false, true),
            (false, false, false, false, false, true),
            (false, false, false, false, false, true),
        ))
        # timeline = NV{k}.((
        #     (true , false, false, true , true , false),
        #     (false, true , false, true,  false, false),
        #     (false, false, false, false, false, false),
        #     (false, false, false, false, false, true ),
        #     (false, false, false, false, false, false),
        #     (false, false, false, false, false, false),
        # ))
        # result = apply_transitions(timeline, logic); collect(result)

        timeline = NV{k}.((
            (true , false, false, false, false, false),
            (true , true , false, false, false, false),
            (true , true , false, false, false, false),
            (false, true , false, false, false, false),
            (true , true , true, false, false, false),
            (false , false, true, true, false, false),
            (false, false, true, false, false, false),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, true, false, false, false, false),
            (true, true, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, true, false, false, false),
            (false, false, true, false, false, false),
            (false, false, true, false, false, false),
        ))

        timeline = NV{k}.((
            (true , false, false, false, false, false),
            (true , true , false, false, false, false),
            (true , true , false, false, false, false),
            (false, true , false, false, false, false),
            (true , true , true, false, false, false),
            (false , false, true, true, false, false),
            (false, false, true, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
        ))
        result = apply_transitions(timeline, logic); collect(result)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, true, false, false, false, false),
            (true, true, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, true, false, false, false),
            (false, false, true, false, false, false),
            (false, false, true, false, false, false),
            (false, false, true, false, false, true),
            (false, false, false, false, false, true),
        ))
    end
end

nothing
