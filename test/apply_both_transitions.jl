using LandscapeChange: apply_both_transitions

transitions = NV(
    native    = NV(native=true,  cleared=false, abandoned=false, urban=false, forestry=false,  water=false),
    cleared   = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    abandoned = NV(native=false, cleared=true,  abandoned=true,  urban=false,  forestry=false,  water=false),
    urban     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=false),
    forestry  = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=true,   water=false),
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=true,  water=true),
)
# from = to
reversed = LandscapeChange.reverse_transitions(transitions) 
indirect = LandscapeChange.indirect_transitions(transitions)
reversed_indirect = LandscapeChange.reverse_transitions(indirect) 
force = map(_ -> false, transitions)
logic = (; transitions, indirect, reversed, reversed_indirect)

# (:native, :cleared, :abandoned, :urban, :forestry, :water)

@testset "apply_both_transitions" begin

    @testset "Certainty wins" begin
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (true, true,  false, false, false, true),
        ))
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, true, false, false, false, true),
            (true, false,  false, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, true, false, false, false, true),
            (true, true,  false, true, true, false),
        ))
        result = apply_both_transitions(timeline, logic)
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
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
        ))
    end

    @testset "Multiple possible paths are combined while uncertainty is removed" begin
        timeline = NV{k}.((
            (true, false, true, false, true, true),
            (true, true,  false, true, false, true),
        ))
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, true, true, false, true, true),
            (true, true, true, false, true, true),
            (false, true,  true, true, false, true),
        ))
        result = apply_both_transitions(timeline, logic)
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
        result = apply_both_transitions(timeline, logic)
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
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (false, false, false, true, false, false),
            (false, false, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, true, false),
            (false, false, false, false, false, true),
        ))
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (false, false, false, false, true, false),
            (false, false, false, false, false, true),
        ))
    end

    @testset "Keep both certainties through impossible transitions" begin
        timeline = NV{k}.((
            (false, true, false, false, false, false),
            (true, false, false, false, false, false),
        ))
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (true, true, false, false, false, false),
            (true, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
        ))
        @test apply_both_transitions(timeline, logic) == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
        ))
    end

    @testset "Transitions take the most direct path through uncertainty" begin
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
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
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (false, false, true, true, false, false),
            (false, true,  false, false, false, true),
        ))
    end

    @testset "Multiple lines" begin
        timeline = NV{k}.((
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (false, true, false, false, false, true),
        ))
        timeline = NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, true, false),
            (false, true, true, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        @test result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
            (false, true, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, false, false, false, true, false),
            (false, false, false, false, false, true),
            (true, false, false, false, false, false),
            (false, true, true, false, false, false),
            (false, true, false, false, true, false),
            (false, true, true, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        collect(result)
        @test_broken result == NV{k}.((
            (true, false, false, false, true, false),
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (false, true, false, false, false, true),
            (false, true, false, false, false, true),
            (false, true, false, false, false, true),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, true,  false, false),
            (true,  false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        collect(result)
        @test_broken result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, true ),
            (true, false, false, true,  false, true ),
            (true, false, false, true,  false, true ),
            (true, false, false, true,  false, true ),
            (true, false, false, true,  false, true ),
        ))
        timeline = NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, true),
            (false, false, false, false, false, false),
            (false, false, true,  false, false, false),
            (true,  false, false, false, false, false),
        ))
        result = apply_both_transitions(timeline, logic)
        collect(result)
        @test_broken result == NV{k}.((
            (true, false, false, false, false, false),
            (true, false, false, false, false, false),
            (true, false, false, false, false, true),
            (true, false, false, false, false, true),
            (true, false, true,  false, false, true),
            (true, false, true,  false, false, true),
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
        result = apply_both_transitions(timeline, logic)
        collect(result)
        @test result == NV{k}.((
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
            (false, false, false, false, false, false),
        ))
        timeline = NV{k}.((
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
        ))
        result = apply_both_transitions(timeline, logic)
        collect(result)
        @test result == NV{k}.((
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
            (true, true, true, true, true, true),
        ))
    end
end
