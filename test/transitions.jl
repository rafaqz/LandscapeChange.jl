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
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=true,  water=true),
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
