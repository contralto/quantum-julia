#=
transform:
- Author: contralto
- Date: 2021-01-10
=#


@testset "Transform (loops vs matrices)" begin
    for n in 1:4
        state1 = init_state(Int(n))
        state2 = init_state(Int(n))

        for t in 0:n - 1
            transform!(state1, t, h)
            transform_with_matrix!(state2, t, h)
            @test state1 == state2
        end
    end
end


