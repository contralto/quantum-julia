#=
gates:
- Author: contralto
- Date: 2021-01-12
=#

@testset "Testing H" begin
    # hadamard across all targets creates equal amplitudes (tests 2 to 16 qubits)
    for val in 1:4
        state = init_state((val))
        for t in 0:val - 1
            transform!(state, Int(t), h)
        end
        state = round.(state, digits = 15)
        @test all(state .== (fill(round(1/sqrt(2)^val, digits = 15), 2^val)))
    end


    # test distance between nonzero amplitudes is correct
    # targeting qubit 0
    state = init_state(2)
    t0 = [1/sqrt(2), 0, 1/sqrt(2), 0]

    transform!(state, 0, h)
    @test all(state .== t0)



    # targeting qubit 1
    state = init_state(2)
    t1 = [1/sqrt(2), 1/sqrt(2), 0, 0]

    transform!(state, 1, h)
    @test all(state .== t1)

end

@testset "Testing X" begin
    # checks that x swapped the first amplitude for each target
    n = 4
    state = init_state(n)
    for t in 0:(n - 1)
        expected = fill(0, 2^n)
        changed = 2 ^ ((n - 1) - t) + 1
        expected[changed] = 1

        transform!(state, t, x)
        @test all(state .== expected)
        transform!(state, t, x)
    end
end

@testset "Testing Z" begin
    n = 5
    for t in 0:(n - 1)
        state = init_state(n)
        transform!(state, t, h)

        expected = fill(0.0, 2^n)
        changed = 2 ^ (n - 1 - t) + 1
        expected[1] = 1/sqrt(2.0)
        expected[changed] = -1/sqrt(2)

        transform!(state, t, z)
        @test all(state .== expected)
    end
end