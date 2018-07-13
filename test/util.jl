# 
# Run these tests by calling include("test/util.jl")
# 
using Base.Test, EarthMoversDistance
import Util

# 
# Test Util.mdpa
# 
NUM_LEVELS = 16 # size of histogram arrays
cityblock = (x, y) -> abs(x - y) # cityblock distance

# simple case 1)
# in each histogram, only one level is 1.0, the others are 0.0
for i in 1:NUM_LEVELS, j in 1:NUM_LEVELS
    
    # construct histograms
    histogram1 = zeros(Float64, NUM_LEVELS)
    histogram1[i] = 1.0
    histogram2 = zeros(Float64, NUM_LEVELS)
    histogram2[j] = 1.0
    
    # test cityblock distance
    mdpa = Util.mdpa(histogram1, histogram2)
    emd  = EarthMoversDistance.emd(histogram1, histogram2, cityblock)
    
    if mdpa != emd
        warn("Case 1 failing for cityblock distance with i = $i, j = $j")
    end
    @test mdpa == emd
    
end

# simple case 2)
# one histogram like before, but the other has two levels that are 0.5
for i in 1:NUM_LEVELS, j in 1:NUM_LEVELS, k in 1:NUM_LEVELS
    
    # construct histograms
    histogram1 = zeros(Float64, NUM_LEVELS)
    histogram1[i] = 1.0
    histogram2 = zeros(Float64, NUM_LEVELS)
    histogram2[j] = 0.5
    histogram2[k] = histogram2[k] + 0.5 # has to be 1 if j == k
    
    # test cityblock distance
    mdpa = Util.mdpa(histogram1, histogram2)
    emd  = EarthMoversDistance.emd(histogram1, histogram2, cityblock)
    
    if mdpa != emd
        warn("Case 2 failing with i = $i, j = $j, k = $k")
    end
    @test mdpa == emd
    
end

# case 3: random histograms with single flow)
for i in 1:100 # test hundred times
    histogram1 = rand(Float64, NUM_LEVELS)
    histogram2 = copy(histogram1)
    
    from = rand(1:NUM_LEVELS)
    to   = rand(setdiff(1:NUM_LEVELS, [from]))
    flow = rand(Float64) * min(histogram2[from], 1 - histogram2[to])
    histogram2[from] = histogram2[from] - flow
    histogram2[to]   = histogram2[to]   + flow
    
    mdpa = Util.mdpa(histogram1, histogram2)
    emd  = EarthMoversDistance.emd(histogram1, histogram2, cityblock)
    if !isapprox(mdpa, emd, atol=1e-6)
        println(histogram1)
        println("$from -> $to: $flow")
        println(histogram2)
        println("Case 3 (i = $i) failing with EMD $emd (MDPA $mdpa, difference $(abs(emd-mdpa)))")
    end
    @test isapprox(mdpa, emd, atol=1e-6)
end

# case 4: random histograms with 10 flows)
for i in 1:100 # test hundred times
    histogram1 = rand(Float64, NUM_LEVELS)
    histogram2 = copy(histogram1)
    
    for j in 1:10
        from = rand(1:NUM_LEVELS)
        to   = rand(setdiff(1:NUM_LEVELS, [from]))
        flow = rand(Float64) * min(histogram2[from], 1 - histogram2[to])
        histogram2[from] = histogram2[from] - flow
        histogram2[to]   = histogram2[to]   + flow
    end
    
    mdpa = Util.mdpa(histogram1, histogram2)
    emd  = EarthMoversDistance.emd(histogram1, histogram2, cityblock)
    if !isapprox(mdpa, emd, atol=1e-6)
        println("Case 4 failing with EMD $emd (upper bound $upperbound, difference $(emd-upperbound))")
    end
    @test isapprox(mdpa, emd, atol=1e-6)
end

info("All tests passed")

