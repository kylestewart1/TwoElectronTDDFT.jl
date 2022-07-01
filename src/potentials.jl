function soft_coulomb(x; Z=1.0)
    @. Z / sqrt(x^2 + 1.0)
end

function harmonic(x; k=1.0)
    @. 0.5 * k * x^2
end

