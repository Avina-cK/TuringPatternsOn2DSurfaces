function fts(x, n)
    # Small values: format as *-**e-*
    if x != 0 && abs(x) < 1.0e-4
        e = floor(Int, log10(abs(x)))
        m = round(x / 10.0^e, digits=n)
        if abs(m) >= 10 # rounding pushed mantissa to 10.0
            m /= 10
            e += 1
        end
        s = string(m)
        i = findfirst('.', s)
        s = rpad(s, i + n, '0')
        return replace(s, "." => "-") * "e" * string(e)
    end

    # Old version (if breaks)
    s = string(round(x, digits=n))
    i = findfirst('.', s)
    s = rpad(s, i + n, '0')
    return replace(s, "." => "-")
end
