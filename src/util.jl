macro swap!(p,q)
    return quote
        begin
            tmp = $(esc(p))
            $(esc(p)) = $(esc(q))
            $(esc(q)) = tmp
        end
    end
end
