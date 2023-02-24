"""
    cover_change(A, B; categories=union(A))

Calculate fraction of changes from each category to
all other categories, between arrays `B` and `B`.

Using a `NamedTuple` or `NamedVector` for `categories`
will give named results.
"""
function cover_change(A, B; categories=union(A))
    map(categories) do c1
        map(categories) do c2
            if c1 == c2
                NaN
            else
                n_c1_start = count(==(c1), A)
                n_c1_to_c2 = count(zip(A, B)) do (a, b)
                    a == c1 && b == c2
                end
                n_c1_to_c2 / n_c1_start
            end
        end
    end
end

"""
    cover_fraction(A; categories=union(A))

Calculate the fractions of the (non-missing valued)
map covered by each category.
"""
function cover_fraction(A; categories=union(A))
    total_count = count(_ -> true, skipmissing(A))
    map(categories) do c
        cover_count = count(==(c), skipmissing(A))
        cover_count / total_count
    end
end

"""
    cover_persistence(A, B; categories=union(A))

Calculate fraction of persistence of each category 
between arrays `B` and `B`.

Using a `NamedTuple` or `NamedVector` for `categories`
will give named results.
"""
function cover_persistence(A, B; categories=union(A))
    map(categories) do c
        nstart = count(==(c), A)
        npersisted = count(zip(A, B)) do (a, b)
            (a == c) & (b == c)
        end
        npersisted / nstart
    end
end
