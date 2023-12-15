"""
    fuzzy_scores(predicted, real; neighborhood, ignored=())

Return a fuzzy score for each cell. Compares each cell of the two arrays.
It finds the shortest distance within the neighborhood where the
value in the `real` array is present in `predict` array.

The longer the distance, the lower score for this cell. Maximum
is 1 and if not present its 0. The weight function is linear.

- `predict`: Array of predicted data
- `real`: Array of real target data

# Keywords

- `neighborhood`: any `Neighborhood` oject, `Window(1)` by default.
- `banlist`: Values in the arrays which will not be evaluated.

Adapted from: https://github.com/johanlahti/urban-lu-model/blob/master/script/Utils.py
"""
function fuzzy_scores(predicted::AbstractDimArray, real::AbstractDimArray; kw...)
    fuzzy = fuzzy_scores(parent(predicted), parent(real); kw...)
    return rebuild(predicted; data=fuzzy, name=:fuzzy)
end
function fuzzy_scores(predicted, real; neighborhood=Moore{2}(), ignored=0)
    broadcast_neighborhood(neighborhood, real, predicted) do h, r, p
        r === p && return _weight(0, h)
        r in ignored && return 0.0

        weight = 0.0
        shortestdist = Inf

        # Search the neighborhood
        for (neighbor, dist) in zip(neighbors(h), distances(h))
            if p === neighbor && dist < shortestdist
                # Store the shortest distance at which we found the value
                shortestdist = dist
                weight = _weight(dist, h)
            end
            shortestdist == 1 && break # shortcut when there are no shorter distances
        end
        return weight
    end
end

_weight(dist, h) = 1 - dist / (length(h) + 1)
